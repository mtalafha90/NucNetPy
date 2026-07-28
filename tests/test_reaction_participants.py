"""Regression tests for issues found by running real JINA/libnucnet data.

Each test here pins behaviour that was previously wrong and that silently
produced a plausible-looking but incorrect answer, which is the failure mode
that matters most for a scientific library.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from nucnetpy import (Network, RateFit, Reaction, Species, Zone, analytic_jacobian,
                      evolve_zone, constant_thermo, jacobian, solve_nse, time_grid)
from nucnetpy.solver import rhs


def _photodisintegration_network():
    """o16 + gamma -> he4 + c12 with a constant rate, plus its inverse."""
    net = Network()
    net.add_species(Species("he4", 2, 4, mass_excess=2.425, spin=0.0))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0, spin=0.0))
    net.add_species(Species("o16", 8, 16, mass_excess=-4.737, spin=0.0))
    net.reactions.add(Reaction.from_names(
        ["o16", "gamma"], ["he4", "c12"], constant_rate=3.0, q_value=-7.162))
    return net


def test_photon_does_not_enter_the_reactant_product():
    """A photodisintegration must flow even though Y(gamma) is never stored.

    The photon appears in the reaction record so that charge and mass number
    balance, but it is part of the thermodynamic state rather than the
    abundance vector.  Counting it as a reactant multiplied every
    photodisintegration flow by zero.
    """
    net = _photodisintegration_network()
    reaction = net.reactions.reactions[0]

    assert [p.species for p in reaction.nuclear_reactants] == ["o16"]
    assert reaction.reactant_order == 1

    abundances = {"o16": 0.05}
    flux = reaction.flux(abundances, t9=2.0, rho=1.0e6)
    assert flux > 0.0
    # One-body process: the flow must not pick up a density factor.
    assert flux == pytest.approx(3.0 * 0.05)
    assert reaction.flux(abundances, t9=2.0, rho=1.0e9) == pytest.approx(flux)


def test_photon_is_not_an_evolved_species():
    net = _photodisintegration_network()
    zone = Zone(abundances={"o16": 1.0 / 16.0})
    net.add_zone(zone)
    result = evolve_zone(net, zone, time_grid(0.0, 1.0e-2, 8),
                         thermo=constant_thermo(t9=2.0, rho=1.0e6), method="bdf")
    assert result.success
    assert "gamma" not in result.species
    # Oxygen is consumed and the alpha/carbon pair produced.
    final = result.final_abundances
    assert final["o16"] < 1.0 / 16.0
    assert final["he4"] > 0.0
    assert final["c12"] > 0.0
    # Baryon number is conserved by the flow.
    xsum = sum(net.species[k].a * v for k, v in final.items())
    assert xsum == pytest.approx(1.0, abs=1e-10)


def test_statistical_factor_ignores_massless_participants():
    net = Network()
    for name, z, a in [("he4", 2, 4), ("c12", 6, 12)]:
        net.add_species(Species(name, z, a))
    triple_alpha = Reaction.from_names(["he4", "he4", "he4"], ["c12", "gamma"])
    # 3! for the three identical alphas; the photon must not contribute.
    assert triple_alpha.statistical_factor() == 6
    assert triple_alpha.reactant_order == 3


def test_beta_decay_conserves_charge_without_leptons_in_the_species_map():
    """Charge balance must not depend on leptons being registered as species."""
    net = Network()
    net.add_species(Species("n", 0, 1))
    net.add_species(Species("h1", 1, 1))
    decay = Reaction.from_names(["n"], ["h1", "electron", "anti-neutrino_e"])
    ok, da, dz = decay.conserves_a_z(net.species)
    assert (da, dz) == (0, 0)
    assert ok


@pytest.mark.parametrize("method", ["bdf", "radau", "lsoda"])
def test_stiff_solvers_are_dispatched_to_scipy(method):
    """Radau was previously passed to SciPy as 'RADAU' and silently rejected."""
    net = _photodisintegration_network()
    zone = Zone(abundances={"o16": 1.0 / 16.0})
    net.add_zone(zone)
    result = evolve_zone(net, zone, time_grid(0.0, 1.0e-2, 8),
                         thermo=constant_thermo(t9=2.0, rho=1.0e6), method=method)
    assert result.success
    assert "fallback" not in result.message


def test_analytic_jacobian_matches_finite_differences():
    net = Network()
    for name, z, a in [("he4", 2, 4), ("c12", 6, 12), ("o16", 8, 16)]:
        net.add_species(Species(name, z, a))
    net.reactions.add(Reaction.from_names(
        ["he4", "he4", "he4"], ["c12"], rate_fits=[RateFit([10, 0, 0, 0, 0, 0, 0])]))
    net.reactions.add(Reaction.from_names(
        ["c12", "he4"], ["o16"], rate_fits=[RateFit([8, 0, 0, 0, 0, 0, 0])]))
    species = ["he4", "c12", "o16"]
    thermo = constant_thermo(t9=2.0, rho=1.0e5)

    ja = analytic_jacobian(net, species, thermo, sparse=False)
    jf = jacobian(net, species, thermo, sparse=False)
    for y in (np.array([0.25, 0.0, 0.0]), np.array([0.1, 0.02, 0.005])):
        A = np.asarray(ja(0.0, y))
        F = np.asarray(jf(0.0, y))
        scale = max(float(np.abs(F).max()), 1e-30)
        assert np.abs(A - F).max() / scale < 1e-5


def test_analytic_jacobian_is_correct_when_an_abundance_vanishes():
    """d/dY of a linear reactant stays finite at Y = 0."""
    net = _photodisintegration_network()
    species = ["he4", "c12", "o16"]
    thermo = constant_thermo(t9=2.0, rho=1.0e6)
    J = np.asarray(analytic_jacobian(net, species, thermo, sparse=False)(
        0.0, np.array([0.0, 0.0, 0.05])))
    # d(dY_he4/dt)/dY_o16 = +3, d(dY_o16/dt)/dY_o16 = -3.
    assert J[species.index("he4"), species.index("o16")] == pytest.approx(3.0)
    assert J[species.index("o16"), species.index("o16")] == pytest.approx(-3.0)


def test_nse_uses_the_ground_state_spin_degeneracy():
    """libnucnet partition tables are normalised, so (2J+1) must be applied.

    Two otherwise identical nuclides differing only in ground-state spin must
    appear in the ratio of their statistical weights.
    """
    net = Network()
    net.add_species(Species("n", 0, 1, mass_excess=8.071, spin=0.5))
    net.add_species(Species("h1", 1, 1, mass_excess=7.289, spin=0.5))
    # Two fictitious mirror nuclides with identical masses but different spins.
    net.add_species(Species("he4", 2, 4, mass_excess=2.425, spin=0.0))
    net.add_species(Species("li7", 3, 7, mass_excess=14.907, spin=1.5))

    from nucnetpy.nse import _log_prefactor
    spin_zero = Species("c12", 6, 12, mass_excess=-5.0, spin=0.0)
    spin_two = Species("c12", 6, 12, mass_excess=-5.0, spin=2.0)
    delta = _log_prefactor(spin_two, 5.0, 1e7) - _log_prefactor(spin_zero, 5.0, 1e7)
    assert delta == pytest.approx(math.log(5.0))

    # A species with unknown spin must not be given an invented weight.
    unknown = Species("c12", 6, 12, mass_excess=-5.0, spin=None)
    assert _log_prefactor(unknown, 5.0, 1e7) == pytest.approx(
        _log_prefactor(spin_zero, 5.0, 1e7))


def test_nse_excludes_species_without_nuclear_data():
    """A placeholder mass excess of zero must not enter an equilibrium solve.

    Species referenced by a reaction file but absent from the nuclide file are
    synthesised with a zero mass excess.  Left in, an unbound nuclide such as
    li5 competes with the iron peak.
    """
    net = Network()
    net.add_species(Species("n", 0, 1, mass_excess=8.071, spin=0.5))
    net.add_species(Species("h1", 1, 1, mass_excess=7.289, spin=0.5))
    net.add_species(Species("he4", 2, 4, mass_excess=2.425, spin=0.0))
    net.add_species(Species("ni56", 28, 56, mass_excess=-53.907, spin=0.0))
    net.add_species(Species("li5", 3, 5, mass_excess=0.0))
    net.metadata["species_without_nuclear_data"] = "li5"

    assert net.species_without_nuclear_data() == ["li5"]
    assert net.validate()["species_without_nuclear_data"] == ["li5"]

    guarded = solve_nse(net, t9=5.0, rho=1.0e7, ye=0.5)
    assert guarded.success
    assert "li5" not in guarded.abundances
    assert guarded.xsum == pytest.approx(1.0, rel=1e-6)

    # Opting in restores the old behaviour, and shows why it is not the default.
    unguarded = solve_nse(net, t9=5.0, rho=1.0e7, ye=0.5, require_nuclear_data=False)
    assert "li5" in unguarded.abundances
    assert 5 * unguarded.abundances["li5"] > 1.0e-3


def test_failed_stiff_solve_is_reported_as_a_failure():
    """A solver that cannot start must not be reported as a success.

    The previous behaviour substituted an explicit fixed-step pass and returned
    success=True, so a diverged stiff integration looked like a valid result.
    """
    net = Network()
    net.add_species(Species("he4", 2, 4))
    net.add_species(Species("c12", 6, 12))
    # An enormous constant rate makes the requested step unreachable.
    net.reactions.add(Reaction.from_names(
        ["he4", "he4", "he4"], ["c12"], constant_rate=1.0e300))
    zone = Zone(abundances={"he4": 0.25})
    net.add_zone(zone)
    result = evolve_zone(net, zone, time_grid(0.0, 1.0e6, 5),
                         thermo=constant_thermo(t9=1.0, rho=1.0e8), method="bdf")
    assert not result.success


def test_flux_saturates_instead_of_raising_on_overflow():
    """A diverging explicit step must not crash the right-hand side."""
    net = Network()
    net.add_species(Species("he4", 2, 4))
    net.add_species(Species("c12", 6, 12))
    reaction = Reaction.from_names(["he4", "he4", "he4"], ["c12"], constant_rate=1.0)
    value = reaction.flux({"he4": 1.0e200}, t9=1.0, rho=1.0e10)
    assert math.isinf(value)
