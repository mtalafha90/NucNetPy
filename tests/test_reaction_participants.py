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


def test_positivity_projection_that_creates_mass_is_reported():
    """Clipping a genuinely negative abundance invents baryon number.

    The right-hand side clips its input to non-negative values, so a component
    that goes negative stops influencing the derivative and its Jacobian column
    vanishes; an implicit solver can then carry it far negative while still
    satisfying its own convergence test.  Projecting that away silently adds
    mass, so the amount added must be measured and reported.
    """
    from nucnetpy.solver import _project_positive

    net = Network()
    net.add_species(Species("he4", 2, 4))
    net.add_species(Species("c12", 6, 12))
    species = ["he4", "c12"]

    # A trajectory whose second row carries a large negative abundance.
    y = np.array([[0.25, 0.0], [-0.25, 0.10]])
    clipped, success, message = _project_positive(net, species, y, True, "ok")
    assert (clipped >= 0.0).all()
    assert not success
    assert "positivity projection created" in message

    # Round-off level negatives must not be reported as a failure.
    y_small = np.array([[0.25, 0.0], [-1.0e-18, 0.25]])
    _, success_small, message_small = _project_positive(
        net, species, y_small, True, "ok")
    assert success_small
    assert message_small == "ok"


def test_flux_saturates_instead_of_raising_on_overflow():
    """A diverging explicit step must not crash the right-hand side."""
    net = Network()
    net.add_species(Species("he4", 2, 4))
    net.add_species(Species("c12", 6, 12))
    reaction = Reaction.from_names(["he4", "he4", "he4"], ["c12"], constant_rate=1.0)
    value = reaction.flux({"he4": 1.0e200}, t9=1.0, rho=1.0e10)
    assert math.isinf(value)


def test_detailed_balance_network_relaxes_to_nse():
    """A network with detailed-balance reverse rates must reproduce NSE.

    Rate libraries fit the forward and reverse directions independently, so
    their ratio is not exactly the equilibrium constant implied by the masses
    and partition functions.  A network built from them therefore relaxes to a
    stationary state that is not the NSE of the same nuclear data.  Rebuilding
    the reverse rates from detailed balance removes the inconsistency.
    """
    from nucnetpy import consistent_reverse_network

    net = Network()
    # Two nuclides and a radiative capture, with a deliberately inconsistent
    # library reverse rate (a factor of three away from detailed balance).
    net.add_species(Species("he4", 2, 4, mass_excess=2.425, spin=0.0))
    net.add_species(Species("c12", 6, 12, mass_excess=0.0, spin=0.0))
    net.add_species(Species("o16", 8, 16, mass_excess=-4.737, spin=0.0))
    forward = Reaction.from_names(["c12", "he4"], ["o16", "gamma"],
                                  constant_rate=1.0e3, q_value=7.162)
    net.reactions.add(forward)
    net.reactions.add(Reaction.from_names(["o16", "gamma"], ["c12", "he4"],
                                          constant_rate=3.0, q_value=-7.162))

    rebuilt = consistent_reverse_network(net)
    assert len(rebuilt.reactions.reactions) == 2
    reverse = [r for r in rebuilt.reactions.reactions
               if r.source == "detailed_balance"]
    assert len(reverse) == 1
    reverse = reverse[0]
    # The exothermic direction is kept as the forward reaction.
    assert {p.species for p in reverse.reactants} == {"o16", "gamma"}
    # Its rate is recomputed rather than copied from the library.
    assert reverse.rate_function is not None
    assert reverse.bare_rate(3.0) != pytest.approx(3.0)

    # The rebuilt reverse rate must satisfy detailed balance by construction.
    from nucnetpy.detailed_balance import reverse_rate
    for t9 in (1.0, 3.0, 5.0):
        assert reverse.bare_rate(t9) == pytest.approx(
            reverse_rate(forward, net.species, t9), rel=1e-12)


def test_reaction_rate_function_is_evaluated():
    reaction = Reaction.from_names(["he4"], ["c12"], rate_function=lambda t9: 2.0 * t9)
    assert reaction.bare_rate(3.0) == pytest.approx(6.0)
    # A callable composes with the other rate representations.
    reaction.constant_rate = 1.0
    assert reaction.bare_rate(3.0) == pytest.approx(7.0)


def test_skynet_screening_recovers_the_weak_limit():
    """At low density mu(Z) must scale as Z^2 and match Salpeter screening."""
    from nucnetpy.screening import SkyNetScreening, weak_screening_factor

    species = {"he4": Species("he4", 2, 4), "c12": Species("c12", 6, 12),
               "o16": Species("o16", 8, 16)}
    screening = SkyNetScreening(species)
    composition = {"he4": 0.25}

    screening.update(composition, t9=1.0, rho=1.0)
    mu2 = screening.chemical_potential(2)
    assert mu2 < 0.0
    # Weak screening is quadratic in charge.
    assert screening.chemical_potential(6) / mu2 == pytest.approx(9.0, rel=1e-6)
    assert screening.chemical_potential(8) / mu2 == pytest.approx(16.0, rel=1e-6)
    assert screening.chemical_potential(0) == 0.0

    # The enhancement of a two-body reaction must agree with the independent
    # Salpeter pairwise formula where weak screening applies.
    reaction = Reaction.from_names(["c12", "he4"], ["o16"])
    ion = sum(species[k].z ** 2 * v for k, v in composition.items())
    for rho in (1.0, 1.0e2, 1.0e4):
        screening.update(composition, t9=1.0, rho=rho)
        assert screening.factor(reaction) == pytest.approx(
            weak_screening_factor(6, 2, 1.0, rho, 0.5, ion), rel=2e-3)


def test_skynet_screening_grows_with_coupling_and_leaves_neutrals_alone():
    from nucnetpy.screening import SkyNetScreening

    species = {"he4": Species("he4", 2, 4), "n": Species("n", 0, 1),
               "si28": Species("si28", 14, 28)}
    screening = SkyNetScreening(species)
    composition = {"si28": 1.0 / 28.0}

    previous = 1.0
    charged = Reaction.from_names(["si28", "he4"], ["n"])
    for t9, rho in [(5.0, 1.0e6), (5.0, 1.0e8), (1.0, 1.0e8), (0.5, 1.0e9)]:
        screening.update(composition, t9=t9, rho=rho)
        factor = screening.factor(charged)
        assert factor >= previous
        previous = factor

    # A reaction with no charged reactant is unscreened.
    neutral = Reaction.from_names(["n", "n"], ["he4"])
    assert screening.factor(neutral) == 1.0


def test_screening_model_is_refreshed_with_the_composition():
    """The network must hand a composition-dependent model the current state."""
    from nucnetpy.reactions import ReactionNetwork

    seen = []

    class Recording:
        def update(self, abundances, t9, rho, ye):
            seen.append((dict(abundances), t9, rho))

        def __call__(self, reaction, t9=0.0, rho=0.0, ye=None):
            return 2.0

    net = ReactionNetwork()
    net.add(Reaction.from_names(["he4"], ["c12"], constant_rate=1.0))
    out = net.ydot({"he4": 0.25}, t9=2.0, rho=1.0e5, screening=Recording())
    # update() is called once for the evaluation, not once per reaction.
    assert len(seen) == 1
    assert seen[0][1] == 2.0 and seen[0][2] == 1.0e5
    # and the screening factor was applied to the flow
    assert out["he4"] == pytest.approx(-2.0 * 0.25)


def test_nuclear_energy_release_from_mass_excesses():
    """Energy release must follow from the change in total mass excess."""
    from nucnetpy.analysis import (nuclear_energy_generation_rate,
                                   nuclear_energy_release)
    from nucnetpy.constants import AVOGADRO, MEV_TO_ERG

    net = Network()
    net.add_species(Species("si28", 14, 28, mass_excess=-21.493, spin=0.0))
    net.add_species(Species("ni56", 28, 56, mass_excess=-53.907, spin=0.0))
    net.add_zone(Zone(abundances={"si28": 1.0 / 28.0}))

    # Burning all of the silicon to nickel releases 2*ME(Si) - ME(Ni) per
    # nickel nucleus formed.
    released = nuclear_energy_release(net, {"si28": 1.0 / 28.0}, {"ni56": 1.0 / 56.0})
    q = 2.0 * (-21.493) - (-53.907)
    assert released == pytest.approx(q / 56.0 * AVOGADRO * MEV_TO_ERG, rel=1e-12)
    assert released > 0.0

    # The rate is the same quantity per unit time, so a supplied derivative
    # integrates to the release above.
    rate = nuclear_energy_generation_rate(
        net, abundances={"si28": 1.0 / 28.0},
        ydot_values={"si28": -1.0 / 28.0, "ni56": 1.0 / 56.0}, t9=5.0, rho=1e8)
    assert rate == pytest.approx(released, rel=1e-12)
