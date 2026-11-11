"""ReacLib electron captures must be multiplied by the electron density.

ReacLib fits its three electron captures (on 3He and 7Be, and the pep
reaction p + p + e-) per unit rho * Ye (Cyburt et al. 2010, ApJS 189, 240).
NucNetPy treated the electron like the photon, whose density a rate already
contains, and so ran these captures too slowly by the factor rho * Ye.
"""
import math

import numpy as np
import pytest

from nucnetpy.analysis import charge_changing_flows
from nucnetpy.core import Network, Zone
from nucnetpy.reactions import RateFit, Reaction, ReactionNetwork, electron_fraction
from nucnetpy.solver import constant_thermo, evolve_zone, rhs
from nucnetpy.species import Species

# 7Be(e-,nu)7Li from ReacLib (label "bkmo"), per unit rho * Ye.
BE7_FIT = [-23.8328, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5]
# p(p e-,nu)d, the pep reaction (label "ec"), per unit rho * Ye.
PEP_FIT = [-43.6499, -0.00246064, -2.7507, -0.424877, 0.015987, -0.000690875, -0.207625]


def _be7_capture():
    return Reaction.from_names(["be7", "electron"], ["li7", "neutrino_e"],
                               rate_fits=[RateFit(BE7_FIT)])


def test_capture_rate_carries_the_electron_density():
    reaction = _be7_capture()
    fit = RateFit(BE7_FIT).rate(0.015)
    assert reaction.electron_count == 1
    assert reaction.rate(0.015, rho=150.0, ye=0.86) == pytest.approx(fit * 150.0 * 0.86,
                                                                     rel=1e-14, abs=0)
    # Without Ye the fit is returned per unit electron density.
    assert reaction.rate(0.015, rho=150.0) == pytest.approx(fit, rel=1e-14, abs=0)


def test_pep_flow_has_two_powers_of_density_and_the_identical_particle_factor():
    pep = Reaction.from_names(["h1", "h1", "electron"], ["h2", "neutrino_e"],
                              rate_fits=[RateFit(PEP_FIT)])
    t9, rho, ye, yp = 0.015, 150.0, 0.86, 0.7
    expected = RateFit(PEP_FIT).rate(t9) * rho * (rho * ye) * yp ** 2 / 2.0
    assert pep.flux({"h1": yp}, t9, rho=rho, ye=ye) == pytest.approx(expected, rel=1e-14,
                                                                     abs=0)


def test_decays_and_other_rates_are_unchanged():
    beta = Reaction.from_names(["n"], ["h1", "electron", "anti-neutrino_e"],
                               rate_fits=[RateFit([-6.78, 0, 0, 0, 0, 0, 0])])
    capture = Reaction.from_names(["c12", "h1"], ["n13", "gamma"],
                                  rate_fits=[RateFit([17.1482, 0.0, -13.692, -0.230881,
                                                      4.44362, -3.15898, -0.666667])])
    for reaction in (beta, capture):
        assert reaction.electron_count == 0
        assert reaction.rate(0.03, rho=150.0, ye=0.86) == reaction.rate(0.03)
        assert reaction.rate(0.03) > 0.0
    # A constant rate is taken as complete, even for an electron capture.
    constant = Reaction.from_names(["be7", "electron"], ["li7", "neutrino_e"],
                                   constant_rate=1e-7)
    assert constant.rate(0.03, rho=150.0, ye=0.86) == pytest.approx(1e-7, abs=0)


def test_ydot_works_out_ye_from_the_abundances():
    network = ReactionNetwork([_be7_capture()])
    abundances = {"h1": 0.7, "he4": 0.07, "be7": 1e-10}
    ye = electron_fraction(abundances)
    assert ye == pytest.approx(0.7 + 2 * 0.07 + 4e-10)
    implicit = network.ydot(abundances, 0.015, rho=150.0)
    explicit = network.ydot(abundances, 0.015, rho=150.0, ye=ye)
    assert implicit["be7"] < 0.0
    assert implicit == pytest.approx(explicit, rel=1e-14, abs=0)
    # ... and that the electron density is in it.
    assert implicit["be7"] == pytest.approx(
        -RateFit(BE7_FIT).rate(0.015) * 150.0 * ye * 1e-10, rel=1e-12, abs=0)


def test_flux_and_diagnostics_work_out_ye_too():
    capture = _be7_capture()
    species = {"h1": Species("h1", 1, 1, 7.289, 0.5),
               "he4": Species("he4", 2, 4, 2.4249, 0.0),
               "be7": Species("be7", 4, 7, 15.769, 1.5),
               "li7": Species("li7", 3, 7, 14.907, 1.5)}
    abundances = {"h1": 0.7, "he4": 0.07, "be7": 1e-10}
    ye = electron_fraction(abundances)
    expected = RateFit(BE7_FIT).rate(0.015) * 150.0 * ye * 1e-10
    assert capture.flux(abundances, 0.015, rho=150.0) == pytest.approx(expected, rel=1e-12,
                                                                       abs=0)
    # The rate of change of Ye, a diagnostic built on the flows, agrees.
    net = Network(species=species, reactions=ReactionNetwork([capture]),
                  zones=[Zone(abundances=abundances)])
    changes = charge_changing_flows(net, 0, t9=0.015, rho=150.0)
    assert changes[capture.string] == pytest.approx(-expected, rel=1e-12, abs=0)
    rates = ReactionNetwork([capture]).rates(0.015, rho=150.0, ye=ye)
    assert rates[capture.string] == pytest.approx(expected / 1e-10, rel=1e-12, abs=0)


def test_be7_decays_at_the_rate_set_by_the_electron_density():
    species = {"h1": Species("h1", 1, 1, 7.289, 0.5),
               "be7": Species("be7", 4, 7, 15.769, 1.5),
               "li7": Species("li7", 3, 7, 14.907, 1.5)}
    net = Network(species=species, reactions=ReactionNetwork([_be7_capture()]))
    t9, rho = 0.015, 150.0
    zone = Zone(abundances={"h1": 0.9999, "be7": 1e-4 / 7})
    # Ye is set by the protons; the 7Be changes it by only 1.4e-5.
    ye = 0.9999 + 4 * 1e-4 / 7
    lam = RateFit(BE7_FIT).rate(t9) * rho * ye
    t_end = 3.0 / lam
    result = evolve_zone(net, zone, np.linspace(0.0, t_end, 7),
                         thermo=constant_thermo(t9, rho), method="bdf",
                         rtol=1e-10, atol=1e-20)
    assert result.success
    be7 = result.y[:, result.species.index("be7")]
    expected = (1e-4 / 7) * np.exp(-lam * result.time)
    np.testing.assert_allclose(be7, expected, rtol=1e-4)
    # The right-hand side NucNetPy builds sees the same rate.
    f = rhs(net, result.species, constant_thermo(t9, rho))
    y0 = np.array([zone.get_abundance(s) for s in result.species])
    dy = f(0.0, y0)[result.species.index("be7")]
    assert dy == pytest.approx(-lam * y0[result.species.index("be7")], rel=1e-4, abs=0)
    assert math.isfinite(dy)
