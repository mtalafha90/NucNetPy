"""The adaptive implicit Euler method of evolve_zone.

``method="adaptive_implicit_euler"`` takes implicit Euler steps with
Newton-Raphson iterations and sizes each step by the fractional change of
the abundances, as NucNet Tools does.  It is first-order accurate and cannot
turn a decaying abundance negative, which keeps it working as a fuel runs out
in large networks, where SciPy's BDF method combined with NucNetPy's clipping
of negative abundances lets the step size collapse.
"""
from pathlib import Path

import numpy as np
import pytest

from nucnetpy.core import Network, Zone
from nucnetpy.io.xml import read_xml
from nucnetpy.reactions import RateFit, Reaction, ReactionNetwork
from nucnetpy.solver import METHODS, constant_thermo, evolve_zone
from nucnetpy.species import Species

REFERENCE = Path(__file__).resolve().parents[1] / "validation" / "reference" / "si_burning_network.xml"
METHOD = "adaptive_implicit_euler"


def _network(reactions, names):
    species = {n: Species.parse(n) for n in names}
    return Network(species=species, reactions=ReactionNetwork(reactions))


def _constant(value):
    return [float(np.log(value)), 0, 0, 0, 0, 0, 0]


def test_method_is_offered():
    assert METHOD in METHODS


def test_decay_chain_is_first_order_accurate():
    # he4 -> 4 h1 and h1 -> n, two decays with an exact solution.  The names
    # only label the abundances; mass conservation does not matter here.
    l1, l2 = 1.0, 0.25
    net = _network([Reaction.from_names(["he4"], ["h1", "h1", "h1", "h1"],
                                        rate_fits=[RateFit(_constant(l1))]),
                    Reaction.from_names(["h1"], ["n"],
                                        rate_fits=[RateFit(_constant(l2))])],
                   ["he4", "h1", "n"])
    times = np.concatenate(([0.0], np.geomspace(1e-3, 8.0, 40)))
    exact_b = 4.0 * l1 / (l2 - l1) * (np.exp(-l1 * times) - np.exp(-l2 * times))

    def error(max_change):
        r = evolve_zone(net, Zone(abundances={"he4": 1.0}), times,
                        thermo=constant_thermo(1.0, 1.0), method=METHOD,
                        max_change=max_change)
        assert r.success, r.message
        assert np.all(r.y >= 0.0)
        b = r.y[:, r.species.index("h1")]
        return np.max(np.abs(b[1:] - exact_b[1:]) / exact_b[1:])

    coarse, fine = error(0.02), error(0.01)
    assert fine < 0.02
    assert coarse / fine == pytest.approx(2.0, rel=0.2)       # first order


def test_burning_to_exhaustion_follows_the_exact_power_law():
    # he4 + he4 -> be8 with flow k Y^2 / 2 gives Y(t) = Y0 / (1 + k Y0 t):
    # the fuel falls through ten decades and must neither stall nor go
    # negative.
    k = 1.0
    net = _network([Reaction.from_names(["he4", "he4"], ["be8"],
                                        rate_fits=[RateFit(_constant(k))])],
                   ["he4", "be8"])
    times = np.concatenate(([0.0], np.geomspace(1e-6, 1e10, 50)))
    r = evolve_zone(net, Zone(abundances={"he4": 0.25}), times,
                    thermo=constant_thermo(1.0, 1.0), method=METHOD)
    assert r.success, r.message
    fuel = r.y[:, r.species.index("he4")]
    exact = 0.25 / (1.0 + k * 0.25 * times)
    np.testing.assert_allclose(fuel, exact, rtol=0.01)
    a = np.array([4.0, 8.0])[[r.species.index("he4"), r.species.index("be8")]]
    masses = r.y[:, [r.species.index("he4"), r.species.index("be8")]] @ a
    np.testing.assert_allclose(masses, 1.0, atol=1e-12)


@pytest.mark.skipif(not REFERENCE.is_file(), reason="reference archive not present")
def test_alpha_chain_reaches_the_same_equilibrium_as_bdf():
    net = read_xml(str(REFERENCE))
    times = np.concatenate(([0.0], np.geomspace(1e-8, 10.0, 30)))
    thermo = constant_thermo(5.0, 1.0e8)
    bdf = evolve_zone(net, net.zone(0), times, thermo=thermo, method="bdf",
                      rtol=1e-8, atol=1e-14)
    # The equilibrium does not depend on the step size, so coarse steps do.
    euler = evolve_zone(net, net.zone(0), times, thermo=thermo, method=METHOD,
                        max_change=0.1)
    assert bdf.success and euler.success, euler.message
    assert euler.steps > 0 and euler.nfev > 0 and euler.njev > 0
    np.testing.assert_array_equal(euler.time, times)
    a = np.array([net.species[s].a for s in euler.species])
    x_bdf, x_euler = bdf.y[-1] * a, euler.y[-1] * a
    big = x_bdf > 1e-6
    np.testing.assert_allclose(x_euler[big], x_bdf[big], rtol=1e-3)
    assert abs(x_euler.sum() - 1.0) < 1e-10


def test_time_dependent_conditions_are_followed():
    # A decay whose rate doubles with T9, while T9 rises linearly: Y = exp(-t^2).
    fit = [0.0, 0, 0, 0, 0, 0, 1.0]                  # lambda = T9
    net = _network([Reaction.from_names(["he4"], ["h1", "h1", "h1", "h1"],
                                        rate_fits=[RateFit(fit)])], ["he4", "h1"])
    times = np.linspace(0.0, 2.0, 21)
    r = evolve_zone(net, Zone(abundances={"he4": 0.25}), times,
                    thermo=lambda t, y: (2.0 * t + 1e-12, 1.0), method=METHOD,
                    max_change=0.002)
    assert r.success, r.message
    fuel = r.y[:, r.species.index("he4")]
    np.testing.assert_allclose(fuel, 0.25 * np.exp(-times ** 2), rtol=0.02)


def test_very_short_times_are_reached_by_real_steps():
    # A decay over 2e-15 s.  Each output time must be reached by integrating
    # up to it, however small the times are.
    lam = 1.0e14
    net = _network([Reaction.from_names(["he4"], ["h1", "h1", "h1", "h1"],
                                        rate_fits=[RateFit(_constant(lam))])],
                   ["he4", "h1"])
    times = np.array([0.0, 1.0e-15, 2.0e-15])
    r = evolve_zone(net, Zone(abundances={"he4": 0.25}), times,
                    thermo=constant_thermo(1.0, 1.0), method=METHOD)
    assert r.success, r.message
    fuel = r.y[:, r.species.index("he4")]
    np.testing.assert_allclose(fuel, 0.25 * np.exp(-lam * times), rtol=2e-3)
    assert r.steps > 20


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.1])
def test_step_limit_must_be_a_fraction(bad):
    net = _network([], ["he4"])
    with pytest.raises(ValueError):
        evolve_zone(net, Zone(abundances={"he4": 0.25}), [0.0, 1.0],
                    method=METHOD, max_change=bad)


def test_times_must_increase():
    net = _network([], ["he4"])
    with pytest.raises(ValueError):
        evolve_zone(net, Zone(abundances={"he4": 0.25}), [0.0, 1.0, 1.0],
                    method=METHOD)
