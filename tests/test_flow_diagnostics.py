"""Flow diagnostics on networks that list reactions in both directions.

JINA networks list each reaction's inverse explicitly.  These tests guard
against counting such a pair twice, against results set by an arbitrary floor
for absent species, and against photons being reported as species.
"""
import copy
import math

import numpy as np
import pytest

from nucnetpy import Network, Species, Reaction, RateFit, Zone, constant_thermo, evolve_zone, time_grid
from nucnetpy.analysis import entropy_generation_rate, integrated_currents, system_timescales
from nucnetpy.detailed_balance import consistent_reverse_network
from nucnetpy.nse import _log_prefactor


def alpha_network():
    net = Network()
    for name, me in [("he4", 2.4249), ("be8", 4.9416), ("c12", 0.0), ("o16", -4.7370)]:
        net.add_species(Species.parse(name, mass_excess=me, spin=0.0))
    for reactants, products, a0 in [(["he4", "he4"], ["be8"], 2.0),
                                     (["be8", "he4"], ["c12"], 3.0),
                                     (["c12", "he4"], ["o16"], 1.0)]:
        net.reactions.add(Reaction.from_names(reactants, products,
                                              rate_fits=[RateFit([a0, 0, 0, 0, 0, 0, 0])]))
    return net


def with_both_directions(net):
    """The same network with every inverse listed explicitly, as in JINA.

    The inverses use the exact detailed-balance rate; a tabulated one would
    differ by its interpolation error.
    """
    from nucnetpy.detailed_balance import _exact_reverse_reaction
    out = copy.deepcopy(net)
    for r in list(net.reactions.reactions):
        out.reactions.add(_exact_reverse_reaction(r, net.species, True))
    return out


def ground_truth(net, abundances, t9, rho):
    """-sum_i (mu_i/kT) dY_i/dt for the network exactly as it is evolved."""
    dy = net.reactions.ydot(abundances, t9=t9, rho=rho)
    return -sum((math.log(abundances[k]) - _log_prefactor(net.species[k], t9, rho)) * v
                for k, v in dy.items())


COMPOSITION = {"he4": 0.2, "be8": 1e-9, "c12": 0.01, "o16": 0.005}


def test_entropy_counts_a_listed_pair_once():
    net = with_both_directions(alpha_network())
    assert len(net.reactions.reactions) == 6
    net.zones = [Zone(abundances=dict(COMPOSITION))]
    truth = ground_truth(net, COMPOSITION, 3.0, 1e6)
    assert entropy_generation_rate(net, 0, t9=3.0, rho=1e6) == pytest.approx(truth, rel=1e-10)
    assert entropy_generation_rate(net, 0, t9=3.0, rho=1e6, use_reverse=False) == pytest.approx(truth, rel=1e-10)


def test_entropy_with_implied_reverses_matches_the_explicit_pair():
    # Listing the detailed-balance inverse explicitly, or leaving it implied,
    # must give the same answer.
    implied = alpha_network()
    implied.zones = [Zone(abundances=dict(COMPOSITION))]
    explicit = with_both_directions(alpha_network())
    explicit.zones = [Zone(abundances=dict(COMPOSITION))]
    a = entropy_generation_rate(implied, 0, t9=3.0, rho=1e6)
    b = entropy_generation_rate(explicit, 0, t9=3.0, rho=1e6)
    assert a == pytest.approx(b, rel=1e-6)
    assert a > 0.0


def test_entropy_does_not_depend_on_a_floor_for_absent_species():
    net = alpha_network()
    absent = dict(COMPOSITION, be8=0.0)
    trace = dict(COMPOSITION, be8=1e-200)
    net.zones = [Zone(abundances=absent), Zone(abundances=trace)]
    s_absent = entropy_generation_rate(net, 0, t9=3.0, rho=1e6)
    s_trace = entropy_generation_rate(net, 1, t9=3.0, rho=1e6)
    assert math.isfinite(s_absent)
    assert s_absent == pytest.approx(s_trace, rel=1e-12)
    # Only c12 + he4 -> o16 involves no absent species.
    only = Network(species=net.species)
    only.reactions.add(net.reactions.reactions[2])
    only.zones = [Zone(abundances=absent)]
    assert s_absent == pytest.approx(entropy_generation_rate(only, 0, t9=3.0, rho=1e6), rel=1e-12)


def test_photons_are_not_reported_as_species():
    net = alpha_network()
    net.reactions.add(Reaction.from_names(["o16", "he4"], ["ne20", "gamma"], constant_rate=1.0))
    net.add_species(Species.parse("ne20", mass_excess=-7.0419))
    net.zones = [Zone(abundances=dict(COMPOSITION, ne20=1e-6))]
    assert "gamma" not in net.reactions.ydot(COMPOSITION, t9=3.0, rho=1e6)
    assert "gamma" not in system_timescales(net, 0, t9=3.0, rho=1e6)


def test_integrated_currents_add_repeated_listings():
    net = Network()
    for name in ["ni56", "co56"]:
        net.add_species(Species.parse(name))
    for _ in range(2):          # the same decay listed twice
        net.reactions.add(Reaction.from_names(["ni56"], ["co56"], constant_rate=1.5))
    zone = Zone(abundances={"ni56": 1e-2, "co56": 0.0})
    thermo = constant_thermo(1.0, 1.0)
    res = evolve_zone(net, zone, time_grid(0.0, 0.3, 601), thermo=thermo, method="rk4")
    currents = integrated_currents(net, res, thermo, use_reverse=False)
    assert currents["ni56 -> co56"] == pytest.approx(res.final_abundances["co56"], rel=1e-6)
