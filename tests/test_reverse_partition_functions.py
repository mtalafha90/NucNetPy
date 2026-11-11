"""Library reverse rates must carry the partition functions ReacLib leaves out.

ReacLib derives each reverse rate from its forward rate by detailed balance
with the ground-state spins only, and expects the user to multiply it by the
ratio of partition functions.  Evaluated without that factor, the alpha chain
at T9 = 5 relaxes to a composition up to 48 per cent away from nuclear
statistical equilibrium.  These tests pin the correction, the detection of the
derived reverse rates it relies on, and the choice of forward member in
:func:`consistent_reverse_network`, which libnucnet files without Q values had
broken.
"""
from pathlib import Path

import numpy as np
import pytest

from nucnetpy.core import Network, Zone
from nucnetpy.detailed_balance import (_q_value, apply_reverse_partition_functions,
                                       consistent_reverse_network, reverse_rate)
from nucnetpy.io.xml import read_xml, write_xml
from nucnetpy.nse import solve_nse
from nucnetpy.reactions import RateFit, Reaction, ReactionNetwork
from nucnetpy.solver import constant_thermo, evolve_zone, time_grid
from nucnetpy.species import Species

REFERENCE = Path(__file__).resolve().parents[1] / "validation" / "reference"
NETWORK = REFERENCE / "si_burning_network.xml"
needs_reference = pytest.mark.skipif(not NETWORK.is_file(),
                                     reason="reference archive not present")

# 52Fe(a,g)56Ni from the archived network, and its ReacLib reverse.
FORWARD = [66.6417, 0.0, -91.6819, -9.51885, -0.533014, 0.0892607, -0.666667]
Q_FORWARD = 7.997


def _reaclib_reverse(coeffs, q, a6_shift, a0_shift=25.0):
    """Apply ReacLib's forward-to-reverse conversion to one set of coefficients."""
    rev = list(coeffs)
    rev[0] += a0_shift
    rev[1] -= 11.6045 * q
    rev[6] += a6_shift
    return rev


def _toy_network(reverse_coeffs, reverse_products=("fe52", "he4")):
    species = {
        "ni56": Species("ni56", 28, 56, -53.904, 0.0,
                        partition={1.0: 1.0, 5.0: 1.2, 10.0: 2.0}),
        "fe52": Species("fe52", 26, 52, -48.332, 0.0,
                        partition={1.0: 1.0, 5.0: 2.1, 10.0: 5.0}),
        "he4": Species("he4", 2, 4, 2.4249, 0.0),
    }
    forward = Reaction.from_names(["fe52", "he4"], ["ni56", "gamma"],
                                  rate_fits=[RateFit(FORWARD)])
    reverse = Reaction.from_names(["ni56", "gamma"], list(reverse_products),
                                  rate_fits=[RateFit(reverse_coeffs)])
    return Network(species=species,
                   reactions=ReactionNetwork([forward, reverse])), forward, reverse


def test_derived_reverse_rate_gets_the_partition_function_ratio():
    rev = _reaclib_reverse(FORWARD, Q_FORWARD, 1.5)
    net, forward, reverse = _toy_network(rev)
    fitted = RateFit(rev).rate(5.0)
    assert apply_reverse_partition_functions(net) == [reverse]
    assert forward.partition_ratio is None
    # G(fe52) G(he4) / G(ni56) at T9 = 5, interpolated as the NSE solver does.
    assert reverse.rate(5.0) == pytest.approx(fitted * 2.1 / 1.2, rel=1e-12)
    assert reverse.rate(7.5) == pytest.approx(RateFit(rev).rate(7.5) * 3.55 / 1.6,
                                              rel=1e-12)


def test_independent_reverse_fit_is_left_alone():
    rev = _reaclib_reverse(FORWARD, Q_FORWARD, 1.5)
    rev[3] += 0.5            # no longer shares a2..a5 with the forward fit
    net, _, reverse = _toy_network(rev)
    assert apply_reverse_partition_functions(net) == []
    assert reverse.partition_ratio is None


def test_wrong_power_of_temperature_is_left_alone():
    net, _, reverse = _toy_network(_reaclib_reverse(FORWARD, Q_FORWARD, 0.0))
    assert apply_reverse_partition_functions(net) == []


def test_reaction_with_a_lepton_is_left_alone():
    species = {"n": Species("n", 0, 1, 8.0713, 0.5),
               "h1": Species("h1", 1, 1, 7.2890, 0.5)}
    decay = Reaction.from_names(["n"], ["h1", "electron", "anti-neutrino_e"],
                                rate_fits=[RateFit([-6.78, 0, 0, 0, 0, 0, 0])])
    capture = Reaction.from_names(["h1", "electron", "anti-neutrino_e"], ["n"],
                                  rate_fits=[RateFit([-16.2, 0, 0, 0, 0, 0, 0])])
    net = Network(species=species, reactions=ReactionNetwork([decay, capture]))
    assert apply_reverse_partition_functions(net) == []


@needs_reference
def test_every_reverse_rate_of_the_reference_network_is_corrected():
    net = read_xml(str(NETWORK))
    corrected = [r for r in net.reactions.reactions if r.partition_ratio is not None]
    assert len(corrected) == 16          # every forward/reverse pair
    assert all(_q_value(r, net.species) < 0 for r in corrected)


@needs_reference
@pytest.mark.parametrize("t9", [2.0, 5.0, 7.0])
def test_corrected_reverse_rates_obey_detailed_balance(t9):
    net = read_xml(str(NETWORK))
    raw = read_xml(str(NETWORK), reverse_partition_functions=False)
    by_key = {r.key: r for r in net.reactions.reactions}
    raw_by_key = {r.key: r for r in raw.reactions.reactions}
    fixed, uncorrected = [], []
    for r in net.reactions.reactions:
        if r.partition_ratio is None:
            continue
        db = reverse_rate(by_key[(r.key[1], r.key[0])], net.species, t9)
        fixed.append(r.rate(t9) / db)
        uncorrected.append(raw_by_key[r.key].rate(t9) / db)
    # What is left is the library's own Q values and the error of its fits.
    assert max(abs(np.log(fixed))) < np.log(1.03)
    if t9 >= 5.0:
        assert max(abs(np.log(uncorrected))) > np.log(1.5)


@needs_reference
def test_correction_survives_an_xml_round_trip(tmp_path):
    net = read_xml(str(NETWORK))
    path = tmp_path / "round_trip.xml"
    write_xml(net, path)
    back = read_xml(str(path))
    rates = {r.key: r.rate(5.0) for r in net.reactions.reactions}
    for r in back.reactions.reactions:
        assert r.rate(5.0) == pytest.approx(rates[r.key], rel=1e-12)


@needs_reference
def test_consistent_reverse_network_keeps_the_exothermic_member():
    net = consistent_reverse_network(read_xml(str(NETWORK)))
    kept = [r for r in net.reactions.reactions if r.source != "detailed_balance"]
    rebuilt = [r for r in net.reactions.reactions if r.source == "detailed_balance"]
    assert len(rebuilt) == 16
    assert all(_q_value(r, net.species) > 0 for r in kept
               if (r.key[1], r.key[0]) in {x.key for x in rebuilt})


@needs_reference
def test_alpha_chain_relaxes_to_nse_within_a_per_cent():
    def median_offset(reverse_partition_functions):
        net = read_xml(str(NETWORK),
                       reverse_partition_functions=reverse_partition_functions)
        zone = net.zone(0)
        ye = zone.ye(net.species)
        result = evolve_zone(net, zone, time_grid(0.0, 10.0, 40),
                             thermo=constant_thermo(5.0, 1.0e8),
                             method="bdf", rtol=1e-8, atol=1e-14)
        assert result.success
        nse = solve_nse(net, t9=5.0, rho=1.0e8, ye=ye,
                        species=[n for n in net.species if n not in {"n", "h1"}])
        final = result.final_abundances
        offsets = [abs(final.get(k, 0.0) - v) / v
                   for k, v in nse.abundances.items()
                   if k not in {"n", "h1"} and net.species[k].a * v >= 1e-6]
        return float(np.median(offsets))

    assert median_offset(True) < 0.01
    assert median_offset(False) > 0.03


# --- Which member of a pair is the derived reverse rate ---------------------

#: A capture whose product is unbound, so the capture is endothermic
#: (Q = -0.75 MeV), like 24O(n,g)25O near the neutron drip line.  ReacLib
#: still fits the capture and derives the photodisintegration from it.
UNBOUND_CAPTURE = [12.0, 0.0, -1.5, 0.8, -0.1, 0.01, 0.4]
Q_UNBOUND = -0.75


def _unbound_capture_pair(capture_label="", reverse_label=""):
    species = {
        "o24": Species("o24", 8, 24, 18.5, 0.0, partition={1.0: 1.0, 5.0: 1.0}),
        "o25": Species("o25", 8, 25, 27.3, 1.5, partition={1.0: 1.0, 5.0: 1.6}),
        "n": Species("n", 0, 1, 8.0713, 0.5),
    }
    capture = Reaction.from_names(["o24", "n"], ["o25", "gamma"],
                                  rate_fits=[RateFit(UNBOUND_CAPTURE, label=capture_label)])
    reverse = Reaction.from_names(
        ["o25", "gamma"], ["o24", "n"],
        rate_fits=[RateFit(_reaclib_reverse(UNBOUND_CAPTURE, Q_UNBOUND, 1.5),
                           label=reverse_label)])
    net = Network(species=species, reactions=ReactionNetwork([capture, reverse]))
    return net, capture, reverse


def test_photodisintegration_is_the_derived_member_even_when_exothermic():
    # The photodisintegration has the higher a1 here, so the old rule (the
    # endothermic member is the derived one) corrected the capture instead.
    net, capture, reverse = _unbound_capture_pair()
    assert apply_reverse_partition_functions(net) == [reverse]
    assert capture.partition_ratio is None


def test_reaclib_flag_decides_when_the_labels_carry_it():
    # A two-body pair whose forward rate is the endothermic one: only the
    # library's v flag can tell which member was derived.
    species = {
        "he4": Species("he4", 2, 4, 2.4249, 0.0),
        "ne22": Species("ne22", 10, 22, -8.0247, 0.0, partition={1.0: 1.0, 5.0: 1.3}),
        "mg25": Species("mg25", 12, 25, -13.9336, 2.5, partition={1.0: 1.0, 5.0: 1.5}),
        "n": Species("n", 0, 1, 8.0713, 0.5),
    }
    forward = [20.0, -5.0, -40.0, 3.0, -0.2, 0.01, -0.67]
    q = -0.478                                     # 22Ne(a,n)25Mg is endothermic
    fitted = Reaction.from_names(["ne22", "he4"], ["mg25", "n"],
                                 rate_fits=[RateFit(forward, label="ka02r")])
    derived = Reaction.from_names(["mg25", "n"], ["ne22", "he4"],
                                  rate_fits=[RateFit(_reaclib_reverse(forward, q, 0.0),
                                                     label="ka02rv")])
    net = Network(species=species, reactions=ReactionNetwork([fitted, derived]))
    assert apply_reverse_partition_functions(net) == [derived]
    assert fitted.partition_ratio is None


def test_without_flags_a_two_body_pair_falls_back_to_the_endothermic_member():
    species = {"he4": Species("he4", 2, 4, 2.4249, 0.0),
               "ne22": Species("ne22", 10, 22, -8.0247, 0.0),
               "mg25": Species("mg25", 12, 25, -13.9336, 2.5),
               "n": Species("n", 0, 1, 8.0713, 0.5)}
    forward = [20.0, -5.0, -40.0, 3.0, -0.2, 0.01, -0.67]
    fitted = Reaction.from_names(["ne22", "he4"], ["mg25", "n"],
                                 rate_fits=[RateFit(forward)])
    derived = Reaction.from_names(["mg25", "n"], ["ne22", "he4"],
                                  rate_fits=[RateFit(_reaclib_reverse(forward, -0.478, 0.0))])
    net = Network(species=species, reactions=ReactionNetwork([fitted, derived]))
    # With no flag and no photon the endothermic member is taken to be the
    # derived one; here that is the fitted rate.  This is the documented limit.
    assert apply_reverse_partition_functions(net) == [fitted]


def test_rounded_a6_coefficients_are_still_recognised():
    # 29Si(g,p)28Al in the 2022 ReacLib snapshot: a6 = 98.8098 for the capture
    # and 100.31 for the derived rate, rounded from 100.3098.
    capture_fit = [82.7866, -3.1339, 114.903, -203.401, 11.5793, -0.6728, 98.8098]
    reverse_fit = [107.692, -146.264, 114.903, -203.401, 11.5793, -0.6728, 100.31]
    species = {"al28": Species("al28", 13, 28, -16.85, 3.0, partition={1.0: 1.0, 5.0: 2.0}),
               "h1": Species("h1", 1, 1, 7.289, 0.5),
               "si29": Species("si29", 14, 29, -21.895, 0.5, partition={1.0: 1.0, 5.0: 1.1})}
    capture = Reaction.from_names(["al28", "h1"], ["si29", "gamma"],
                                  rate_fits=[RateFit(capture_fit)])
    reverse = Reaction.from_names(["si29", "gamma"], ["al28", "h1"],
                                  rate_fits=[RateFit(reverse_fit)])
    net = Network(species=species, reactions=ReactionNetwork([capture, reverse]))
    assert apply_reverse_partition_functions(net) == [reverse]
