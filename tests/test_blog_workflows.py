"""Tests for the analysis workflows from the NucNet Tools blog.

Each test corresponds to one or more posts on Bradley Meyer's SourceForge
blog (https://sourceforge.net/u/mbradle/blog/); see docs/BLOG_COVERAGE.md for the
full post-to-feature map.
"""
import math

import numpy as np
import pytest

from nucnetpy import (
    Network, Species, Zone, Reaction, RateFit,
    solve_nse, evolve_zone, time_grid, constant_thermo,
    log_equilibrium_constant, reverse_rate, reverse_reaction, net_flows,
    charge_changing_flows, system_timescales, heavy_nuclei_abundance,
    neutron_exposure, entropy_generation_rate, separation_energy,
    fission_reaction, select_species, limit_network,
)
from nucnetpy.constants import KB_MEV


def alpha_network():
    net = Network()
    for name, me in [("he4", 2.4249), ("be8", 4.9416), ("c12", 0.0), ("o16", -4.7370)]:
        net.add_species(Species.parse(name, mass_excess=me))
    net.reactions.add(Reaction.from_names(["he4", "he4"], ["be8"],
                                          rate_fits=[RateFit([2.0, 0, 0, 0, 0, 0, 0])], q_value=-0.0918))
    net.reactions.add(Reaction.from_names(["be8", "he4"], ["c12"],
                                          rate_fits=[RateFit([3.0, 0, 0, 0, 0, 0, 0])], q_value=7.3666))
    net.reactions.add(Reaction.from_names(["c12", "he4"], ["o16"],
                                          rate_fits=[RateFit([1.0, 0, 0, 0, 0, 0, 0])], q_value=7.1616))
    return net


# --- Comparing forward and reverse reaction rates / computing reaction flows ---

def test_equilibrium_constant_contains_q_value():
    # d(ln K)/d(1/kT) = Q: check exp(Q/kT) dominates the T-dependence ratio.
    net = alpha_network()
    rxn = net.reactions.reactions[1]  # be8 + he4 -> c12, Q = 7.3666 MeV
    t9a, t9b = 3.0, 6.0
    ka = log_equilibrium_constant(rxn, net.species, t9a, rho=1e6)
    kb = log_equilibrium_constant(rxn, net.species, t9b, rho=1e6)
    q = sum(net.species[p.species].mass_excess * -p.count for p in rxn.products)
    q += sum(net.species[r.species].mass_excess * r.count for r in rxn.reactants)
    expected = q / (KB_MEV * 1e9) * (1.0 / t9a - 1.0 / t9b)
    # the translational/partition part varies slowly compared with exp(Q/kT)
    assert (ka - kb) == pytest.approx(expected, rel=0.15)


def test_reverse_rate_is_density_independent():
    net = alpha_network()
    rxn = net.reactions.reactions[2]
    r1 = reverse_rate(rxn, net.species, 4.0, rho=1.0)
    r2 = reverse_rate(rxn, net.species, 4.0, rho=1e9)
    assert r1 == pytest.approx(r2, rel=1e-10)


def test_net_flows_vanish_at_nse():
    # Detailed balance and the NSE solver must agree: at NSE abundances the
    # net flux of every Z,A-conserving reaction is zero.
    net = alpha_network()
    t9, rho = 4.0, 1.0e6
    nse = solve_nse(net, t9=t9, rho=rho, ye=0.5)
    assert nse.success
    result = net_flows(net, nse.abundances, t9=t9, rho=rho)
    for reaction, (fwd, rev, net_flux) in result.items():
        assert fwd > 0.0, reaction
        assert net_flux == pytest.approx(0.0, abs=1e-8 * fwd), reaction


def test_reverse_reaction_object():
    net = alpha_network()
    rxn = net.reactions.reactions[1]              # be8 + he4 -> c12
    rev = reverse_reaction(rxn, net.species)      # c12 -> be8 + he4
    assert rev.q_value == -rxn.q_value
    assert {p.species for p in rev.reactants} == {"c12"}
    assert {p.species for p in rev.products} == {"be8", "he4"}
    assert rev.rate(4.0) == pytest.approx(reverse_rate(rxn, net.species, 4.0), rel=0.05)


# --- Computing charge-changing flows -----------------------------------------

def test_charge_changing_flows():
    net = Network()
    for name in ["n", "h1", "he4"]:
        net.add_species(Species.parse(name))
    net.reactions.add(Reaction.from_names(["n"], ["h1"], constant_rate=2.0))       # beta: dZ=+1
    net.reactions.add(Reaction.from_names(["he4"], ["he4"], constant_rate=5.0))    # no-op: dZ=0
    z = Zone(abundances={"n": 0.3, "h1": 0.1, "he4": 0.15})
    net.add_zone(z)
    contributions = charge_changing_flows(net, 0, t9=1.0, rho=1.0)
    assert list(contributions) == ["n -> h1"]
    assert contributions["n -> h1"] == pytest.approx(2.0 * 0.3)
    assert sum(contributions.values()) == pytest.approx(0.6)  # dYe/dt


# --- Computing system timescales ----------------------------------------------

def test_system_timescales():
    net = Network()
    for name in ["n", "h1"]:
        net.add_species(Species.parse(name))
    net.reactions.add(Reaction.from_names(["n"], ["h1"], constant_rate=4.0))
    z = Zone(abundances={"n": 0.2, "h1": 0.0})
    net.add_zone(z)
    ts = system_timescales(net, 0, t9=1.0, rho=1.0)
    assert ts["n"] == pytest.approx(0.25)  # Y/|dY/dt| = 0.2/(4*0.2)
    assert ts["h1"] == pytest.approx(0.0)  # Y=0, growing


# --- Computing the number of heavy nuclei --------------------------------------

def test_heavy_nuclei_abundance():
    net = alpha_network()
    z = Zone(abundances={"he4": 0.2, "be8": 0.01, "c12": 0.003, "o16": 0.001})
    assert heavy_nuclei_abundance(z, net.species, zmin=3) == pytest.approx(0.014)
    assert heavy_nuclei_abundance(z, net.species, zmin=6) == pytest.approx(0.004)


# --- Computing the s-process neutron exposure ----------------------------------

def test_neutron_exposure_constant_state():
    # With constant Y_n, rho, T the exposure is n_n * v_T * dt analytically.
    net = Network()
    for name in ["n", "fe56"]:
        net.add_species(Species.parse(name))
    z = Zone(abundances={"n": 1e-8, "fe56": 1e-3})
    net.add_zone(z)
    times = time_grid(0.0, 100.0, 11)
    thermo = constant_thermo(0.3, 3.0e3)
    res = evolve_zone(net, z, times, thermo=thermo, method="rk4")  # no reactions: static
    tau = neutron_exposure(res, thermo)
    from nucnetpy.constants import AVOGADRO, KB_CGS, MN_G
    n_n = 3.0e3 * AVOGADRO * 1e-8
    v_t = math.sqrt(2.0 * KB_CGS * 0.3e9 / MN_G)
    assert tau == pytest.approx(n_n * v_t * 100.0 * 1e-27, rel=1e-9)


# --- Computing the entropy generation rate --------------------------------------

def test_entropy_generation_rate_positive_toward_equilibrium():
    net = alpha_network()
    z = Zone(abundances={"he4": 0.25, "be8": 1e-12, "c12": 1e-12, "o16": 1e-12})
    net.add_zone(z)
    # far from equilibrium, flowing downhill: entropy generation must be > 0
    assert entropy_generation_rate(net, 0, t9=3.0, rho=1e6) > 0.0


# --- Computing separation energies ----------------------------------------------

def test_separation_energy_neutron_and_proton():
    species_map = {}
    for name, me in [("n", 8.0713), ("h1", 7.2890), ("ni56", -53.90),
                     ("ni55", -45.33), ("co55", -54.03)]:
        species_map[name] = Species.parse(name, mass_excess=me)
    s_n = separation_energy("ni56", species_map, particle="n")
    s_p = separation_energy("ni56", species_map, particle="p")
    assert s_n == pytest.approx(-45.33 + 8.0713 - (-53.90))   # ME(ni55)+ME(n)-ME(ni56)
    assert s_p == pytest.approx(-54.03 + 7.2890 - (-53.90))   # ME(co55)+ME(h1)-ME(ni56)
    # the proton daughter must be Z-1 (co55), not the same element
    assert s_p != pytest.approx(s_n)


# --- Adding fission to an r-process calculation ----------------------------------

def test_fission_reaction_conserves_a_z():
    species_map = {n: Species.parse(n) for n in ["n", "u238", "zr100", "te136"]}
    # 238 = 100 + 136 + 2 neutrons; 92 = 40 + 52
    rxn = fission_reaction("u238", ["zr100", "te136"], neutrons=2, rate=1e-3)
    ok, da, dz = rxn.conserves_a_z(species_map)
    assert ok, (da, dz)
    assert rxn.constant_rate == pytest.approx(1e-3)
    assert rxn.source == "fission"


# --- Selecting an input network (nuclide XPath) -----------------------------------

def test_select_species_and_limit():
    net = alpha_network()
    net.add_species(Species.parse("ni56"))
    sel = select_species(net, zmin=4, zmax=8)
    assert sel == ["be8", "c12", "o16"]
    sel_el = select_species(net, elements=["he", "ni"])
    assert sel_el == ["he4", "ni56"]
    limit_network(net, select_species(net, zmax=6))
    assert set(net.species) == {"he4", "be8", "c12"}
    assert all(all(p.species in net.species for p in r.reactants + r.products)
               for r in net.reactions.reactions)
