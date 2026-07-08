"""Tests for features ported from the NucNet Tools C++ source (r647).

Covers the previously "not ported" items from BLOG_COVERAGE.md, implemented
from user/nse_corr.cpp, user/flow_utilities.cpp, and the libnuceq cluster
workflow of examples/analysis/compare_equil.cpp.
"""
import math

import numpy as np
import pytest

from nucnetpy import (
    Network, Species, Zone, Reaction, RateFit,
    solve_nse, solve_qse, QSECluster, cluster_abundance, cluster_ydot,
    coulomb_nse_correction, gamma_e,
    species_coulomb_chemical_potential, species_coulomb_energy, species_coulomb_entropy,
    entropy_generation_rate, integrated_currents,
    evolve_zone, time_grid, constant_thermo, heavy_nuclei_abundance,
)


def iron_network():
    net = Network()
    for name, me in [("n", 8.0713), ("h1", 7.2890), ("he4", 2.4249),
                     ("si28", -21.4928), ("fe56", -60.6054), ("ni56", -53.9042)]:
        net.add_species(Species.parse(name, mass_excess=me))
    return net


# --- Coulomb corrections to NSE (user/nse_corr.cpp) ---------------------------

def test_coulomb_branches_match_at_gamma_one():
    # find (t9, rho) putting Gamma_i exactly at 1 for Z=1, then check the
    # two Bravo & Garcia-Senz branches agree in value there (they are matched
    # analytically in value and first derivative).
    t9, ye = 1.0, 0.5
    lo, hi = 1e0, 1e12
    for _ in range(200):
        rho = math.sqrt(lo * hi)
        if gamma_e(t9, rho, ye) > 1.0:
            hi = rho
        else:
            lo = rho
    rho = math.sqrt(lo * hi)
    g = gamma_e(t9, rho, ye)
    assert g == pytest.approx(1.0, rel=1e-6)
    below = species_coulomb_chemical_potential(1, t9, rho * 0.999, ye)
    above = species_coulomb_chemical_potential(1, t9, rho * 1.001, ye)
    assert below == pytest.approx(above, rel=1e-2)


def test_coulomb_neutron_uncorrected_and_sign():
    t9, rho, ye = 3.0, 1e9, 0.5
    assert species_coulomb_chemical_potential(0, t9, rho, ye) == 0.0
    assert species_coulomb_energy(0, t9, rho, ye) == 0.0
    assert species_coulomb_entropy(0, t9, rho, ye) == 0.0
    # in the strong-coupling regime the Coulomb chemical potential is negative
    # (plasma binding), so the NSE correction factor is positive.
    mu = species_coulomb_chemical_potential(26, t9, rho, ye)
    assert mu < 0.0
    sp = Species.parse("fe56")
    assert coulomb_nse_correction(sp, t9, rho, ye) == pytest.approx(-mu)


def test_coulomb_correction_shifts_nse_toward_heavies():
    # Coulomb binding lowers the free energy of high-Z nuclei, so switching the
    # correction on at high density must increase the iron-peak share.
    net = iron_network()
    t9, rho, ye = 5.5, 1.0e9, 0.5
    plain = solve_nse(net, t9, rho, ye)
    corr = solve_nse(net, t9, rho, ye, nse_correction=coulomb_nse_correction)
    assert plain.success and corr.success
    assert abs(corr.xsum - 1.0) < 1e-6 and abs(corr.computed_ye - ye) < 1e-6
    heavies_plain = plain.abundances["ni56"] + plain.abundances["fe56"]
    heavies_corr = corr.abundances["ni56"] + corr.abundances["fe56"]
    assert heavies_corr > heavies_plain


# --- QSE cluster equilibria (libnuceq) -----------------------------------------

def test_qse_no_clusters_equals_nse():
    net = iron_network()
    t9, rho, ye = 5.0, 1e8, 0.5
    nse = solve_nse(net, t9, rho, ye)
    qse = solve_qse(net, t9, rho, ye, clusters=[])
    assert qse.success
    for name, y in nse.abundances.items():
        assert qse.abundances[name] == pytest.approx(y, rel=1e-6)


def test_qse_cluster_constraint_is_enforced():
    net = iron_network()
    t9, rho, ye = 5.0, 1e8, 0.5
    nse = solve_nse(net, t9, rho, ye)
    cluster_species = ["si28", "fe56", "ni56"]
    y_nse = cluster_abundance(nse.abundances, cluster_species)
    # constrain the heavy cluster to half its NSE value
    target = 0.5 * y_nse
    qse = solve_qse(net, t9, rho, ye,
                    clusters=[QSECluster(cluster_species, target, label="heavies")])
    assert qse.success
    assert cluster_abundance(qse.abundances, cluster_species) == pytest.approx(target, rel=1e-6)
    assert qse.xsum == pytest.approx(1.0, abs=1e-6)
    assert qse.computed_ye == pytest.approx(ye, abs=1e-6)
    # the cluster multiplier must be nonzero when the constraint binds
    assert abs(qse.lambdas[0]) > 1e-3


def test_qse_with_nse_constraint_recovers_nse():
    net = iron_network()
    t9, rho, ye = 5.0, 1e8, 0.5
    nse = solve_nse(net, t9, rho, ye)
    cluster_species = ["si28", "fe56", "ni56"]
    y_nse = cluster_abundance(nse.abundances, cluster_species)
    qse = solve_qse(net, t9, rho, ye, clusters=[QSECluster(cluster_species, y_nse)])
    assert qse.success
    assert qse.lambdas[0] == pytest.approx(0.0, abs=1e-6)
    for name in cluster_species:
        assert qse.abundances[name] == pytest.approx(nse.abundances[name], rel=1e-5)


# --- Cluster flows (compute_Ycdot) ----------------------------------------------

def test_cluster_ydot_counts_boundary_crossings_only():
    net = Network()
    for name in ["he4", "be8", "c12"]:
        net.add_species(Species.parse(name))
    net.reactions.add(Reaction.from_names(["he4", "he4"], ["be8"], constant_rate=2.0))
    net.reactions.add(Reaction.from_names(["be8", "he4"], ["c12"], constant_rate=3.0))
    abund = {"he4": 0.1, "be8": 0.01, "c12": 0.001}
    contributions = cluster_ydot(net, abund, ["be8", "c12"], t9=1.0, rho=1.0)
    # 2 he4 -> be8 enters the cluster (+1); be8 + he4 -> c12 stays inside (0)
    assert set(contributions) == {"he4 + he4 -> be8"}
    expected = 1 * 2.0 * 1.0 * 0.1 ** 2 / 2.0  # delta * rate * rho^1 * Y^2 / 2!
    assert contributions["he4 + he4 -> be8"] == pytest.approx(expected)


# --- Entropy generation with reverse flows (flow_utilities.cpp) ------------------

def test_entropy_generation_vanishes_at_nse_and_positive_off():
    net = Network()
    for name, me in [("he4", 2.4249), ("be8", 4.9416), ("c12", 0.0)]:
        net.add_species(Species.parse(name, mass_excess=me))
    net.reactions.add(Reaction.from_names(["he4", "he4"], ["be8"],
                                          rate_fits=[RateFit([2, 0, 0, 0, 0, 0, 0])]))
    net.reactions.add(Reaction.from_names(["be8", "he4"], ["c12"],
                                          rate_fits=[RateFit([3, 0, 0, 0, 0, 0, 0])]))
    t9, rho = 4.0, 1e6
    nse = solve_nse(net, t9, rho, ye=0.5)
    net.add_zone(Zone(abundances=dict(nse.abundances)))
    at_nse = entropy_generation_rate(net, 0, t9=t9, rho=rho)
    scale = sum(r.flux(nse.abundances, t9=t9, rho=rho) for r in net.reactions.reactions)
    assert abs(at_nse) < 1e-6 * scale
    net.add_zone(Zone(abundances={"he4": 0.25, "be8": 1e-12, "c12": 1e-12}))
    assert entropy_generation_rate(net, 1, t9=t9, rho=rho) > 0.0


# --- Integrated currents (update_flow_currents) -----------------------------------

def test_integrated_currents_match_abundance_change():
    # single irreversible decay: the integrated current must equal the change
    # in daughter abundance.
    net = Network()
    for name in ["ni56", "co56"]:
        net.add_species(Species.parse(name))
    lam = 3.0
    net.reactions.add(Reaction.from_names(["ni56"], ["co56"], constant_rate=lam))
    z = Zone(abundances={"ni56": 1e-2, "co56": 0.0})
    net.add_zone(z)
    thermo = constant_thermo(1.0, 1.0)
    times = time_grid(0.0, 0.3, 601)
    res = evolve_zone(net, z, times, thermo=thermo, method="rk4")
    currents = integrated_currents(net, res, thermo, use_reverse=False)
    dy_co56 = res.final_abundances["co56"] - 0.0
    assert currents["ni56 -> co56"] == pytest.approx(dy_co56, rel=1e-6)
