"""Properties of the shared NSE/QSE solver that must not depend on SciPy's version.

Each test guards a defect that once made results differ between machines:
non-unique chemical potentials reported as if they were unique, a success flag
looser than the requested tolerance, and trace abundances left undetermined by
a charge residual below floating-point resolution.
"""
import math

import pytest

from nucnetpy import Network, Species, solve_nse, solve_qse, QSECluster, cluster_abundance

ALPHA_CHAIN = [("he4", 2.4249), ("c12", 0.0), ("o16", -4.7370), ("ne20", -7.0419),
               ("mg24", -13.9336), ("si28", -21.4928), ("s32", -26.0157),
               ("ar36", -30.2315), ("ca40", -34.8463), ("ti44", -37.5485),
               ("cr48", -42.8215), ("fe52", -48.3320), ("ni56", -53.9042)]
NUCLEONS = [("n", 8.0713), ("h1", 7.2890)]


def network(entries):
    net = Network()
    for name, mass_excess in entries:
        net.add_species(Species.parse(name, mass_excess=mass_excess, spin=0.0))
    return net


def test_potentials_are_minimum_norm_when_not_unique():
    # With only N = Z nuclei, every abundance depends on mu_p + mu_n alone.
    # The split is arbitrary, so the solver must report a reproducible one.
    res = solve_nse(network(ALPHA_CHAIN), t9=5.0, rho=1e8, ye=0.5)
    assert res.success
    assert res.mu_p == pytest.approx(res.mu_n, abs=1e-12)
    assert "not unique" in res.message


def test_potentials_are_left_alone_when_unique():
    res = solve_nse(network(ALPHA_CHAIN + NUCLEONS), t9=5.0, rho=1e8, ye=0.5)
    assert res.success
    assert "not unique" not in res.message
    assert abs(res.mu_p - res.mu_n) > 1e-3


@pytest.mark.parametrize("t9, rho", [(2.0, 1e9), (3.5, 1e9), (5.0, 1e8)])
def test_trace_nucleons_balance_charge_at_ye_one_half(t9, rho):
    # At Ye = 0.5 the bulk is in N = Z nuclei and the only charged-asymmetric
    # species are n and p, so exact charge balance requires Y(n) = Y(p).
    # Their effect on sum(Z Y) - Ye can be ~1e-18, below double precision, so
    # a linear charge residual leaves them undetermined; this pins the fix.
    res = solve_nse(network(ALPHA_CHAIN + NUCLEONS), t9=t9, rho=rho, ye=0.5)
    assert res.success
    assert res.abundances["n"] / res.abundances["h1"] == pytest.approx(1.0, rel=1e-9)


def test_success_means_the_requested_tolerance():
    # A constraint set to the unconstrained NSE value must give lambda = 0 to
    # far better than the old 1e-6 acceptance threshold.
    net = network(NUCLEONS + [("he4", 2.4249), ("si28", -21.4928),
                              ("fe56", -60.6054), ("ni56", -53.9042)])
    nse = solve_nse(net, t9=5.0, rho=1e8, ye=0.5)
    heavy = ["si28", "fe56", "ni56"]
    qse = solve_qse(net, t9=5.0, rho=1e8, ye=0.5,
                    clusters=[QSECluster(heavy, cluster_abundance(nse.abundances, heavy))])
    assert qse.success
    assert abs(qse.lambdas[0]) < 1e-9


def test_unreachable_ye_is_reported_as_failure():
    # Ye = 0.55 cannot be reached without free protons.
    res = solve_nse(network(ALPHA_CHAIN + [("n", 8.0713)]), t9=5.0, rho=1e8, ye=0.55)
    assert not res.success


def test_qse_rejects_a_cluster_with_no_solvable_members():
    with pytest.raises(ValueError, match="none of the species"):
        solve_qse(network(ALPHA_CHAIN), t9=5.0, rho=1e8, ye=0.5,
                  clusters=[QSECluster(["u238"], 1e-3)])
