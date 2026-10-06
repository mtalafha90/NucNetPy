"""The electron gas, the physical entropy and its inversion, and the automatic
electron term of SkyNet screening.

The checks use limits with closed forms (classical, degenerate and
pair-dominated gases), thermodynamic identities that hold everywhere, and
constants typed in independently of the package.
"""
import math

import numpy as np
import pytest

from nucnetpy import Species, electrons, thermo
from nucnetpy.screening import SkyNetScreening, weak_screening_factor

# CODATA 2018, cgs, typed in independently of nucnetpy.constants
K_B = 1.380649e-16
HBAR = 1.054571817e-27
C = 2.99792458e10
M_E = 9.1093837015e-28
M_U = 1.66053906660e-24
N_A = 6.02214076e23
A_RAD = 7.565733250e-15


def _beta(t9):
    return K_B * t9 * 1e9 / (M_E * C * C)


# ---------------------------------------------------------------------------
# Electron gas
# ---------------------------------------------------------------------------

def test_quadrature_matches_adaptive_integration():
    from scipy.integrate import quad
    from scipy.special import expit

    for eta in (-300.0, -30.0, -5.0, -1.0, 0.0, 0.5, 3.0, 10.0, 40.0, 200.0, 1000.0):
        for beta in (1e-4, 1e-2, 0.17, 1.0, 10.0):
            log_scale, ours = electrons._integrals(eta, beta)
            ours = ours * math.exp(log_scale)
            rel = lambda x: math.sqrt(1.0 + 0.5 * beta * x)
            f = lambda x: expit(eta - x)
            df = lambda x: expit(eta - x) * expit(x - eta)
            kw = dict(epsabs=0.0, epsrel=1e-13, limit=500,
                      points=[eta] if eta > 0 else None)
            top = max(eta, 0.0) + 80.0
            ref = [quad(lambda x: math.sqrt(x) * rel(x) * (1 + beta * x) * f(x), 0, top, **kw)[0],
                   quad(lambda x: math.sqrt(x) * rel(x) * (1 + beta * x) * df(x), 0, top, **kw)[0],
                   quad(lambda x: x ** 1.5 * rel(x) ** 3 * f(x), 0, top, **kw)[0],
                   quad(lambda x: x ** 1.5 * rel(x) * (1 + beta * x) * f(x), 0, top, **kw)[0]]
            np.testing.assert_allclose(ours, ref, rtol=1e-10, atol=0.0)


def test_quadrature_matches_the_sommerfeld_expansion_when_degenerate():
    from scipy.integrate import quad

    for beta in (1e-4, 1.0):
        for eta in (1e4, 1e5):
            log_scale, ours = electrons._integrals(eta, beta)
            g = lambda x: math.sqrt(x) * math.sqrt(1 + 0.5 * beta * x) * (1 + beta * x)
            h = 1e-3 * eta
            sommerfeld = (quad(g, 0, eta, epsabs=0.0, epsrel=1e-13, limit=200)[0]
                          + math.pi ** 2 / 6.0 * (g(eta + h) - g(eta - h)) / (2 * h))
            assert ours[0] == pytest.approx(sommerfeld, rel=1e-10)


@pytest.mark.parametrize("t9", [1e-3, 1e-2, 0.1, 1.0, 5.0, 30.0, 100.0])
@pytest.mark.parametrize("rho", [1e-12, 1e-4, 1.0, 1e3, 1e8, 1e12, 1e14])
def test_chemical_potential_reproduces_the_net_electron_density(t9, rho):
    for ye in (0.05, 0.5, 1.0):
        xi, beta = electrons._solve_xi(t9, rho, ye)
        net, _ = electrons._net_density(xi, beta)
        assert net == pytest.approx(rho * N_A * ye, rel=1e-10)


def test_classical_limit():
    for t9, rho in [(0.0005, 1e-8), (0.005, 1e-5)]:
        gas = electrons.electron_gas(t9, rho, 0.5)
        n = rho * N_A * 0.5
        temperature = t9 * 1e9
        n_quantum = (M_E * K_B * temperature / (2 * math.pi * HBAR ** 2)) ** 1.5
        sackur_tetrode = 2.5 + math.log(2 * n_quantum / n)
        # The corrections are relativistic, of order beta.
        assert gas.pressure / (n * K_B * temperature) == pytest.approx(1.0, abs=3 * gas.beta)
        assert gas.entropy_density / n == pytest.approx(sackur_tetrode, rel=3 * gas.beta)
        assert gas.screening_term == pytest.approx(0.5, rel=1e-4)
        assert gas.n_positrons == 0.0


def test_degenerate_limit():
    for t9, rho in [(0.001, 1e5), (0.01, 1e9)]:
        gas = electrons.electron_gas(t9, rho, 0.5)
        n = rho * N_A * 0.5
        temperature = t9 * 1e9
        x_f = HBAR * (3 * math.pi ** 2 * n) ** (1 / 3) / (M_E * C)
        # Relativistic Sommerfeld result for the entropy per electron.
        expected = (math.pi ** 2 * K_B * temperature * math.sqrt(1 + x_f ** 2)
                    / (x_f ** 2 * M_E * C * C))
        assert gas.entropy_density / n == pytest.approx(expected, rel=1e-3)
        # Degenerate electrons hardly screen.
        assert gas.screening_term < 0.01 * 0.5
        assert gas.degeneracy_factor < 0.01


def test_hot_pair_plasma_has_seven_quarters_of_the_photon_entropy():
    gas = electrons.electron_gas(2000.0, 1e-3, 0.5)
    temperature = 2000.0e9
    photons = 4 * A_RAD * temperature ** 3 / (3 * K_B)
    assert gas.entropy_density / photons == pytest.approx(7.0 / 4.0, rel=1e-5)
    assert gas.screening_term > 1e6 * 0.5          # pairs dominate the screening


@pytest.mark.parametrize("t9, rho", [(0.01, 1e9), (5.0, 1e8), (3.0, 1.0), (1.0, 1e5)])
def test_entropy_and_density_follow_from_the_pressure(t9, rho):
    """s = dP/dT at fixed mu, and n = dP/dmu at fixed T, with mu the full
    chemical potential (electrons mu, positrons -mu)."""
    gas = electrons.electron_gas(t9, rho, 0.5)
    mu = gas.mu_mev / electrons.ELECTRON_MASS_MEV          # in m_e c^2

    def pressure(t9_, mu_):
        beta = _beta(t9_)
        return (electrons._species((mu_ - 1.0) / beta, beta)[2]
                + electrons._species((-mu_ - 1.0) / beta, beta)[2])

    # Five-point differences: degenerate matter has a small thermal pressure
    # on top of a large cold one, so the step cannot be very small.
    h = 1e-3 * t9
    dp_dt = (8 * (pressure(t9 + h, mu) - pressure(t9 - h, mu))
             - (pressure(t9 + 2 * h, mu) - pressure(t9 - 2 * h, mu))) / (12 * h * 1e9)
    assert dp_dt / K_B == pytest.approx(gas.entropy_density, rel=1e-6)

    if gas.degeneracy_factor < 10.0:
        # Where pairs dominate, n- - n+ is a tiny difference of the pressure
        # and cannot be resolved this way; the net density is tested directly
        # in test_chemical_potential_reproduces_the_net_electron_density.
        dmu = 1e-4 * abs(mu)
        dp_dmu = (8 * (pressure(t9, mu + dmu) - pressure(t9, mu - dmu))
                  - (pressure(t9, mu + 2 * dmu) - pressure(t9, mu - 2 * dmu))) / (12 * dmu * M_E * C * C)
        assert dp_dmu == pytest.approx(gas.n_electrons - gas.n_positrons, rel=1e-7)

    def net(eta):
        beta = gas.beta
        return electrons._species(eta, beta)[0] - electrons._species(-eta - 2 / beta, beta)[0]

    d_eta = 1e-6 * max(1.0, abs(gas.eta))
    assert (net(gas.eta + d_eta) - net(gas.eta - d_eta)) / (2 * d_eta) == pytest.approx(
        gas.dn_deta, rel=1e-6)


def test_invalid_state_is_rejected():
    with pytest.raises(ValueError):
        electrons.electron_gas(0.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        electrons.electron_gas(1.0, 1.0, 0.0)


# ---------------------------------------------------------------------------
# Entropy
# ---------------------------------------------------------------------------

HE4 = {"he4": Species("he4", 2, 4, spin=0.0)}


def test_ion_entropy_matches_sackur_tetrode():
    t9, rho, y = 5.0, 1e8, 0.25
    temperature = t9 * 1e9
    n_quantum = (4 * M_U * K_B * temperature / (2 * math.pi * HBAR ** 2)) ** 1.5
    expected = y * (2.5 + math.log(n_quantum / (rho * N_A * y)))
    assert thermo.ion_entropy_per_nucleon({"he4": y}, HE4, t9, rho) == pytest.approx(expected, rel=1e-9)
    assert expected == pytest.approx(3.243, abs=1e-3)


def test_photon_entropy():
    # The textbook figure: 1.213 T9^3 / (rho / 1e5 g cm^-3).
    assert thermo.photon_entropy_per_nucleon(1.0, 1e5) == pytest.approx(1.2133, rel=1e-4)
    assert thermo.photon_entropy_per_nucleon(5.0, 1e8) == pytest.approx(
        4 * A_RAD * 5e9 ** 3 / (3 * 1e8 * N_A * K_B), rel=1e-12)


def test_excited_states_enter_through_the_partition_function():
    """The ion entropy must equal -d(free energy)/dT, including d ln G / dT."""
    from nucnetpy.nse import _log_quantum_abundance

    fe56 = Species("fe56", 26, 56, spin=0.0,
                   partition={1.0: 1.0, 4.0: 1.2, 6.0: 2.0, 8.0: 4.5})
    species = {"fe56": fe56}
    y, rho = 1.0 / 56.0, 1e8

    def free_energy_over_k(t9):          # per nucleon, divided by k, T in GK
        log_q = _log_quantum_abundance(fe56, t9, rho)
        return -t9 * 1e9 * y * (1.0 + log_q - math.log(y))

    for t9 in (2.5, 5.0, 7.0):
        h = 1e-6
        s_fd = -(free_energy_over_k(t9 + h) - free_energy_over_k(t9 - h)) / (2 * h * 1e9)
        assert thermo.ion_entropy_per_nucleon({"fe56": y}, species, t9, rho) == pytest.approx(s_fd, rel=1e-7)
    # Without excited states the entropy is lower.
    assert (thermo.ion_entropy_per_nucleon({"fe56": y}, species, 7.0, rho, include_partition=False)
            < thermo.ion_entropy_per_nucleon({"fe56": y}, species, 7.0, rho))


def test_entropy_per_nucleon_adds_the_parts():
    s = thermo.entropy_per_nucleon({"he4": 0.25, "gamma": 1.0, "c12": 0.0}, HE4, 5.0, 1e8)
    assert s.ions == pytest.approx(3.24296, rel=1e-5)
    assert s.electrons == pytest.approx(electrons.electron_gas(5.0, 1e8, 0.5).entropy_per_nucleon)
    assert s.photons == pytest.approx(thermo.photon_entropy_per_nucleon(5.0, 1e8))
    assert s.coulomb == 0.0
    assert s.total == pytest.approx(s.ions + s.electrons + s.photons)

    with_coulomb = thermo.entropy_per_nucleon({"he4": 0.25}, HE4, 5.0, 1e8, coulomb=True)
    from nucnetpy.coulomb import coulomb_entropy_per_nucleon
    assert with_coulomb.coulomb == pytest.approx(
        coulomb_entropy_per_nucleon({"he4": 0.25}, HE4, 5.0, 1e8, 0.5))

    ions_only = thermo.entropy_per_nucleon({"he4": 0.25}, HE4, 5.0, 1e8,
                                           electrons=False, photons=False)
    assert ions_only.total == pytest.approx(s.ions)


def test_entropy_accepts_a_network():
    from nucnetpy.core import Network
    network = Network(species=dict(HE4))
    assert thermo.entropy_per_nucleon({"he4": 0.25}, network, 5.0, 1e8).total == pytest.approx(
        thermo.entropy_per_nucleon({"he4": 0.25}, HE4, 5.0, 1e8).total)


@pytest.mark.parametrize("t9, rho", [(5.0, 1e-2), (5.0, 1e3), (5.0, 1e8), (0.5, 1e12), (0.05, 1e6)])
@pytest.mark.parametrize("coulomb", [False, True])
def test_density_and_temperature_inversions_round_trip(t9, rho, coulomb):
    composition = {"he4": 0.2, "c12": 0.1 / 12.0, "n": 1e-3}
    species = {"he4": Species("he4", 2, 4, spin=0.0), "c12": Species("c12", 6, 12, spin=0.0),
               "n": Species("n", 0, 1, spin=0.5)}
    s = thermo.entropy_per_nucleon(composition, species, t9, rho, coulomb=coulomb).total
    assert thermo.density_for_entropy(s, composition, species, t9, coulomb=coulomb) == pytest.approx(rho, rel=1e-10)
    assert thermo.t9_for_entropy(s, composition, species, rho, coulomb=coulomb) == pytest.approx(t9, rel=1e-10)


def test_inversion_reports_an_unreachable_entropy():
    with pytest.raises(ValueError, match="outside the range"):
        thermo.t9_for_entropy(1e30, {"he4": 0.25}, HE4, 1e8)


def test_old_proxies_warn():
    for call in (lambda: thermo.entropy_ideal_ions(1e8, 5e9, 0.25),
                 lambda: thermo.density_from_entropy(4.0, 5e9, 0.25),
                 lambda: thermo.temperature_from_entropy(4.0, 1e8, 0.25)):
        with pytest.warns(DeprecationWarning):
            call()


def test_constant_entropy_expansion():
    from nucnetpy import evolve_zone, time_grid
    from nucnetpy.core import Zone
    from nucnetpy.io.xml import read_xml

    network = read_xml("tests/golden/golden_network.xml")
    zone = Zone(abundances={"he4": 0.2, "c12": 0.1 / 12.0})
    # Conditions where little burns, so the expansion sets the temperature.
    # (Burning helium to carbon at fixed entropy would heat the matter: fewer
    # particles must share the same entropy.)
    entropy = thermo.entropy_per_nucleon(zone.abundances, network, 0.5, 1e4).total
    density = lambda t: 1e4 * math.exp(-t / 0.5)
    thermo_fn = thermo.constant_entropy_thermo(entropy, density, network)

    t9_0, rho_0 = thermo_fn(0.0, zone.abundances)
    assert (t9_0, rho_0) == (pytest.approx(0.5, rel=1e-10), 1e4)

    result = evolve_zone(network, zone, time_grid(0.0, 1.0, 6), thermo=thermo_fn)
    assert result.success
    previous = math.inf
    for t, y in zip(result.time, result.y):
        abundances = dict(zip(result.species, y))
        t9, rho = thermo_fn(float(t), abundances)
        assert thermo.entropy_per_nucleon(abundances, network, t9, rho).total == pytest.approx(entropy, rel=1e-10)
        assert t9 < previous                        # the matter cools as it expands
        previous = t9


# ---------------------------------------------------------------------------
# SkyNet screening with the electron term computed automatically
# ---------------------------------------------------------------------------

SCREEN_SPECIES = {"he4": Species("he4", 2, 4), "c12": Species("c12", 6, 12),
                  "o16": Species("o16", 8, 16)}


def _exponent(screening, composition, t9, rho):
    from nucnetpy.reactions import Reaction
    screening.update(composition, t9=t9, rho=rho)
    return math.log(screening.factor(Reaction.from_names(["c12", "he4"], ["o16"])))


def test_automatic_electron_term_recovers_salpeter_for_non_degenerate_electrons():
    composition = {"he4": 0.25}
    auto = SkyNetScreening(SCREEN_SPECIES)
    for t9, rho in [(0.05, 1e-2), (0.1, 1.0)]:
        gas = electrons.electron_gas(t9, rho, 0.5)
        assert gas.screening_term == pytest.approx(0.5, rel=1e-3)
        assert _exponent(auto, composition, t9, rho) == pytest.approx(
            math.log(weak_screening_factor(6, 2, t9, rho, 0.5, 1.0)), rel=1e-3)


def test_automatic_electron_term_vanishes_for_degenerate_electrons():
    composition = {"c12": 1.0 / 12.0}
    auto = SkyNetScreening(SCREEN_SPECIES)
    ions_only = SkyNetScreening(SCREEN_SPECIES, pair_term=0.0)
    auto.update(composition, t9=0.01, rho=1e9)
    assert auto.electron_term(0.01, 1e9, 0.5) < 0.01 * 0.5
    assert _exponent(auto, composition, 0.01, 1e9) == pytest.approx(
        _exponent(ions_only, composition, 0.01, 1e9), rel=1e-2)


def test_electron_term_changes_smoothly_and_override_still_works():
    auto = SkyNetScreening(SCREEN_SPECIES)
    terms = [auto.electron_term(0.1, rho, 0.5) for rho in np.geomspace(1.0, 1e10, 41)]
    assert all(b < a for a, b in zip(terms, terms[1:]))       # more degenerate, less screening
    assert terms[0] == pytest.approx(0.5, rel=1e-3) and terms[-1] < 0.05
    assert max(b / a for a, b in zip(terms, terms[1:])) > 0.7  # no jumps

    assert SkyNetScreening(SCREEN_SPECIES, pair_term=0.25).electron_term(0.1, 1.0, 0.5) == 0.25
    assert auto.electron_term(3.0, 1e2, 0.5) > 0.5            # positron pairs add screening


def test_electron_state_is_cached():
    auto = SkyNetScreening(SCREEN_SPECIES)
    calls = []
    original = electrons.electron_gas

    def counting(*args):
        calls.append(args)
        return original(*args)

    electrons.electron_gas = counting
    try:
        for _ in range(5):
            auto.update({"he4": 0.25}, t9=1.0, rho=1e6)
        auto.update({"he4": 0.25 * (1 + 1e-14)}, t9=1.0, rho=1e6)
        assert len(calls) == 1
        auto.update({"he4": 0.25}, t9=1.1, rho=1e6)
        assert len(calls) == 2
    finally:
        electrons.electron_gas = original
