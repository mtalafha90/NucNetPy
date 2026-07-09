"""Golden-output numerical-identity tests.

These tests implement the validation order recommended in
``docs/PURE_PYTHON_PORT_STATUS.md``:

1. XML round-trip tests,
2. rate-by-rate ReacLib comparisons on a T9 grid,
3. screening-factor comparisons,
4. weak-rate table comparisons,
5. one-zone RHS ``ydot`` comparisons at fixed state,
6. full trajectory comparisons on identical time grids.

The golden files in ``tests/golden/`` are snapshots frozen by
``validation/generate_golden.py``.  Out of the box they pin nucnetpy against
itself, so any change to the numerics fails loudly here.  To validate against
an original C++ NucNet Tools build, regenerate the golden ``data`` blocks from
that build for the same fixture inputs; each file carries its own
``rtol``/``atol`` so the required agreement is data-driven, not hard-coded.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from nucnetpy import read_xml, write_xml, read_weak_table, evolve_zone, constant_thermo
from nucnetpy.screening import weak_screening_factor, graboske_intermediate_factor

GOLDEN_DIR = Path(__file__).parent / "golden"
FIXTURE_XML = GOLDEN_DIR / "golden_network.xml"


def load_golden(name: str) -> dict:
    return json.loads((GOLDEN_DIR / name).read_text())


def assert_close(actual: float, expected: float, rtol: float, atol: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(
            f"{label}: got {actual!r}, golden {expected!r} (rtol={rtol}, atol={atol})")


@pytest.fixture(scope="module")
def net():
    return read_xml(FIXTURE_XML)


# 1. XML round-trip -----------------------------------------------------------

def test_xml_round_trip(net, tmp_path):
    out = tmp_path / "roundtrip.xml"
    write_xml(net, out)
    net2 = read_xml(out)

    assert set(net2.species) == set(net.species)
    for name, sp in net.species.items():
        sp2 = net2.species[name]
        assert (sp2.z, sp2.a) == (sp.z, sp.a), name
        assert sp2.mass_excess == sp.mass_excess, name
        assert sp2.partition == sp.partition, name

    assert len(net2.reactions.reactions) == len(net.reactions.reactions)
    by_key = {r.key: r for r in net2.reactions.reactions}
    for r in net.reactions.reactions:
        r2 = by_key[r.key]
        assert r2.constant_rate == r.constant_rate, r.string
        assert len(r2.rate_fits) == len(r.rate_fits), r.string
        for f, f2 in zip(r.rate_fits, r2.rate_fits):
            assert list(f2.coefficients) == list(f.coefficients), r.string

    assert len(net2.zones) == len(net.zones)
    for z, z2 in zip(net.zones, net2.zones):
        assert z2.abundances == z.abundances
        merged = {**z.properties, **z.optional_properties}
        merged2 = {**z2.properties, **z2.optional_properties}
        assert merged2 == merged


def test_round_trip_is_a_fixed_point(net, tmp_path):
    # write -> read -> write must be byte-identical: the writer is canonical.
    p1, p2 = tmp_path / "a.xml", tmp_path / "b.xml"
    write_xml(net, p1)
    write_xml(read_xml(p1), p2)
    assert p1.read_bytes() == p2.read_bytes()


# 2. ReacLib rates on a T9 grid ----------------------------------------------

def test_reaclib_rates_on_t9_grid(net):
    g = load_golden("rates_reaclib.json")
    by_label = {r.label or r.string: r for r in net.reactions.reactions}
    assert set(by_label) == set(g["data"])
    for label, expected in g["data"].items():
        r = by_label[label]
        for t9_str, val in expected.items():
            assert_close(r.rate(float(t9_str)), val, g["rtol"], g["atol"],
                         f"{label} @ T9={t9_str}")


# 3. Screening factors ---------------------------------------------------------

def test_screening_factors():
    g = load_golden("screening.json")
    for row in g["data"]:
        args = (row["z1"], row["z2"], row["t9"], row["rho"], row["ye"])
        assert_close(weak_screening_factor(*args), row["weak"], g["rtol"], g["atol"],
                     f"weak screening {args}")
        assert_close(graboske_intermediate_factor(*args), row["graboske"],
                     g["rtol"], g["atol"], f"graboske screening {args}")
        assert row["weak"] >= 1.0 and row["graboske"] >= 1.0


# 4. Weak-rate table interpolation ---------------------------------------------

def test_weak_rate_table_interpolation():
    g = load_golden("weak_rates.json")
    table = read_weak_table(GOLDEN_DIR / g["table"])
    for row in g["data"]:
        assert_close(table.rate(row["t9"], row["rho_ye"]), row["rate"],
                     g["rtol"], g["atol"], f"rate @ ({row['t9']}, {row['rho_ye']})")
        assert_close(table.neutrino_loss(row["t9"], row["rho_ye"]), row["nu_loss"],
                     g["rtol"], g["atol"], f"nu_loss @ ({row['t9']}, {row['rho_ye']})")


# 5. One-zone RHS ydot at fixed state -------------------------------------------

def test_ydot_at_fixed_states(net):
    g = load_golden("ydot.json")
    zone = net.zone(0)
    assert dict(sorted(zone.abundances.items())) == g["abundances"]
    for state, expected in g["data"].items():
        t9, rho = (float(v) for v in state.split("|"))
        actual = net.reactions.ydot(zone.abundances, t9=t9, rho=rho)
        assert set(actual) == set(expected), state
        for name, val in expected.items():
            assert_close(actual[name], val, g["rtol"], g["atol"], f"ydot[{name}] @ {state}")
        # a valid RHS conserves nucleon number: sum(A_i * ydot_i) == 0
        drift = sum(net.species[n].a * v for n, v in actual.items())
        scale = max(abs(v) for v in actual.values()) or 1.0
        assert abs(drift) <= 1e-10 * scale, f"nucleon-number drift at {state}: {drift}"


# 6. Full trajectory on an identical time grid ----------------------------------

def test_full_trajectory(net):
    g = load_golden("trajectory.json")
    zone = net.zone(0)
    times = np.linspace(g["t0"], g["t1"], g["steps"])
    res = evolve_zone(net, zone, times,
                      thermo=constant_thermo(g["t9"], g["rho"]), method=g["method"])
    assert res.species == g["species"]
    idx = {float(t): k for k, t in enumerate(res.time)}
    for t, yrow in zip(g["data"]["time"], g["data"]["y"]):
        k = idx[float(t)]
        for j, expected in enumerate(yrow):
            assert_close(float(res.y[k, j]), expected, g["rtol"], g["atol"],
                         f"Y[{g['species'][j]}] @ t={t}")


def test_trajectory_solver_cross_check(net):
    # Not a golden comparison: the adaptive BDF result must agree with the
    # deterministic rk4 golden path to solver accuracy, whatever SciPy version.
    g = load_golden("trajectory.json")
    zone = net.zone(0)
    times = np.linspace(g["t0"], g["t1"], g["steps"])
    thermo = constant_thermo(g["t9"], g["rho"])
    rk4 = evolve_zone(net, zone, times, thermo=thermo, method="rk4")
    bdf = evolve_zone(net, zone, times, thermo=thermo, method="bdf", rtol=1e-9, atol=1e-18)
    assert bdf.success
    ya, yb = rk4.y[-1], bdf.y[-1]
    denom = np.maximum(np.abs(ya), 1e-12)
    assert float(np.max(np.abs(ya - yb) / denom)) < 1e-3
