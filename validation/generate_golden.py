#!/usr/bin/env python
"""Generate golden-output files for the numerical-identity test suite.

The files written here are consumed by ``tests/test_golden_identity.py``.
They freeze the numerical output of the current nucnetpy implementation for:

1. ReacLib/tabular/constant reaction rates on a T9 grid,
2. electron-screening factors on a (T9, rho, Ye) grid,
3. two-dimensional weak-rate table interpolation,
4. one-zone RHS ``ydot`` at fixed thermodynamic states,
5. a full fixed-step trajectory on an exact time grid.

(The first stage of the recommended validation order, XML round-tripping, is a
pure identity test and needs no golden data.)

By default the goldens are *self-consistent regression snapshots*: they detect
any future change in nucnetpy's numerics.  To validate against an original C++
NucNet Tools build instead, overwrite the ``data`` blocks of these JSON files
with values produced by the original tools for the same inputs (same fixture
network, same grids), set ``source`` to a label describing that build, and
loosen the per-file ``rtol``/``atol`` to the agreement you require.

Usage::

    python validation/generate_golden.py                 # refresh snapshots
    python validation/generate_golden.py --source my-cpp-build --rtol 1e-6
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from nucnetpy import read_xml, read_weak_table, evolve_zone, constant_thermo
from nucnetpy.screening import weak_screening_factor, graboske_intermediate_factor

REPO_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"
FIXTURE_XML = GOLDEN_DIR / "golden_network.xml"
FIXTURE_WEAK = GOLDEN_DIR / "golden_weak_rates.txt"

T9_GRID = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
SCREEN_GRID = [
    # (z1, z2, t9, rho, ye)
    (1, 6, 0.05, 1.0e2, 0.5),
    (1, 6, 0.5, 1.0e5, 0.5),
    (2, 6, 0.5, 1.0e5, 0.5),
    (2, 8, 1.0, 1.0e6, 0.5),
    (6, 6, 2.0, 1.0e7, 0.5),
    (8, 10, 3.0, 1.0e8, 0.45),
]
WEAK_POINTS = [
    # on-grid, interior, and clamped-outside points
    (0.2, 1.0e3), (1.0, 1.0e5), (5.0, 1.0e7),
    (0.6, 3.0e4), (2.5, 5.0e5), (4.0, 2.0e6),
    (0.05, 1.0e2), (9.0, 1.0e8),
]
YDOT_STATES = [
    # (t9, rho) evaluated with the fixture zone's abundances
    (0.2, 1.0e4), (0.5, 1.0e5), (1.0, 1.0e6), (3.0, 1.0e7),
]
TRAJ_T0, TRAJ_T1, TRAJ_STEPS = 0.0, 2.0e-5, 201
TRAJ_SAMPLE_EVERY = 20


def _write(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {path.relative_to(REPO_ROOT)}")


def generate(source: str, rtol: float, atol: float) -> None:
    net = read_xml(FIXTURE_XML)
    common = {"source": source, "fixture": FIXTURE_XML.name, "rtol": rtol, "atol": atol}

    # 2. rate-by-rate ReacLib/tabular/constant comparisons on a T9 grid
    rates = {
        r.label or r.string: {f"{t9:g}": r.rate(t9) for t9 in T9_GRID}
        for r in net.reactions.reactions
    }
    _write(GOLDEN_DIR / "rates_reaclib.json", {
        **common,
        "description": "reaction.rate(t9) per reaction label on the T9 grid",
        "t9_grid": T9_GRID,
        "data": rates,
    })

    # 3. screening-factor comparisons
    screening = [
        {
            "z1": z1, "z2": z2, "t9": t9, "rho": rho, "ye": ye,
            "weak": weak_screening_factor(z1, z2, t9, rho, ye),
            "graboske": graboske_intermediate_factor(z1, z2, t9, rho, ye),
        }
        for (z1, z2, t9, rho, ye) in SCREEN_GRID
    ]
    _write(GOLDEN_DIR / "screening.json", {
        **common,
        "description": "weak and graboske screening factors on a (Z1,Z2,T9,rho,Ye) grid",
        "data": screening,
    })

    # 4. weak-rate table comparisons
    table = read_weak_table(FIXTURE_WEAK)
    weak = [
        {
            "t9": t9, "rho_ye": rho_ye,
            "rate": table.rate(t9, rho_ye),
            "nu_loss": table.neutrino_loss(t9, rho_ye),
        }
        for (t9, rho_ye) in WEAK_POINTS
    ]
    _write(GOLDEN_DIR / "weak_rates.json", {
        **common,
        "description": "2-D weak-rate table interpolation at on/off/outside-grid points",
        "table": FIXTURE_WEAK.name,
        "data": weak,
    })

    # 5. one-zone RHS ydot at fixed states
    zone = net.zone(0)
    ydots = {
        f"{t9:g}|{rho:g}": {
            name: val
            for name, val in sorted(net.reactions.ydot(zone.abundances, t9=t9, rho=rho).items())
        }
        for (t9, rho) in YDOT_STATES
    }
    _write(GOLDEN_DIR / "ydot.json", {
        **common,
        "description": "network ydot for the fixture zone abundances at fixed (T9, rho)",
        "abundances": dict(sorted(zone.abundances.items())),
        "data": ydots,
    })

    # 6. full trajectory on an identical time grid (deterministic fixed-step rk4)
    times = np.linspace(TRAJ_T0, TRAJ_T1, TRAJ_STEPS)
    res = evolve_zone(net, zone, times,
                      thermo=constant_thermo(zone.temperature9(), zone.density()),
                      method="rk4")
    sample = list(range(0, TRAJ_STEPS, TRAJ_SAMPLE_EVERY))
    if sample[-1] != TRAJ_STEPS - 1:
        sample.append(TRAJ_STEPS - 1)
    _write(GOLDEN_DIR / "trajectory.json", {
        **common,
        "description": "fixed-step rk4 one-zone evolution sampled on the exact time grid",
        "method": "rk4",
        "t0": TRAJ_T0, "t1": TRAJ_T1, "steps": TRAJ_STEPS,
        "t9": zone.temperature9(), "rho": zone.density(),
        "species": res.species,
        "data": {
            "time": [float(res.time[i]) for i in sample],
            "y": [[float(v) for v in res.y[i]] for i in sample],
        },
    })


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source", default="nucnetpy-python",
                    help="label recorded in each golden file (e.g. a C++ build id)")
    ap.add_argument("--rtol", type=float, default=1e-10,
                    help="relative tolerance the test suite applies to these goldens")
    ap.add_argument("--atol", type=float, default=1e-300,
                    help="absolute tolerance the test suite applies to these goldens")
    args = ap.parse_args(argv)
    generate(args.source, args.rtol, args.atol)


if __name__ == "__main__":
    main()
