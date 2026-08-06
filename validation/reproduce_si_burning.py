"""Reproduce the silicon-burning result of the accompanying article.

This script needs nothing but the package and the files in
validation/reference/. The third-party JINA/libnucnet database is not
required: the 15-nuclide, 32-reaction subnetwork the result depends on is
archived alongside its initial zone, and validation/reference/
si_burning_provenance.json records the checksums of the database it came from.

    python validation/reproduce_si_burning.py

The exit status is 0 if every archived value is reproduced within the
tolerances recorded in si_burning_expected.json, and 1 otherwise.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from nucnetpy.io.xml import read_xml
from nucnetpy.nse import solve_nse
from nucnetpy.solver import constant_thermo, evolve_zone, time_grid

REFERENCE_DIR = Path(__file__).resolve().parent / "reference"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def compare(label, got, want, rtol, atol, failures):
    ok = bool(np.isclose(got, want, rtol=rtol, atol=atol))
    if not ok:
        failures.append(f"{label}: got {got!r}, expected {want!r}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference-dir", default=str(REFERENCE_DIR))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    ref = Path(args.reference_dir)
    spec = json.loads((ref / "si_burning_expected.json").read_text())
    network_path = ref / spec["network"]["file"]

    failures: list[str] = []

    # The archived network is the input; confirm it is the one the expectations
    # were produced from before trusting any comparison against them.
    digest = sha256(network_path)
    if digest != spec["network"]["sha256"]:
        failures.append(f"network checksum: got {digest}, "
                        f"expected {spec['network']['sha256']}")

    net = read_xml(str(network_path))
    if len(net.species) != spec["network"]["species"]:
        failures.append(f"species: {len(net.species)} != {spec['network']['species']}")
    if len(net.reactions.reactions) != spec["network"]["reactions"]:
        failures.append(f"reactions: {len(net.reactions.reactions)} != "
                        f"{spec['network']['reactions']}")
    if not net.zones:
        failures.append("the archived network carries no initial zone")
        print("\n".join(failures), file=sys.stderr)
        return 1

    c = spec["conditions"]
    tol = spec["tolerances"]
    exp = spec["expected"]
    zone = net.zone(0)
    ye = zone.ye(net.species)

    result = evolve_zone(net, zone, time_grid(0.0, c["t_end"], c["steps"]),
                         thermo=constant_thermo(c["t9"], c["rho"]),
                         method=c["method"], rtol=c["rtol"], atol=c["atol"])
    if not result.success:
        failures.append(f"evolution failed: {result.message}")

    nse = solve_nse(net, t9=c["t9"], rho=c["rho"], ye=ye)
    if not nse.success:
        failures.append("NSE solve failed")

    final = result.final_abundances
    xsum = sum(net.species[k].a * v for k, v in final.items() if k in net.species)

    compare("initial Ye", ye, c["initial_ye"], tol["rtol"], tol["atol"], failures)
    compare("sum A_i Y_i", xsum, exp["final_xsum"], 0.0, tol["xsum_atol"], failures)
    compare("NSE mu_p", nse.mu_p, exp["nse_mu_p"], tol["rtol"], tol["atol"], failures)
    compare("NSE mu_n", nse.mu_n, exp["nse_mu_n"], tol["rtol"], tol["atol"], failures)

    if not args.quiet:
        print(f"network: {len(net.species)} species, "
              f"{len(net.reactions.reactions)} reactions, checksum verified")
        print(f"conditions: T9={c['t9']}, rho={c['rho']:g} g/cm^3, "
              f"t_end={c['t_end']} s, {c['method']}, "
              f"rtol={c['rtol']:g}, atol={c['atol']:g}")
        print(f"\n{'species':>8}  {'Y (network)':>16}  {'Y (expected)':>16}  "
              f"{'X (network)':>14}")

    for name, want in sorted(exp["final_abundances"].items(), key=lambda kv: -kv[1]):
        got = float(final.get(name, 0.0))
        ok = compare(f"final Y({name})", got, want, tol["rtol"], tol["atol"], failures)
        if not args.quiet and (want > 1e-12 or not ok):
            a = net.species[name].a if name in net.species else 0
            flag = "" if ok else "   <-- MISMATCH"
            print(f"{name:>8}  {got:16.8e}  {want:16.8e}  {a*got:14.6e}{flag}")

    for name, want in exp["nse_abundances"].items():
        compare(f"NSE Y({name})", float(nse.abundances.get(name, 0.0)), want,
                tol["rtol"], tol["atol"], failures)

    if failures:
        print(f"\nFAILED: {len(failures)} comparison(s) did not reproduce",
              file=sys.stderr)
        for f in failures[:20]:
            print(f"  {f}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(f"\nOK: every archived value reproduced within "
              f"rtol={tol['rtol']:g}, atol={tol['atol']:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
