"""Measure NucNetPy performance across a ladder of network sizes.

The script builds a sequence of progressively larger networks from a JINA /
libnucnet database by applying charge and mass-number cuts, then measures the
costs that dominate a reaction-network calculation: reading the database,
evaluating the right-hand side, forming the Jacobian, and integrating a
one-zone trajectory.

It writes a JSON record and a LaTeX table so that the performance table of the
accompanying manuscript can be regenerated from a single command.

Usage
-----
python validation/benchmark_performance.py \
    --nuclides /path/to/nuclides.xml \
    --reactions /path/to/reaction_data.xml \
    --zone /path/to/zone.xml \
    --json benchmark.json --latex benchmark.tex

Every measurement is taken on the machine that runs the script; no timing is
inferred or copied from another source.
"""
from __future__ import annotations

import argparse
import copy
import json
import platform
import resource
import statistics
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from nucnetpy.io.jina import read_jina_xml
from nucnetpy.network_limiter import limit_network, select_species
from nucnetpy.species import is_massless
from nucnetpy.solver import (analytic_jacobian, constant_thermo, evolve_zone,
                             jacobian, rhs, time_grid)

# (label, zmax, amax).  ``None`` keeps the complete database.
DEFAULT_LADDER = [
    (r"Light ($Z\leq2$)", 2, 4),
    (r"CNO ($Z\leq8$)", 8, 20),
    (r"NeNa ($Z\leq10$)", 10, 24),
    (r"MgAl ($Z\leq14$)", 14, 32),
    (r"Iron group ($Z\leq20$)", 20, 44),
    (r"Extended ($Z\leq30$)", 30, 70),
    ("Full database", None, None),
]


def peak_memory_mb() -> float:
    """Return peak resident set size in MiB (Linux reports KiB).

    This is ``ru_maxrss``, the high-water mark for the whole process.  It never
    decreases, and the complete database is loaded before any network is cut,
    so the value recorded against a small case is dominated by that load and
    says nothing about the case itself.  It is kept as a ceiling on the run,
    not as a per-case measurement.  For an isolated figure, run
    ``validation/reproduce_si_burning.py``, which starts from the archived
    15-species network in a clean process and reports its own increment.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _time_repeated(fn, min_repeats: int = 3, budget: float = 2.0) -> float:
    """Return the median wall-clock time of ``fn`` in seconds.

    Repeats until either ``min_repeats`` samples are collected or ``budget``
    seconds have elapsed, so that cheap calls are averaged over many samples
    while expensive ones are not repeated needlessly.
    """
    samples: List[float] = []
    start = time.perf_counter()
    while True:
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
        if len(samples) >= min_repeats and (time.perf_counter() - start) > budget:
            break
        if len(samples) >= 200:
            break
    return float(statistics.median(samples))


def environment_record() -> Dict[str, str]:
    import scipy
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "processor": platform.processor() or platform.machine(),
    }


def benchmark_case(base_net, label: str, zmax: Optional[int], amax: Optional[int],
                   t9: float, rho: float, t_end: float, steps: int,
                   rtol: float, atol: float, evolve_limit: int,
                   evolve_timeout: float) -> Dict[str, object]:
    """Benchmark one network size and return a record of the measurements."""
    net = copy.deepcopy(base_net)
    if zmax is not None:
        # The photon has to survive the cut, otherwise every photodisintegration
        # record is discarded along with it and the reduced network is not a
        # representative sample of the database.
        limit_network(net, select_species(net, zmax=zmax, amax=amax) + ["gamma"])
    # Count and evolve nuclides only; photons and leptons are participants in
    # the reaction records but are not part of the abundance vector.
    species = [s for s in net.species_names() if not is_massless(s)]
    n_species = len(species)
    n_reactions = len(net.reactions.reactions)
    record: Dict[str, object] = {
        "label": label,
        "zmax": zmax,
        "amax": amax,
        "species": n_species,
        "reactions": n_reactions,
    }

    if not net.zones:
        record["error"] = "no zone in database"
        return record
    zone = net.zones[0]
    y0 = np.array([zone.get_abundance(s) for s in species], dtype=float)
    thermo = constant_thermo(t9, rho)

    # Right-hand side.
    f = rhs(net, species, thermo)
    f(0.0, y0)  # warm up caches before timing
    record["rhs_seconds"] = _time_repeated(lambda: f(0.0, y0))

    # Jacobians.  The finite-difference matrix costs N+1 right-hand sides, so it
    # is skipped once that becomes the dominant cost of the whole benchmark.
    ja = analytic_jacobian(net, species, thermo, sparse=True)
    ja(0.0, y0)
    record["jacobian_analytic_seconds"] = _time_repeated(lambda: ja(0.0, y0), min_repeats=5, budget=2.0)
    if n_species <= evolve_limit:
        jf = jacobian(net, species, thermo, sparse=True)
        jf(0.0, y0)  # warm up before timing, as for the analytic matrix
        # Both Jacobians are timed the same way.  A single sample of either is
        # not reproducible: the analytic matrix is fast enough that scheduling
        # noise moves it by over 50 per cent between runs, which propagates
        # straight into the reported speedup.
        record["jacobian_fd_seconds"] = _time_repeated(lambda: jf(0.0, y0), min_repeats=3, budget=2.0)
        # Accuracy of the analytic matrix against finite differences.
        A = np.asarray(ja(0.0, y0).todense() if hasattr(ja(0.0, y0), "todense") else ja(0.0, y0))
        F = np.asarray(jf(0.0, y0).todense() if hasattr(jf(0.0, y0), "todense") else jf(0.0, y0))
        scale = max(float(np.abs(F).max()), 1e-300)
        record["jacobian_max_rel_diff"] = float(np.abs(A - F).max() / scale)
    else:
        record["jacobian_fd_seconds"] = None
        record["jacobian_max_rel_diff"] = None

    # Baryon-conservation residual of the right-hand side at the initial state.
    dy = f(0.0, y0)
    record["baryon_residual"] = float(sum(
        net.species[s].a * v for s, v in zip(species, dy) if s in net.species))

    # One-zone integration.  Large networks are reported without an evolution
    # rather than allowed to run unbounded.
    if n_species <= evolve_limit:
        times = time_grid(0.0, t_end, steps)
        t0 = time.perf_counter()
        result = evolve_zone(net, zone, times, thermo=thermo, method="bdf",
                             rtol=rtol, atol=atol)
        record["evolve_seconds"] = time.perf_counter() - t0
        record["evolve_success"] = bool(result.success)
        record["evolve_message"] = str(result.message)
        record["nfev"] = int(result.nfev)
        record["njev"] = int(result.njev)
        if result.success and len(result.y):
            final = result.final_abundances
            record["final_xsum"] = float(sum(
                net.species[k].a * v for k, v in final.items() if k in net.species))
    else:
        record["evolve_seconds"] = None
        record["evolve_success"] = None
        record["evolve_message"] = f"skipped: {n_species} species exceeds --evolve-limit"

    # Process high-water mark, not the cost of this case; see peak_memory_mb.
    record["process_peak_memory_mb"] = peak_memory_mb()
    return record


def latex_table(records: List[Dict[str, object]], env: Dict[str, str]) -> str:
    """Render the measurements as a LaTeX tabular body."""
    def fmt(value, spec: str, missing: str = "--"):
        return missing if value is None else format(value, spec)

    lines = [
        "% Generated by validation/benchmark_performance.py -- do not edit by hand.",
        "%% Environment: Python %s, NumPy %s, SciPy %s" % (env["python"], env["numpy"], env["scipy"]),
        "%% Platform: %s" % env["platform"],
        r"\small",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Case & Species & Reactions & RHS & $J$ analytic & $J$ numerical & "
        r"Evolution\\",
        r" & & & (ms) & (ms) & (ms) & (s)\\",
        r"\midrule",
    ]
    for r in records:
        if "error" in r:
            continue
        lines.append(
            "%s & %d & %d & %s & %s & %s & %s\\\\" % (
                r["label"], r["species"], r["reactions"],
                fmt(r.get("rhs_seconds") and r["rhs_seconds"] * 1e3, ".2f"),
                fmt(r.get("jacobian_analytic_seconds") and r["jacobian_analytic_seconds"] * 1e3, ".2f"),
                fmt(r.get("jacobian_fd_seconds") and r["jacobian_fd_seconds"] * 1e3, ".1f"),
                fmt(r.get("evolve_seconds"), ".1f"),
            )
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nuclides", required=True)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--zone")
    ap.add_argument("--t9", type=float, default=0.2)
    ap.add_argument("--rho", type=float, default=1.0e4)
    ap.add_argument("--t-end", type=float, default=1.0e2)
    ap.add_argument("--steps", type=int, default=60)
    ap.add_argument("--rtol", type=float, default=1.0e-6)
    ap.add_argument("--atol", type=float, default=1.0e-12)
    ap.add_argument("--evolve-limit", type=int, default=400,
                    help="skip integration and finite-difference Jacobian above this species count")
    ap.add_argument("--evolve-timeout", type=float, default=600.0)
    ap.add_argument("--json", default="benchmark.json")
    ap.add_argument("--latex", default="benchmark.tex")
    args = ap.parse_args()

    env = environment_record()
    print("environment:", json.dumps(env), flush=True)

    t0 = time.perf_counter()
    base = read_jina_xml(args.nuclides, args.reactions, args.zone)
    load_seconds = time.perf_counter() - t0
    validation = base.validate()
    print(f"loaded {len(base.species)} species, "
          f"{len(base.reactions.reactions)} reactions in {load_seconds:.2f} s", flush=True)

    records: List[Dict[str, object]] = []
    for label, zmax, amax in DEFAULT_LADDER:
        print(f"--- {label} ---", flush=True)
        rec = benchmark_case(base, label, zmax, amax, args.t9, args.rho,
                             args.t_end, args.steps, args.rtol, args.atol,
                             args.evolve_limit, args.evolve_timeout)
        records.append(rec)
        print(json.dumps(rec), flush=True)

    payload = {
        "environment": env,
        "database": {
            "nuclides_xml": args.nuclides,
            "reactions_xml": args.reactions,
            "zones_xml": args.zone,
            "species": len(base.species),
            "reactions": len(base.reactions.reactions),
            "zones": len(base.zones),
            "missing_species": len(validation["missing_species"]),
            "invalid_reactions": len(validation["invalid_reactions"]),
            "load_seconds": load_seconds,
        },
        "notes": {
            "process_peak_memory_mb": "ru_maxrss for the whole process, which "
            "never decreases and already includes the full database load; it "
            "is not the memory cost of the individual case",
        },
        "conditions": {"t9": args.t9, "rho": args.rho, "t_end": args.t_end,
                       "steps": args.steps, "rtol": args.rtol, "atol": args.atol},
        "cases": records,
    }
    Path(args.json).write_text(json.dumps(payload, indent=2))
    Path(args.latex).write_text(latex_table(records, env))
    print(f"wrote {args.json} and {args.latex}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
