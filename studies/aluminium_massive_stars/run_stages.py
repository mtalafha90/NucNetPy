"""Run the six core burning stages of the paper with NucNetPy.

Each stage holds the temperature and density of the paper's Table 1 for the
stage's duration and starts from the composition the previous stage left.
Hydrogen burning starts from the Lodders (2003) composition built by
build_inputs.py.

By default each stage is integrated by implicit Euler steps with
Newton-Raphson iterations, the scheme of NucNet Tools, built from NucNetPy's
right-hand side and analytic Jacobian (see stepper.py).  ``--solver bdf``
uses ``nucnetpy.evolve_zone`` with SciPy's BDF method instead, which fails on
this network at hydrogen exhaustion (see README.md).

For every stage the script writes ``results/<stage>.npz`` (times in seconds,
species names and abundances on a logarithmic time grid) and adds the stage's
diagnostics to ``results/run_summary.json``.

    python run_stages.py                        # all six stages
    python run_stages.py --stages neon oxygen   # resume from saved stages
    python run_stages.py --screening skynet --tag screened
    python run_stages.py --max-change 0.01 --tag half_steps
    python run_stages.py --stages silicon --silicon-unit yr --tag silicon_year \\
        --from-results results

``--screening skynet`` applies NucNetPy's SkyNet plasma screening; the paper
does not mention screening, so the default is none.  ``--max-change`` sets
the largest fractional change of an abundance in one step, and so the
accuracy.  ``--silicon-unit yr`` runs silicon burning for 0.733 years, the
duration in the paper's Table 1, instead of the 0.733 days its figures show;
``--from-results`` takes the composition the stage starts from out of
another results directory.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from common import RESULTS, STAGES, RateFactors, load_composition, load_network

#: Abundances kept in the saved trajectories.  A species is saved if its mass
#: fraction exceeds this anywhere in the stage.
SAVE_FLOOR = 1.0e-30


def assign_reverse_partition_functions(net):
    """Put the partition-function factor on the rates ReacLib flags as reverse.

    NucNetPy 1.0.0 recognises a ReacLib reverse rate from its coefficients
    and assumes that the reverse member of a pair is the endothermic one.
    For 168 pairs in this network ReacLib's forward rate is the endothermic
    one (a capture into a particle-unbound nucleus near the drip lines), so
    NucNetPy puts the factor on the capture instead of the photodisintegration.
    Three more pairs are missed altogether because the library rounds their
    large ``a6`` coefficients; one of them is 29Si(g,p)28Al, whose factor is
    1.64 at T9 = 3.6.  The network file keeps ReacLib's own flag (a trailing
    ``v`` on each fit label), so the factor is placed on the flagged member
    here.  Returns the number of reactions whose assignment changed.
    """
    from nucnetpy.detailed_balance import _partition_table
    from nucnetpy.reactions import PartitionRatio
    from nucnetpy.species import is_massless

    changed = 0
    for r in net.reactions.reactions:
        flagged = bool(r.rate_fits) and all(f.label.endswith("v") for f in r.rate_fits)
        if flagged and r.partition_ratio is None:
            num = tuple(_partition_table(net.species, p.species) + (p.count,)
                        for p in r.products if not is_massless(p.species))
            den = tuple(_partition_table(net.species, p.species) + (p.count,)
                        for p in r.reactants if not is_massless(p.species))
            r.partition_ratio = PartitionRatio(num, den)
            changed += 1
        elif not flagged and r.partition_ratio is not None:
            r.partition_ratio = None
            changed += 1
    return changed


def output_times(stage, points: int) -> np.ndarray:
    """t = 0 followed by a logarithmic grid from the paper's first plotted time."""
    from common import UNITS
    unit = UNITS[stage.unit]
    grid = np.geomspace(stage.first_output * unit, stage.seconds, points)
    return np.concatenate(([0.0], grid))


class Progress:
    """Report how far a long integration has got, every ``every`` calls."""

    def __init__(self, every: int):
        self.every, self.calls = int(every), 0
        self.start = time.perf_counter()

    def report(self, t, what):
        self.calls += 1
        if self.calls % self.every == 0:
            print(f"    t = {t:.3e} s after {self.calls} {what}, "
                  f"{time.perf_counter() - self.start:.0f} s", flush=True)

    def thermo(self, t9, rho):
        """A constant-state function for evolve_zone that also reports progress."""
        def state(t, abundances):
            self.report(t, "evaluations")
            return float(t9), float(rho)
        return state


def run_stage(net, stage, abundances, args, screening):
    times = output_times(stage, args.points)
    t0 = time.perf_counter()
    if args.solver == "bdf":
        from nucnetpy.core import Zone
        from nucnetpy.solver import evolve_zone
        result = evolve_zone(net, Zone(abundances=dict(abundances)), times,
                             thermo=Progress(2000).thermo(stage.t9, stage.rho),
                             method="bdf", rtol=args.rtol, atol=args.atol,
                             screening=screening)
        if not result.success or len(result.time) != len(times):
            raise RuntimeError(f"{stage.key}: {result.message}")
        stats = {"nfev": int(result.nfev), "njev": int(result.njev)}
    else:
        from stepper import evolve
        progress = Progress(500)
        result = evolve(net, abundances, times, stage.t9, stage.rho, screening=screening,
                        max_change=args.max_change, x_min=args.x_min,
                        progress=lambda t, n: progress.report(t, "steps"))
        if not result.success:
            raise RuntimeError(f"{stage.key}: {result.message}")
        stats = {"steps": result.steps, "rejected_steps": result.rejected,
                 "nfev": result.nfev, "njev": result.njev,
                 "max_negative_mass_fraction_clipped": result.max_negative_clipped}
    seconds = time.perf_counter() - t0
    a = np.array([net.species[s].a for s in result.species], dtype=float)
    z = np.array([net.species[s].z for s in result.species], dtype=float)
    x = result.y * a[None, :]
    keep = np.max(x, axis=0) > SAVE_FLOOR
    info = {
        "t9": stage.t9, "rho": stage.rho, "duration": stage.duration,
        "unit": stage.unit, "seconds": stage.seconds,
        "success": bool(result.success), "message": result.message, **stats,
        "wall_seconds": round(seconds, 1),
        "max_abs_xsum_minus_1": float(np.max(np.abs(x.sum(axis=1) - 1.0))),
        "ye_initial": float(result.y[0] @ z), "ye_final": float(result.y[-1] @ z),
        "species_saved": int(keep.sum()),
    }
    final = {s: float(v) for s, v in zip(result.species, result.y[-1])}
    return result, keep, info, final


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stages", nargs="+", choices=[s.key for s in STAGES],
                    default=[s.key for s in STAGES])
    ap.add_argument("--screening", choices=["none", "skynet"], default="none")
    ap.add_argument("--solver", choices=["implicit-euler", "bdf"], default="implicit-euler")
    ap.add_argument("--max-change", type=float, default=0.02,
                    help="implicit Euler: largest fractional change of an abundance per step")
    ap.add_argument("--x-min", type=float, default=1.0e-15,
                    help="implicit Euler: mass fraction below which an abundance does "
                         "not limit the step")
    ap.add_argument("--rtol", type=float, default=1.0e-6, help="BDF only")
    ap.add_argument("--atol", type=float, default=1.0e-20, help="BDF only")
    ap.add_argument("--points", type=int, default=300,
                    help="output times per stage, besides t = 0")
    ap.add_argument("--tag", default="", help="suffix for the output directory")
    ap.add_argument("--silicon-unit", choices=["d", "yr"], default="d",
                    help="unit of the 0.733 duration of silicon burning")
    ap.add_argument("--from-results", type=Path, default=None,
                    help="results directory holding the stage before the first one run")
    args = ap.parse_args()

    outdir = RESULTS if not args.tag else RESULTS.parent / f"results_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)
    stages = {s.key: s for s in STAGES}
    if args.silicon_unit != "d":
        from dataclasses import replace
        stages["silicon"] = replace(stages["silicon"], unit=args.silicon_unit)
    summary_path = outdir / "run_summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}

    net = load_network()
    reassigned = assign_reverse_partition_functions(net)
    print(f"{len(net.species)} species, {len(net.reactions.reactions)} reactions; "
          f"reverse-rate factor reassigned on {reassigned} reactions", flush=True)
    screening = None
    if args.screening == "skynet":
        from nucnetpy.screening import SkyNetScreening
        screening = SkyNetScreening(net.species)
    hook = RateFactors(screening)

    order = [s.key for s in STAGES]
    first = order.index(args.stages[0])
    if first == 0:
        abundances = {k: v / net.species[k].a for k, v in load_composition().items()}
    else:
        source = args.from_results or outdir
        previous = np.load(source / f"{order[first - 1]}.npz")
        abundances = dict(zip(previous["final_species"], previous["final_y"]))

    import nucnetpy
    import scipy
    summary["settings"] = {
        "solver": args.solver, "screening": args.screening,
        **({"rtol": args.rtol, "atol": args.atol} if args.solver == "bdf" else
           {"max_change": args.max_change, "x_min": args.x_min}),
        "points": args.points, "reverse_rate_factor_reassigned": reassigned,
        "electron_captures_times_rho_ye": True, "silicon_unit": args.silicon_unit,
        "started_from": str(args.from_results or outdir.name),
        "nucnetpy": nucnetpy.__version__, "numpy": np.__version__,
        "scipy": scipy.__version__, "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for key in order[first:]:
        if key not in args.stages:
            break
        stage = stages[key]
        print(f"{stage.title}: T9 = {stage.t9}, rho = {stage.rho:g} g/cm^3, "
              f"{stage.duration} {stage.unit} ...", flush=True)
        result, keep, info, final = run_stage(net, stage, abundances, args, hook)
        names = np.array(result.species)
        np.savez_compressed(
            outdir / f"{key}.npz",
            time=result.time, species=names[keep],
            y=result.y[:, keep].astype(np.float32),
            final_species=names, final_y=result.y[-1])
        summary[key] = info
        summary_path.write_text(json.dumps(summary, indent=1) + "\n")
        print(f"  done in {info['wall_seconds']} s, "
              f"|sum X - 1| <= {info['max_abs_xsum_minus_1']:.1e}", flush=True)
        abundances = final
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
