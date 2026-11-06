"""Run the six core burning stages of the paper with NucNetPy.

Each stage holds the temperature and density of the paper's Table 1 for the
stage's duration and starts from the composition the previous stage left.
Hydrogen burning starts from the Lodders (2003) composition built by
build_inputs.py.  The network is integrated with ``nucnetpy.evolve_zone``
(SciPy's BDF method with NucNetPy's analytic Jacobian).

For every stage the script writes ``results/<stage>.npz`` (times in seconds,
species names and abundances on a logarithmic time grid) and adds the stage's
diagnostics to ``results/run_summary.json``.

    python run_stages.py                        # all six stages
    python run_stages.py --stages neon oxygen   # resume from saved stages
    python run_stages.py --screening skynet --tag screened

``--screening skynet`` applies NucNetPy's SkyNet plasma screening; the paper
does not mention screening, so the default is none.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from common import (RESULTS, STAGES, STAGE_BY_KEY, RateFactors, load_composition,
                    load_network)

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
    """Constant temperature and density, reporting progress as it is asked for them.

    NucNetPy asks the thermodynamic function for the state once per
    right-hand-side or Jacobian evaluation, which makes it a convenient place
    to report how far a long integration has got.
    """

    def __init__(self, t9: float, rho: float, every: int = 2000):
        self.t9, self.rho, self.every = float(t9), float(rho), int(every)
        self.calls = 0
        self.start = time.perf_counter()

    def __call__(self, t, abundances):
        self.calls += 1
        if self.calls % self.every == 0:
            print(f"    t = {t:.3e} s after {self.calls} evaluations, "
                  f"{time.perf_counter() - self.start:.0f} s", flush=True)
        return self.t9, self.rho


def run_stage(net, stage, abundances, args, screening):
    from nucnetpy.core import Zone
    from nucnetpy.solver import evolve_zone

    zone = Zone(abundances=dict(abundances))
    times = output_times(stage, args.points)
    t0 = time.perf_counter()
    result = evolve_zone(net, zone, times, thermo=Progress(stage.t9, stage.rho),
                         method="bdf", rtol=args.rtol, atol=args.atol,
                         screening=screening)
    seconds = time.perf_counter() - t0
    if not result.success or len(result.time) != len(times):
        raise RuntimeError(f"{stage.key}: {result.message}")
    a = np.array([net.species[s].a for s in result.species], dtype=float)
    z = np.array([net.species[s].z for s in result.species], dtype=float)
    x = result.y * a[None, :]
    keep = np.max(x, axis=0) > SAVE_FLOOR
    xsum = x.sum(axis=1)
    info = {
        "t9": stage.t9, "rho": stage.rho, "duration": stage.duration,
        "unit": stage.unit, "seconds": stage.seconds,
        "success": bool(result.success), "message": result.message,
        "nfev": int(result.nfev), "njev": int(result.njev),
        "wall_seconds": round(seconds, 1),
        "max_abs_xsum_minus_1": float(np.max(np.abs(xsum - 1.0))),
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
    ap.add_argument("--rtol", type=float, default=1.0e-6)
    ap.add_argument("--atol", type=float, default=1.0e-20)
    ap.add_argument("--points", type=int, default=300,
                    help="output times per stage, besides t = 0")
    ap.add_argument("--tag", default="", help="suffix for the output directory")
    args = ap.parse_args()

    outdir = RESULTS if not args.tag else RESULTS.parent / f"results_{args.tag}"
    outdir.mkdir(parents=True, exist_ok=True)
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
        previous = np.load(outdir / f"{order[first - 1]}.npz")
        abundances = dict(zip(previous["species"], previous["final_y"]))

    import nucnetpy
    import scipy
    summary["settings"] = {
        "screening": args.screening, "rtol": args.rtol, "atol": args.atol,
        "points": args.points, "reverse_rate_factor_reassigned": reassigned,
        "electron_captures_times_rho_ye": True,
        "nucnetpy": nucnetpy.__version__, "numpy": np.__version__,
        "scipy": scipy.__version__, "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for key in order[first:]:
        if key not in args.stages:
            break
        stage = STAGE_BY_KEY[key]
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
        print(f"  done in {info['wall_seconds']} s, {info['nfev']} evaluations, "
              f"|sum X - 1| <= {info['max_abs_xsum_minus_1']:.1e}", flush=True)
        abundances = final
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
