"""Reproduce the demonstration calculation of the NucNetPy manuscript.

The script runs a silicon-burning problem drawn entirely from a JINA /
libnucnet database and checks the result against a separately solved nuclear
statistical equilibrium (NSE) composition.  As run here the two sides share no
numerical input: the network integrates the library's own ReacLib forward and
reverse fits in time, while the NSE solve builds the Saha equations from the
mass table.  Their disagreement is therefore informative about the data, and
is what the 4.5 per cent result measures.

That independence does not survive rebuilding the reverse rates with
consistent_reverse_network, which derives them from the same equilibrium
prefactor solve_nse uses.  The parts-per-million agreement that follows
measures how consistently the integrator, the stoichiometry and the
equilibrium solver treat one shared formulation; it does not test the
formulation itself.

It writes a JSON record, LaTeX tables, and figures, so that every demonstration
number and plot in the manuscript can be regenerated with one command.

Usage
-----
python validation/demonstration.py \
    --nuclides /path/to/nuclides.xml \
    --reactions /path/to/reaction_data.xml \
    --outdir results/

Add ``--zone /path/to/zone.xml`` to additionally run the hydrogen-burning case
on the composition stored in that zone file.
"""
from __future__ import annotations

import argparse
import copy
import json
import platform
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from nucnetpy.core import Zone
from nucnetpy.io.jina import read_jina_xml
from nucnetpy.detailed_balance import consistent_reverse_network
from nucnetpy.network_limiter import limit_network, select_species
from nucnetpy.nse import solve_nse
from nucnetpy.solver import constant_thermo, evolve_zone, rhs, time_grid

#: Classical alpha chain plus free nucleons and the photon.  The photon must be
#: retained so that photodisintegration records survive the network cut.
ALPHA_CHAIN = ["gamma", "n", "h1", "he4", "c12", "o16", "ne20", "mg24", "si28",
               "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56"]


def environment_record() -> Dict[str, str]:
    import scipy
    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }


def mass_fractions(net, abundances: Dict[str, float]) -> Dict[str, float]:
    return {k: net.species[k].a * v for k, v in abundances.items() if k in net.species}


def silicon_burning(base_net, t9: float, rho: float, t_end: float, steps: int,
                    rtol: float, atol: float) -> Dict[str, object]:
    """Burn pure silicon-28 towards the iron group and compare with NSE."""
    net = copy.deepcopy(base_net)
    limit_network(net, ALPHA_CHAIN)
    zone = Zone(abundances={"si28": 1.0 / 28.0})
    net.zones = [zone]
    ye = zone.ye(net.species)

    times = time_grid(0.0, t_end, steps)
    t0 = time.perf_counter()
    result = evolve_zone(net, zone, times, thermo=constant_thermo(t9, rho),
                         method="bdf", rtol=rtol, atol=atol)
    evolve_seconds = time.perf_counter() - t0

    t0 = time.perf_counter()
    nse = solve_nse(net, t9=t9, rho=rho, ye=ye)
    nse_seconds = time.perf_counter() - t0

    record: Dict[str, object] = {
        "case": "silicon burning",
        "t9": t9, "rho": rho, "t_end": t_end, "steps": steps,
        "rtol": rtol, "atol": atol,
        "species": len(net.species), "reactions": len(net.reactions.reactions),
        "initial_ye": ye,
        "evolve_seconds": evolve_seconds, "nse_seconds": nse_seconds,
        "evolve_success": bool(result.success), "evolve_message": str(result.message),
        "nfev": int(result.nfev), "njev": int(result.njev),
        "nse_success": bool(nse.success),
        "nse_xsum": float(nse.xsum), "nse_ye": float(nse.computed_ye),
        "mu_p": float(nse.mu_p), "mu_n": float(nse.mu_n),
    }
    if not result.success:
        return record

    final = result.final_abundances
    x_net = mass_fractions(net, final)
    x_nse = mass_fractions(net, nse.abundances)
    record["final_xsum"] = float(sum(x_net.values()))
    record["final_ye"] = float(sum(net.species[k].z * v for k, v in final.items()
                                   if k in net.species))

    # Compare species that the restricted network is able to populate.  Free
    # nucleons are reported separately: an alpha chain has no reaction that
    # liberates a nucleon, so it relaxes to a constrained equilibrium in which
    # Y(n) and Y(p) stay at their initial values.
    comparison = []
    for name in sorted(x_nse, key=lambda k: -x_nse[k]):
        if name in {"n", "h1"}:
            continue
        if x_nse[name] < 1.0e-6:
            continue
        net_x = x_net.get(name, 0.0)
        comparison.append({
            "species": name,
            "x_network": float(net_x),
            "x_nse": float(x_nse[name]),
            "rel_diff": float(abs(net_x - x_nse[name]) / x_nse[name]),
        })
    record["nse_comparison"] = comparison
    if comparison:
        record["max_rel_diff_vs_nse"] = max(c["rel_diff"] for c in comparison)
        record["median_rel_diff_vs_nse"] = float(
            np.median([c["rel_diff"] for c in comparison]))

    # Baryon conservation of the right-hand side along the trajectory.
    f = rhs(net, result.species, constant_thermo(t9, rho))
    residuals = []
    for i in range(0, len(result.time), max(1, len(result.time) // 20)):
        dy = f(float(result.time[i]), result.y[i])
        residuals.append(abs(sum(net.species[s].a * v
                                 for s, v in zip(result.species, dy)
                                 if s in net.species)))
    record["max_baryon_residual"] = float(max(residuals))

    record["_trajectory"] = {
        "time": result.time.tolist(),
        "species": result.species,
        "mass_fractions": result.mass_fraction_history(net.species).tolist(),
        "nse_mass_fractions": {k: float(v) for k, v in x_nse.items()},
    }
    return record


def detailed_balance_study(base_net, t9: float, rho: float, t_end: float,
                           steps: int, rtol: float, atol: float) -> Dict[str, object]:
    """Repeat the silicon burn with reverse rates rebuilt from detailed balance.

    Two variants are recorded.  The first attaches each reverse rate as a
    function of temperature, which is exact but not writable to XML.  The
    second tabulates it on a logarithmic grid, which is what makes a network
    serialisable and costs accuracy.

    The agreement this produces is not an independent check: the reverse rates
    are derived from the same equilibrium prefactor ``solve_nse`` uses, so the
    two sides share their nuclear data.  It measures how consistently the
    integrator, the stoichiometry and the equilibrium solver treat one shared
    formulation.
    """
    out: Dict[str, object] = {"t9": t9, "rho": rho, "t_end": t_end,
                              "steps": steps, "rtol": rtol, "atol": atol}
    for label, tabulate in (("function", False), ("tabulated", True)):
        net = copy.deepcopy(base_net)
        limit_network(net, ALPHA_CHAIN)
        net = consistent_reverse_network(net, tabulate=tabulate)
        zone = Zone(abundances={"si28": 1.0 / 28.0})
        net.zones = [zone]
        ye = zone.ye(net.species)

        t0 = time.perf_counter()
        result = evolve_zone(net, zone, time_grid(0.0, t_end, steps),
                             thermo=constant_thermo(t9, rho), method="bdf",
                             rtol=rtol, atol=atol)
        seconds = time.perf_counter() - t0
        rec: Dict[str, object] = {
            "tabulated": tabulate,
            "seconds": seconds,
            "evolve_success": bool(result.success),
            "evolve_message": str(result.message),
            "reactions": len(net.reactions.reactions),
        }
        if result.success and len(result.y):
            nse = solve_nse(net, t9=t9, rho=rho, ye=ye)
            x_net = mass_fractions(net, result.final_abundances)
            x_nse = mass_fractions(net, nse.abundances)
            diffs = [abs(x_net.get(k, 0.0) - v) / v
                     for k, v in x_nse.items()
                     if k not in {"n", "h1"} and v >= 1.0e-6]
            if diffs:
                rec["median_rel_diff_vs_nse"] = float(np.median(diffs))
                rec["max_rel_diff_vs_nse"] = float(max(diffs))
            rec["final_xsum"] = float(sum(x_net.values()))
        out[label] = rec
    return out


def solver_comparison(base_net, t9: float, rho: float, t_end: float, steps: int,
                      methods: List[str], rtol: float, atol: float) -> List[Dict[str, object]]:
    """Integrate the same problem with several solvers and compare endpoints."""
    net = copy.deepcopy(base_net)
    limit_network(net, ALPHA_CHAIN)
    zone = Zone(abundances={"si28": 1.0 / 28.0})
    net.zones = [zone]
    times = time_grid(0.0, t_end, steps)

    rows: List[Dict[str, object]] = []
    reference: Optional[Dict[str, float]] = None
    for method in methods:
        t0 = time.perf_counter()
        result = evolve_zone(net, zone, times, thermo=constant_thermo(t9, rho),
                             method=method, rtol=rtol, atol=atol)
        seconds = time.perf_counter() - t0
        row: Dict[str, object] = {
            "method": method, "seconds": seconds,
            "success": bool(result.success), "message": str(result.message),
            "nfev": int(result.nfev), "njev": int(result.njev),
        }
        if result.success and len(result.y):
            final = result.final_abundances
            row["xsum"] = float(sum(mass_fractions(net, final).values()))
            if reference is None:
                reference = final
                row["max_rel_diff_vs_reference"] = 0.0
            else:
                diffs = [abs(final[k] - reference[k]) / max(abs(reference[k]), 1e-20)
                         for k in reference if reference[k] > 1e-12]
                row["max_rel_diff_vs_reference"] = float(max(diffs)) if diffs else 0.0
        rows.append(row)
    return rows


def tolerance_study(base_net, t9: float, rho: float, t_end: float, steps: int,
                    tolerances: List[float]) -> List[Dict[str, object]]:
    """Tighten the solver tolerance and track the change in the endpoint."""
    net = copy.deepcopy(base_net)
    limit_network(net, ALPHA_CHAIN)
    zone = Zone(abundances={"si28": 1.0 / 28.0})
    net.zones = [zone]
    times = time_grid(0.0, t_end, steps)

    rows: List[Dict[str, object]] = []
    previous: Optional[Dict[str, float]] = None
    for rtol in tolerances:
        t0 = time.perf_counter()
        result = evolve_zone(net, zone, times, thermo=constant_thermo(t9, rho),
                             method="bdf", rtol=rtol, atol=rtol * 1e-6)
        seconds = time.perf_counter() - t0
        row: Dict[str, object] = {
            "rtol": rtol, "atol": rtol * 1e-6, "seconds": seconds,
            "success": bool(result.success), "nfev": int(result.nfev),
            "njev": int(result.njev),
        }
        if result.success and len(result.y):
            final = result.final_abundances
            row["xsum"] = float(sum(mass_fractions(net, final).values()))
            if previous is not None:
                diffs = [abs(final[k] - previous[k]) / max(abs(previous[k]), 1e-20)
                         for k in previous if previous[k] > 1e-12]
                row["max_rel_change"] = float(max(diffs)) if diffs else 0.0
            previous = final
        rows.append(row)
    return rows


def hydrogen_burning(base_net, t9: float, rho: float, t_end: float, steps: int,
                     zmax: int, amax: int, rtol: float, atol: float) -> Dict[str, object]:
    """Burn the composition supplied in the zone file on a reduced network."""
    net = copy.deepcopy(base_net)
    limit_network(net, select_species(net, zmax=zmax, amax=amax) + ["gamma"])
    if not net.zones:
        return {"case": "hydrogen burning", "error": "no zone supplied"}
    zone = net.zones[0]
    initial = mass_fractions(net, zone.abundances)

    times = time_grid(0.0, t_end, steps)
    t0 = time.perf_counter()
    result = evolve_zone(net, zone, times, thermo=constant_thermo(t9, rho),
                         method="bdf", rtol=rtol, atol=atol)
    seconds = time.perf_counter() - t0

    record: Dict[str, object] = {
        "case": "hydrogen burning",
        "t9": t9, "rho": rho, "t_end": t_end, "steps": steps,
        "zmax": zmax, "amax": amax, "rtol": rtol, "atol": atol,
        "species": len(net.species), "reactions": len(net.reactions.reactions),
        "initial_xsum": float(sum(initial.values())),
        "initial_ye": float(zone.ye(net.species)),
        "evolve_seconds": seconds,
        "evolve_success": bool(result.success), "evolve_message": str(result.message),
        "nfev": int(result.nfev), "njev": int(result.njev),
    }
    # A run that the conservation check has flagged still carries useful
    # diagnostics, so only a run that produced no output at all is abandoned.
    if not len(result.y):
        return record
    final = result.final_abundances
    x_final = mass_fractions(net, final)
    record["final_xsum"] = float(sum(x_final.values()))
    record["final_ye"] = float(sum(net.species[k].z * v for k, v in final.items()
                                   if k in net.species))
    record["largest_final"] = sorted(
        ({"species": k, "x": float(v)} for k, v in x_final.items() if v > 1e-8),
        key=lambda d: -d["x"])[:15]
    # How much baryon number the positivity projection had to invent.  Anything
    # above the requested relative tolerance means the trajectory is not
    # trustworthy, however cleanly the integrator reported finishing.
    if "positivity projection created" in result.message:
        record["projection_created"] = float(
            result.message.split("positivity projection created ")[1].split()[0])
    else:
        record["projection_created"] = 0.0
    # Repeat without positivity projection.  This is what shows the mechanism:
    # baryon number is then conserved to round-off, but a component is carried
    # far negative, and it is projecting that away which invents the mass.
    unclipped = evolve_zone(net, zone, times, thermo=constant_thermo(t9, rho),
                            method="bdf", rtol=rtol, atol=atol,
                            project_positive=False)
    if len(unclipped.y):
        y_end = unclipped.y[-1]
        i = int(np.argmin(y_end))
        record["unclipped"] = {
            "evolve_success": bool(unclipped.success),
            "min_abundance": float(y_end[i]),
            "min_abundance_species": unclipped.species[i],
            "xsum": float(sum(net.species[s].a * v
                              for s, v in zip(unclipped.species, y_end)
                              if s in net.species)),
        }

    record["_trajectory"] = {
        "time": result.time.tolist(),
        "species": result.species,
        "mass_fractions": result.mass_fraction_history(net.species).tolist(),
    }
    return record


def _label(name: str) -> str:
    """Render a species name as a nuclear symbol, e.g. ni56 -> $^{56}$Ni."""
    import re
    m = re.match(r"^([a-z]+)(\d+)$", name)
    if not m:
        return name
    return "$^{%s}$%s" % (m.group(2), m.group(1).capitalize())


def write_figures(silicon: Dict[str, object], hydrogen: Optional[Dict[str, object]],
                  outdir: Path) -> List[str]:
    """Write the manuscript figures; returns the filenames produced."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"matplotlib unavailable ({exc}); skipping figures")
        return []

    written: List[str] = []
    traj = silicon.get("_trajectory")
    if traj:
        t = np.asarray(traj["time"], dtype=float)
        names = traj["species"]
        X = np.asarray(traj["mass_fractions"], dtype=float)
        nse_x = traj["nse_mass_fractions"]
        # Plot the species that end up dominant.
        order = np.argsort(-X[-1])
        keep = [i for i in order if X[-1, i] > 1e-4][:8]

        # Almost all of the burning happens in the first fraction of a second,
        # so a logarithmic time axis is needed for the approach to equilibrium
        # to be visible at all.  The t=0 sample cannot be shown on it.
        pos = t > 0.0

        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(6.6, 6.2), sharex=True,
            gridspec_kw={"height_ratios": [2.0, 1.0], "hspace": 0.08})
        for i in keep:
            line, = ax.plot(t[pos], X[pos, i], lw=1.7, label=_label(names[i]))
            if names[i] in nse_x:
                ax.axhline(nse_x[names[i]], color=line.get_color(),
                           ls=":", lw=1.1, alpha=0.85)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(t[pos][0], t[-1])
        ax.set_ylim(1e-4, 3.0)
        ax.set_ylabel("mass fraction $X_i$")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(ncol=4, fontsize=8.5, frameon=False,
                  loc="lower center", bbox_to_anchor=(0.5, 1.02))

        # Lower panel: relative distance from the NSE composition.
        for i in keep:
            if names[i] in nse_x and nse_x[names[i]] > 0:
                axr.plot(t[pos],
                         np.abs(X[pos, i] - nse_x[names[i]]) / nse_x[names[i]],
                         lw=1.4)
        axr.set_xscale("log")
        axr.set_yscale("log")
        axr.set_ylim(1e-4, 5e2)
        axr.grid(alpha=0.25, lw=0.5)
        axr.set_xlabel("time (s)")
        axr.set_ylabel(r"$|X_i-X_i^{\rm NSE}|\,/\,X_i^{\rm NSE}$")
        fig.tight_layout()
        path = outdir / "fig_silicon_burning.pdf"
        fig.savefig(path)
        plt.close(fig)
        written.append(path.name)

    if hydrogen and hydrogen.get("_trajectory"):
        traj = hydrogen["_trajectory"]
        t = np.asarray(traj["time"], dtype=float)
        names = traj["species"]
        X = np.asarray(traj["mass_fractions"], dtype=float)
        order = np.argsort(-X[-1])
        keep = [i for i in order if X[-1, i] > 1e-6][:10]
        pos = t > 0.0
        fig, ax = plt.subplots(figsize=(6.6, 4.3))
        for i in keep:
            ax.plot(t[pos], X[pos, i], lw=1.5, label=_label(names[i]))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(t[pos][0], t[-1])
        ax.set_xlabel("time (s)")
        ax.set_ylabel("mass fraction $X_i$")
        ax.grid(alpha=0.25, lw=0.5)
        ax.legend(ncol=5, fontsize=8.5, frameon=False,
                  loc="lower center", bbox_to_anchor=(0.5, 1.02))
        fig.tight_layout()
        path = outdir / "fig_hydrogen_burning.pdf"
        fig.savefig(path)
        plt.close(fig)
        written.append(path.name)
    return written


def _sci(value: float) -> str:
    """Render a number as LaTeX scientific notation, e.g. $5.1\\times10^{-12}$."""
    import math
    if value is None:
        return "--"
    if value == 0.0 or not math.isfinite(value):
        return "0" if value == 0.0 else "--"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10.0 ** exponent)
    return r"$%.1f\times10^{%d}$" % (mantissa, exponent)


def _pow10(value: float) -> str:
    """Render an exact power of ten as $10^{n}$."""
    import math
    exponent = int(round(math.log10(value)))
    return r"$10^{%d}$" % exponent


def latex_nse_table(record: Dict[str, object]) -> str:
    rows = record.get("nse_comparison") or []
    lines = [
        "% Generated by validation/demonstration.py -- do not edit by hand.",
        r"\begin{tabular}{lrrr}",
        r"\toprule",
        r"Species & $X_i$ (network) & $X_i$ (NSE) & Relative difference\\",
        r"\midrule",
    ]
    for row in rows:
        lines.append("%s & %.4e & %.4e & %.2e\\\\" % (
            row["species"], row["x_network"], row["x_nse"], row["rel_diff"]))
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def latex_solver_table(rows: List[Dict[str, object]]) -> str:
    lines = [
        "% Generated by validation/demonstration.py -- do not edit by hand.",
        r"\small",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        r"Method & Time & $n_{\rm fev}$ & $n_{\rm jev}$ & $1-\sum_i A_iY_i$ & "
        r"Max.\ rel.\\",
        r" & (s) & & & & difference\\",
        r"\midrule",
    ]
    for r in rows:
        if not r.get("success"):
            lines.append("%s & \\multicolumn{5}{l}{failed: %s}\\\\" % (
                r["method"].upper(), str(r.get("message", ""))[:40]))
            continue
        diff = r.get("max_rel_diff_vs_reference")
        lines.append("%s & %.2f & %d & %d & %s & %s\\\\" % (
            r["method"].upper(), r["seconds"], r["nfev"], r["njev"],
            _sci(1.0 - r["xsum"]) if "xsum" in r else "--",
            "reference" if diff == 0.0 else (_sci(diff) if diff is not None else "--")))
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def latex_tolerance_table(rows: List[Dict[str, object]]) -> str:
    lines = [
        "% Generated by validation/demonstration.py -- do not edit by hand.",
        r"\small",
        r"\begin{tabular}{llrrrl}",
        r"\toprule",
        r"$r_{\rm tol}$ & $a_{\rm tol}$ & Time & $n_{\rm fev}$ & "
        r"$1-\sum_i A_iY_i$ & Change vs.\\",
        r" & & (s) & & & previous\\",
        r"\midrule",
    ]
    for r in rows:
        if not r.get("success"):
            lines.append(r"%s & %s & \multicolumn{4}{l}{failed}\\" % (
                _pow10(r["rtol"]), _pow10(r["atol"])))
            continue
        change = r.get("max_rel_change")
        lines.append("%s & %s & %.2f & %d & %s & %s\\\\" % (
            _pow10(r["rtol"]), _pow10(r["atol"]), r["seconds"], r["nfev"],
            _sci(1.0 - r["xsum"]) if "xsum" in r else "--",
            "--" if change is None else _sci(change)))
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nuclides", required=True)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--zone")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--t9", type=float, default=5.0)
    ap.add_argument("--rho", type=float, default=1.0e8)
    ap.add_argument("--t-end", type=float, default=10.0)
    ap.add_argument("--steps", type=int, default=120)
    ap.add_argument("--rtol", type=float, default=1.0e-8)
    ap.add_argument("--atol", type=float, default=1.0e-14)
    ap.add_argument("--h-t9", type=float, default=0.2)
    ap.add_argument("--h-rho", type=float, default=1.0e4)
    ap.add_argument("--h-t-end", type=float, default=100.0)
    ap.add_argument("--h-t-end-long", type=float, default=1.0e6,
                    help="second, longer hydrogen-burning integration")
    ap.add_argument("--h-zmax", type=int, default=8)
    ap.add_argument("--h-amax", type=int, default=20)
    ap.add_argument("--skip-hydrogen", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    env = environment_record()
    print("environment:", json.dumps(env), flush=True)

    t0 = time.perf_counter()
    base = read_jina_xml(args.nuclides, args.reactions, args.zone)
    load_seconds = time.perf_counter() - t0
    validation = base.validate()
    print(f"loaded {len(base.species)} species, {len(base.reactions.reactions)} "
          f"reactions, {len(base.zones)} zones in {load_seconds:.2f} s", flush=True)

    print("--- silicon burning ---", flush=True)
    silicon = silicon_burning(base, args.t9, args.rho, args.t_end, args.steps,
                              args.rtol, args.atol)
    print(json.dumps({k: v for k, v in silicon.items() if not k.startswith("_")},
                     indent=2), flush=True)

    print("--- detailed balance ---", flush=True)
    detailed = detailed_balance_study(base, args.t9, args.rho, args.t_end,
                                      args.steps, args.rtol, args.atol)
    print(json.dumps(detailed, indent=2), flush=True)

    print("--- solver comparison ---", flush=True)
    solvers = solver_comparison(base, args.t9, args.rho, args.t_end, args.steps,
                                ["bdf", "radau", "lsoda"], args.rtol, args.atol)
    print(json.dumps(solvers, indent=2), flush=True)

    print("--- tolerance study ---", flush=True)
    tolerances = tolerance_study(base, args.t9, args.rho, args.t_end, args.steps,
                                 [1e-4, 1e-6, 1e-8, 1e-10])
    print(json.dumps(tolerances, indent=2), flush=True)

    hydrogen = None
    hydrogen_long = None
    if args.zone and not args.skip_hydrogen:
        print("--- hydrogen burning ---", flush=True)
        hydrogen = hydrogen_burning(base, args.h_t9, args.h_rho, args.h_t_end,
                                    args.steps, args.h_zmax, args.h_amax,
                                    1e-6, 1e-12)
        print(json.dumps({k: v for k, v in hydrogen.items() if not k.startswith("_")},
                         indent=2), flush=True)
        print("--- hydrogen burning, long integration ---", flush=True)
        hydrogen_long = hydrogen_burning(base, args.h_t9, args.h_rho,
                                         args.h_t_end_long, args.steps,
                                         args.h_zmax, args.h_amax, 1e-6, 1e-12)
        print(json.dumps({k: v for k, v in hydrogen_long.items()
                          if not k.startswith("_")}, indent=2), flush=True)

    figures = write_figures(silicon, hydrogen, outdir)

    (outdir / "table_nse_comparison.tex").write_text(latex_nse_table(silicon))
    (outdir / "table_solver_comparison.tex").write_text(latex_solver_table(solvers))
    (outdir / "table_tolerance_study.tex").write_text(latex_tolerance_table(tolerances))

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
        "silicon_burning": {k: v for k, v in silicon.items() if not k.startswith("_")},
        "detailed_balance": detailed,
        "solver_comparison": solvers,
        "tolerance_study": tolerances,
        "figures": figures,
    }
    if hydrogen is not None:
        payload["hydrogen_burning"] = {k: v for k, v in hydrogen.items()
                                       if not k.startswith("_")}
    if hydrogen_long is not None:
        payload["hydrogen_burning_long"] = {k: v for k, v in hydrogen_long.items()
                                            if not k.startswith("_")}
    (outdir / "demonstration.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote demonstration.json, LaTeX tables, and {len(figures)} figures "
          f"to {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
