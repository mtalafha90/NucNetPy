"""Compare the NucNetPy results with the numbers the paper quotes.

Writes ``results/comparison.json`` and ``results/comparison.md``:

* the final mass fraction of every species the paper quotes, beside the
  paper's value, its ratio to it, and the value the paper quotes from
  Iliadis (2015) where it gives one;
* the ten most abundant nuclides at the end of each stage;
* the ratio R_Al = X(26Al)/X(27Al) averaged over each advanced stage.

The paper does not say how it averaged R_Al.  Three averages are reported:
over time, ``(1/T) * integral R dt``; over the logarithm of time, which is
what a plain mean over output steps that grow geometrically gives; and the
value at the end of the stage.

    python summarise.py
    python summarise.py --results results_screened
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from common import (ADVANCED, ILIADIS_FINAL, PAPER_FINAL, PAPER_RAL_MEAN,
                    PAPER_RAL_OVERALL, RESULTS, STAGES, STAGE_BY_KEY, UNITS)


def mass_number(name: str) -> int:
    return 1 if name == "n" else int(re.fullmatch(r"[a-z]+(\d+)", name).group(1))


def pretty(name: str) -> str:
    """'al26' -> '26Al'."""
    if name == "n":
        return "n"
    m = re.fullmatch(r"([a-z]+)(\d+)", name)
    return f"{m.group(2)}{m.group(1).capitalize()}"


def final_mass_fractions(results: Path, key: str):
    d = np.load(results / f"{key}.npz")
    return {str(s): float(y) * mass_number(str(s))
            for s, y in zip(d["final_species"], d["final_y"])}


def ratio_history(results: Path, key: str):
    d = np.load(results / f"{key}.npz")
    names = [str(s) for s in d["species"]]
    y = d["y"].astype(float)
    t = d["time"]
    al26 = y[:, names.index("al26")] * 26 if "al26" in names else np.zeros(len(t))
    al27 = y[:, names.index("al27")] * 27
    return t, al26, al27


def ral_averages(t, al26, al27, stage):
    r = np.where(al27 > 0, al26 / np.where(al27 > 0, al27, 1.0), 0.0)
    time_avg = float(np.trapezoid(r, t) / (t[-1] - t[0]))
    logt = np.log(t[1:])              # the log-spaced output grid, t > 0
    log_avg = float(np.trapezoid(r[1:], logt) / (logt[-1] - logt[0]))
    return {"time_average": time_avg, "log_time_average": log_avg,
            "final": float(r[-1]), "maximum": float(r.max()),
            "integral_seconds": float(np.trapezoid(r, t))}


def fmt(v: float) -> str:
    if v == 0:
        return "0"
    if 1e-3 <= abs(v) < 10:
        return f"{v:.4f}".rstrip("0").rstrip(".") if abs(v) >= 0.01 else f"{v:.2e}"
    return f"{v:.2e}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=RESULTS)
    args = ap.parse_args()

    out = {"final": {}, "top": {}, "ral": {}, "run": {}}
    summary_path = args.results / "run_summary.json"
    if summary_path.exists():
        out["run"] = json.loads(summary_path.read_text())

    md = ["# NucNetPy results beside the paper's values", ""]
    md += ["Final mass fractions. The ratio column is NucNetPy divided by the paper.", ""]
    md += ["| Stage | Nuclide | Paper | NucNetPy | Ratio | Iliadis (2015), as quoted |",
           "|---|---|---|---|---|---|"]
    for stage in STAGES:
        if not (args.results / f"{stage.key}.npz").exists():
            continue
        x = final_mass_fractions(args.results, stage.key)
        rows = {}
        for name, paper in PAPER_FINAL.get(stage.key, {}).items():
            ours = x.get(name, 0.0)
            ili = ILIADIS_FINAL.get(stage.key, {}).get(name)
            rows[name] = {"paper": paper, "nucnetpy": ours, "ratio": ours / paper,
                          "iliadis": ili}
            md.append(f"| {stage.title} | {pretty(name)} | {fmt(paper)} | {fmt(ours)} | "
                      f"{ours / paper:.2f} | {fmt(ili) if ili is not None else ''} |")
        for name, ili in ILIADIS_FINAL.get(stage.key, {}).items():
            if name not in rows:
                ours = x.get(name, 0.0)
                rows[name] = {"paper": None, "nucnetpy": ours, "ratio": None, "iliadis": ili}
                md.append(f"| {stage.title} | {pretty(name)} | | {fmt(ours)} | | {fmt(ili)} |")
        out["final"][stage.key] = rows
        top = sorted(x.items(), key=lambda kv: -kv[1])[:10]
        out["top"][stage.key] = dict(top)
        out["final"][stage.key]["_sum_x"] = sum(x.values())

    md += ["", "Ten most abundant nuclides at the end of each stage (mass fractions).", ""]
    for key, top in out["top"].items():
        md.append(f"* **{STAGE_BY_KEY[key].title}**: " +
                  ", ".join(f"{pretty(k)} {fmt(v)}" for k, v in top.items()))

    if all((args.results / f"{k}.npz").exists() for k in ADVANCED):
        md += ["", "R_Al = X(26Al)/X(27Al) averaged over each advanced stage.", "",
               "| Stage | Paper | Time average | Log-time average | Final value |",
               "|---|---|---|---|---|"]
        integral = duration = 0.0
        for key in ADVANCED:
            stage = STAGE_BY_KEY[key]
            t, al26, al27 = ratio_history(args.results, key)
            avg = ral_averages(t, al26, al27, stage)
            out["ral"][key] = avg
            integral += avg["integral_seconds"]
            duration += t[-1] - t[0]
            md.append(f"| {stage.title} | {PAPER_RAL_MEAN[key]:.2e} | "
                      f"{avg['time_average']:.2e} | {avg['log_time_average']:.2e} | "
                      f"{avg['final']:.2e} |")
        out["ral"]["advanced_time_average"] = integral / duration
        md += ["", f"Averaged over the time of all four advanced stages: "
                   f"{integral / duration:.2e} (the paper's abstract gives "
                   f"{PAPER_RAL_OVERALL:.2e} for the star's lifetime, without saying how "
                   "it was formed)."]

    (args.results / "comparison.json").write_text(json.dumps(out, indent=1) + "\n")
    (args.results / "comparison.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
