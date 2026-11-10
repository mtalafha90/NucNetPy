"""Compare the NucNetPy results with the numbers the paper quotes.

Writes ``comparison.json`` and ``comparison.md`` into the results directory:

* the final mass fraction of every species the paper quotes, beside the
  paper's value, their ratio, the value with plasma screening (when
  ``results_screened/`` exists) and the value the paper quotes from Iliadis
  (2015) where it gives one;
* silicon burning run for 0.733 years instead of days (when
  ``results_silicon_year/`` exists);
* the ten most abundant nuclides at the end of each stage;
* the ratio R_Al = X(26Al)/X(27Al) averaged over each advanced stage;
* with ``--check DIR``, the largest difference from the final compositions of
  another run, for example one made with different solver settings.

The paper does not say how it averaged R_Al.  Three measures are reported:
the average over time, ``(1/T) * integral of R dt``; the average over the
logarithm of time, which is what a plain mean over output steps that grow
geometrically amounts to; and the value at the end of the stage.

    python summarise.py
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

from common import (ADVANCED, HERE, ILIADIS_FINAL, PAPER_FINAL, PAPER_RAL_MEAN,
                    PAPER_RAL_OVERALL, RESULTS, STAGES, STAGE_BY_KEY)

SCREENED = HERE / "results_screened"
SILICON_YEAR = HERE / "results_silicon_year"


def mass_number(name: str) -> int:
    return 1 if name == "n" else int(re.fullmatch(r"[a-z]+(\d+)", name).group(1))


def pretty(name: str) -> str:
    """'al26' -> '26Al'."""
    if name == "n":
        return "n"
    m = re.fullmatch(r"([a-z]+)(\d+)", name)
    return f"{m.group(2)}{m.group(1).capitalize()}"


def has(results: Path, key: str) -> bool:
    return (results / f"{key}.npz").exists()


def final_mass_fractions(results: Path, key: str):
    d = np.load(results / f"{key}.npz")
    return {str(s): float(y) * mass_number(str(s))
            for s, y in zip(d["final_species"], d["final_y"])}


def ral_averages(results: Path, key: str):
    d = np.load(results / f"{key}.npz")
    names = [str(s) for s in d["species"]]
    y = d["y"].astype(float)
    t = d["time"]
    al26 = y[:, names.index("al26")] * 26 if "al26" in names else np.zeros(len(t))
    al27 = y[:, names.index("al27")] * 27
    r = np.divide(al26, al27, out=np.zeros_like(al26), where=al27 > 0)
    logt = np.log(t[1:])              # the logarithmic output grid, t > 0
    return {"time_average": float(np.trapezoid(r, t) / (t[-1] - t[0])),
            "log_time_average": float(np.trapezoid(r[1:], logt) / (logt[-1] - logt[0])),
            "final": float(r[-1]), "maximum": float(r.max()),
            "integral_seconds": float(np.trapezoid(r, t)),
            "seconds": float(t[-1] - t[0])}


def fmt(v) -> str:
    if v is None:
        return ""
    if v == 0:
        return "0"
    if 0.01 <= abs(v) < 10:
        return f"{v:.4f}".rstrip("0").rstrip(".")
    return f"{v:.2e}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--check", type=Path, nargs="*", default=[],
                    help="other results directories to compare the final compositions with")
    args = ap.parse_args()
    res = args.results

    out = {"final": {}, "top": {}, "ral": {}, "tolerance_check": {}, "silicon_year": {}}
    summary_path = res / "run_summary.json"
    if summary_path.exists():
        out["run"] = json.loads(summary_path.read_text())

    md = ["# NucNetPy results beside the paper's values", "",
          "Final mass fractions. *Ratio* is NucNetPy divided by the paper. "
          "*Screened* repeats the calculation with NucNetPy's SkyNet screening.", "",
          "| Stage | Nuclide | Paper | NucNetPy | Ratio | Screened | Iliadis (2015), as the paper quotes it |",
          "|---|---|---|---|---|---|---|"]
    for stage in STAGES:
        if not has(res, stage.key):
            continue
        x = final_mass_fractions(res, stage.key)
        xs = final_mass_fractions(SCREENED, stage.key) if has(SCREENED, stage.key) else {}
        quoted = dict(PAPER_FINAL.get(stage.key, {}))
        for name in ILIADIS_FINAL.get(stage.key, {}):
            quoted.setdefault(name, None)
        rows = {}
        for name, paper in quoted.items():
            ours = x.get(name, 0.0)
            ili = ILIADIS_FINAL.get(stage.key, {}).get(name)
            ratio = ours / paper if paper else None
            rows[name] = {"paper": paper, "nucnetpy": ours, "ratio": ratio,
                          "screened": xs.get(name), "iliadis": ili}
            md.append(f"| {stage.title} | {pretty(name)} | {fmt(paper)} | {fmt(ours)} | "
                      f"{'' if ratio is None else f'{ratio:.2f}'} | {fmt(xs.get(name))} | "
                      f"{fmt(ili)} |")
        out["final"][stage.key] = rows
        out["top"][stage.key] = dict(sorted(x.items(), key=lambda kv: -kv[1])[:10])

    if has(SILICON_YEAR, "silicon") and has(res, "silicon"):
        days = final_mass_fractions(res, "silicon")
        years = final_mass_fractions(SILICON_YEAR, "silicon")
        md += ["", "Silicon burning for 0.733 days (the paper's figures) and for "
                   "0.733 years (its Table 1).", "",
               "| Nuclide | Paper | 0.733 d | 0.733 yr |", "|---|---|---|---|"]
        for name, paper in PAPER_FINAL["silicon"].items():
            out["silicon_year"][name] = {"paper": paper, "days": days.get(name, 0.0),
                                         "years": years.get(name, 0.0)}
            md.append(f"| {pretty(name)} | {fmt(paper)} | {fmt(days.get(name, 0.0))} | "
                      f"{fmt(years.get(name, 0.0))} |")
        out["silicon_year"]["_top_years"] = dict(sorted(years.items(),
                                                        key=lambda kv: -kv[1])[:10])

    md += ["", "Ten most abundant nuclides at the end of each stage (mass fractions).", ""]
    for key, top in out["top"].items():
        md.append(f"* **{STAGE_BY_KEY[key].title}**: " +
                  ", ".join(f"{pretty(k)} {fmt(v)}" for k, v in top.items()))

    if all(has(res, k) for k in ADVANCED):
        screened = all(has(SCREENED, k) for k in ADVANCED)
        md += ["", "R_Al = X(26Al)/X(27Al) averaged over each advanced stage.", "",
               "| Stage | Paper | Time average | Log-time average | Final value | "
               "Time average, screened |", "|---|---|---|---|---|---|"]
        integral = duration = 0.0
        for key in ADVANCED:
            avg = ral_averages(res, key)
            out["ral"][key] = avg
            integral += avg["integral_seconds"]
            duration += avg["seconds"]
            scr = ral_averages(SCREENED, key)["time_average"] if screened else None
            if scr is not None:
                avg["screened_time_average"] = scr
            md.append(f"| {STAGE_BY_KEY[key].title} | {PAPER_RAL_MEAN[key]:.2e} | "
                      f"{avg['time_average']:.2e} | {avg['log_time_average']:.2e} | "
                      f"{avg['final']:.2e} | {fmt(scr)} |")
        out["ral"]["advanced_time_average"] = integral / duration
        md += ["", f"Averaged over the whole time of the four advanced stages: "
                   f"{integral / duration:.2e}. The paper's abstract gives "
                   f"{PAPER_RAL_OVERALL:.2e} over the star's lifetime without saying "
                   "how it was formed."]

    for other in args.check:
        checks = [k for k in (s.key for s in STAGES) if has(other, k) and has(res, k)]
        if not checks:
            continue
        md += ["", f"Compared with `{other.name}`: the largest relative difference in "
                   "any final mass fraction above 1e-12.", ""]
        for key in checks:
            a, b = final_mass_fractions(res, key), final_mass_fractions(other, key)
            worst, where = 0.0, ""
            for name, v in a.items():
                if v > 1e-12:
                    d = abs(b.get(name, 0.0) - v) / v
                    if d > worst:
                        worst, where = d, name
            out["tolerance_check"].setdefault(other.name, {})[key] = {
                "max_relative_difference": worst, "nuclide": where}
            md.append(f"* {STAGE_BY_KEY[key].title}: {worst:.1e} ({pretty(where)})")

    (res / "comparison.json").write_text(json.dumps(out, indent=1) + "\n")
    (res / "comparison.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
