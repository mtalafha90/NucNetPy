"""Redraw Figs. 1-8 of the paper from the NucNetPy results.

Figs. 1-6 show the mass fractions of the paper's chosen species during each
burning stage, on logarithmic axes.  An open circle at the end of a stage
marks each final mass fraction the paper quotes in its text, so the two
calculations can be compared by eye.  Figs. 7 and 8 show 26Al, 27Al and
their ratio during the four advanced stages, against linear time, as in the
paper.  The figures are written to ``figures/``.

    python make_figures.py
    python make_figures.py --results results_screened --outdir figures_screened
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from common import (ADVANCED, FIGURES, PAPER_FINAL, RESULTS, STAGE_BY_KEY,  # noqa: E402
                    STAGES, UNITS)

# Categorical palette (fixed order) and chart ink, light surface.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
           "#6250d6", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
DASHED = (0, (4.0, 1.6))

#: Lower edge of the mass-fraction axis of Figs. 1-6, as in the paper.
YMIN = {"hydrogen": 1e-10, "helium": 3e-5, "carbon": 1e-5, "neon": 1e-4,
        "oxygen": 1e-6, "silicon": 1e-6}


def nuclide_label(name: str) -> str:
    """'al26' -> '$^{26}$Al', the form the paper uses."""
    if name == "n":
        return "n"
    m = re.fullmatch(r"([a-z]+)(\d+)", name)
    return rf"$^{{{m.group(2)}}}$" + m.group(1).capitalize()


def mass_number(name: str) -> int:
    return 1 if name == "n" else int(re.fullmatch(r"[a-z]+(\d+)", name).group(1))


def load(results: Path, key: str):
    d = np.load(results / f"{key}.npz")
    names = [str(s) for s in d["species"]]
    y = d["y"].astype(float)
    x = {s: y[:, i] * mass_number(s) for i, s in enumerate(names)}
    return d["time"], x


def series(x, name, n):
    return x.get(name, np.zeros(n))


def style(ax):
    ax.set_facecolor(SURFACE)
    for spine in ax.spines.values():
        spine.set_color(AXIS)
        spine.set_linewidth(0.8)
    ax.grid(True, which="major", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK2, labelsize=9, which="both")
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def line_style(i):
    """Colour in fixed palette order; past eight series, dashed lines carry identity too."""
    return PALETTE[i % len(PALETTE)], ("-" if i < len(PALETTE) else DASHED)


def stage_figure(stage, results: Path, outdir: Path):
    t, x = load(results, stage.key)
    t_unit = t / UNITS[stage.unit]
    keep = t_unit > 0
    fig, ax = plt.subplots(figsize=(7.6, 4.7))
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    handles = []
    finals = []
    quoted = PAPER_FINAL.get(stage.key, {})
    for i, name in enumerate(stage.species):
        colour, ls = line_style(i)
        xi = series(x, name, len(t))
        ax.plot(t_unit[keep], np.where(xi[keep] > 0, xi[keep], np.nan), color=colour,
                ls=ls, lw=1.8, solid_capstyle="round", solid_joinstyle="round")
        if name in quoted:
            ax.plot([stage.duration], [quoted[name]], ls="none", marker="o", ms=7,
                    mfc=SURFACE, mec=colour, mew=1.6, clip_on=False, zorder=5)
        handles.append(Line2D([], [], color=colour, ls=ls, lw=1.8,
                              label=nuclide_label(name)))
        finals.append(xi[-1])
    # Legend in the order the curves end, top to bottom.
    order = np.argsort(finals)[::-1]
    handles = [handles[i] for i in order]
    if quoted:
        handles.append(Line2D([], [], ls="none", marker="o", ms=7, mfc=SURFACE,
                              mec=INK2, mew=1.4, label="Final value\nquoted by the paper"))
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                    frameon=False, fontsize=9, handlelength=2.2, labelcolor=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(stage.first_output, stage.duration * 1.6)
    ax.set_ylim(YMIN[stage.key], 2.0)
    ax.set_xlabel(f"Time ({stage.unit})")
    ax.set_ylabel("Mass fraction")
    ax.set_title(f"Fig. {stage.figure}. {stage.title}: $T_9$ = {stage.t9:g}, "
                 rf"$\rho$ = {_sci(stage.rho)} g cm$^{{-3}}$, {stage.duration:g} {stage.unit}",
                 loc="left", fontsize=10.5, color=INK)
    fig.tight_layout()
    path = outdir / f"fig{stage.figure}_{stage.key}.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_extra_artists=[leg],
                bbox_inches="tight")
    plt.close(fig)
    return path


def _sci(v: float) -> str:
    if 0.01 <= v < 1000:
        return f"{v:g}"
    e = int(np.floor(np.log10(v)))
    m = v / 10 ** e
    return (rf"{m:g}$\times$10$^{{{e}}}$" if abs(m - 1) > 1e-9 else rf"10$^{{{e}}}$")


def aluminium(results: Path, key: str):
    t, x = load(results, key)
    stage = STAGE_BY_KEY[key]
    return stage, t / UNITS[stage.unit], series(x, "al26", len(t)), series(x, "al27", len(t))


def time_average(t, f):
    """(1/T) * integral of f over the stage, by the trapezoidal rule."""
    return float(np.trapezoid(f, t) / (t[-1] - t[0]))


def figure7(results: Path, outdir: Path):
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.4))
    fig.patch.set_facecolor(SURFACE)
    for ax, key, tag in zip(axes.flat, ADVANCED, "abcd"):
        stage, tu, al26, al27 = aluminium(results, key)
        style(ax)
        for f, colour, name in ((al26, PALETTE[0], "al26"), (al27, PALETTE[1], "al27")):
            ax.plot(tu, np.where(f > 0, f, np.nan), color=colour, lw=1.8)
            ax.annotate(nuclide_label(name), xy=(tu[-1], f[-1]), xytext=(5, 0),
                        textcoords="offset points", va="center", fontsize=9, color=INK)
        ax.set_yscale("log")
        ax.set_xlim(0, stage.duration * 1.12)
        ax.set_xlabel(f"Time ({stage.unit})")
        ax.set_ylabel("Mass fraction")
        ax.set_title(f"({tag}) {stage.title}", loc="left", fontsize=10, color=INK)
    handles = [Line2D([], [], color=PALETTE[0], lw=1.8, label=nuclide_label("al26")),
               Line2D([], [], color=PALETTE[1], lw=1.8, label=nuclide_label("al27"))]
    fig.legend(handles=handles, loc="upper right", ncol=2, frameon=False, fontsize=9.5,
               labelcolor=INK)
    fig.suptitle("Fig. 7. Mass fractions of $^{26}$Al and $^{27}$Al in the advanced stages",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = outdir / "fig7_al26_al27.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def figure8(results: Path, outdir: Path):
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.4))
    fig.patch.set_facecolor(SURFACE)
    for ax, key, tag in zip(axes.flat, ADVANCED, "abcd"):
        stage, tu, al26, al27 = aluminium(results, key)
        style(ax)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(al27 > 0, al26 / al27, np.nan)
        ax.plot(tu, ratio, color=PALETTE[5], lw=1.8)
        mean = time_average(tu, np.nan_to_num(ratio))
        ax.axhline(mean, color=MUTED, lw=0.8)
        ax.annotate(f"time average {mean:.2e}", xy=(stage.duration, mean), xytext=(0, 4),
                    textcoords="offset points", ha="right", va="bottom", fontsize=8.5,
                    color=INK2)
        ax.set_yscale("log")
        ax.set_xlim(0, stage.duration)
        ax.set_xlabel(f"Time ({stage.unit})")
        ax.set_ylabel(r"$R_\mathrm{Al}$ = $X(^{26}$Al$)/X(^{27}$Al$)$")
        ax.set_title(f"({tag}) {stage.title}", loc="left", fontsize=10, color=INK)
    fig.suptitle(r"Fig. 8. The ratio $R_\mathrm{Al}$ of the $^{26}$Al and $^{27}$Al "
                 "mass fractions", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = outdir / "fig8_ratio.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--outdir", type=Path, default=FIGURES)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.default": "regular"})
    for stage in STAGES:
        if (args.results / f"{stage.key}.npz").exists():
            print("wrote", stage_figure(stage, args.results, args.outdir).name)
    if all((args.results / f"{k}.npz").exists() for k in ADVANCED):
        print("wrote", figure7(args.results, args.outdir).name)
        print("wrote", figure8(args.results, args.outdir).name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
