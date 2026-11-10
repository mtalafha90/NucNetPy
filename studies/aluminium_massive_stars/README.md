# Aluminium in massive stars: a NucNetPy reproduction

This directory redoes, with NucNetPy, every calculation and figure of

> A. I. Abdallah and N. M. Ershaidat, "A Study of the Abundance of Aluminium in
> Massive Stars", *Jordan Journal of Physics* **17** (3), 273–282 (2024),
> doi:[10.47011/17.3.2](https://doi.org/10.47011/17.3.2).

The paper used NucNet Tools. It follows the core of a 25 solar-mass star
through six burning stages. Each stage is held at a fixed temperature and
density for a fixed time, and each one starts from the composition that the
previous stage left behind. The paper pays particular attention to the two
aluminium isotopes, ²⁶Al and ²⁷Al, and to the ratio of their mass fractions,
R<sub>Al</sub> = X(²⁶Al)/X(²⁷Al).

Results are added to this file when the calculation finishes.

## What is calculated

The six stages are those of the paper's Table 1.

| Stage | Temperature (GK) | Density (g cm⁻³) | Duration |
|---|---|---|---|
| Hydrogen burning | 0.0381 | 3.81 | 6.7 Myr |
| Helium burning | 0.196 | 7.6 × 10² | 0.839 Myr |
| Carbon burning | 0.9 | 1 × 10⁵ | 3.2 kyr |
| Neon burning | 1.5 | 5 × 10⁶ | 0.891 yr |
| Oxygen burning | 2.2 | 3 × 10⁶ | 0.412 yr |
| Silicon burning | 3.6 | 3 × 10⁷ | 0.733 d |

The paper's Table 1 gives 0.733 yr for silicon burning, but its Figs. 6, 7(d)
and 8(d) all run to 0.733 days, so 0.733 d is used here. A year is taken as
365.25 days.

**Nuclear network.** Every JINA ReacLib rate among the nuclides from hydrogen to
nickel (atomic number Z ≤ 28): 891 nuclides and 10,787 reactions. The rates
come from the ReacLib snapshot of 29 March 2022, as distributed with pynucastro
3.2.0. The paper was submitted on 15 March 2022 and does not say which snapshot
it used; this is the closest one available. Every reverse rate that ReacLib
derives by detailed balance is multiplied by the ratio of partition functions
that ReacLib leaves out, as NucNetPy does by default.

**Nuclear data.** Masses from the 2020 Atomic Mass Evaluation (AME2020), spins
from NUBASE2020, and partition functions from Rauscher and Thielemann (2000),
all as distributed with pynucastro. For the 279 nuclides near the drip lines
that AME2020 does not list, the mass is worked out from a ReacLib Q value.
These masses do not enter the calculation, which uses fitted rates, but a
nuclide file needs one.

**Initial composition.** The solar composition of Lodders (2003, ApJ 591, 1220),
isotope by isotope, as tabulated in the file `sollo03.dat` of the StarFit
package. Hydrogen, helium and the heavier elements are each scaled to the
present-day photospheric totals that the paper quotes (X = 0.7491, Y = 0.2377,
Z = 0.0133). Those totals add up to 1.0001 because of rounding, so the whole
composition is then divided by 1.0001.

**Integration.** Implicit (backward) Euler steps solved by Newton–Raphson
iterations, the scheme NucNet Tools uses, built in `stepper.py` from
NucNetPy's own right-hand side and analytic Jacobian. Each step is sized so
that no mass fraction above 10⁻¹⁵ changes by more than about 2 per cent. The
next section explains why NucNetPy's own integrator, `evolve_zone`, could not
be used. Results are recorded at 300 times per stage, spaced evenly in the
logarithm of time from the start of the paper's time axis. A second run with
steps twice as large (`results_coarse/`) shows how much the step size matters.

**Screening.** The paper does not mention electron screening, so the main
calculation leaves it out. A second run with NucNetPy's SkyNet screening is in
`results_screened/`, to show how much it matters.

## How to reproduce

From the repository root, with NucNetPy installed (`pip install -e .`), and
matplotlib for the figures:

```bash
cd studies/aluminium_massive_stars
python run_stages.py          # the six stages, a few hours in total
python summarise.py           # results/comparison.md and .json
python make_figures.py        # figures/fig1_hydrogen.png ... fig8_ratio.png
```

`python run_stages.py --screening skynet --tag screened` repeats the
calculation with screening, writing to `results_screened/`.

The network and the initial composition are stored in `data/`, so the
calculation needs nothing but NucNetPy. `build_inputs.py` regenerates them; it
needs pynucastro 3.2.0 and the file `sollo03.dat` from StarFit 0.19.1 (see the
instructions at the top of that script).

## Why the stages are not integrated with `evolve_zone`

`nucnetpy.evolve_zone` passes the network to SciPy's BDF integrator, and
NucNetPy's right-hand side sets any negative abundance to zero before using
it. On this 891-nuclide network that combination fails in the very first
stage. As hydrogen runs out, the higher-order BDF formula predicts a slightly
negative hydrogen abundance, the right-hand side treats it as zero and stops
every proton capture, and the step size collapses. With an absolute tolerance
of 10⁻²⁰ the integration took more than two hours to pass 1.8 of the 6.7 Myr
and then stopped, the clipping having created 2.5 × 10⁻³ of the baryon number.
With 10⁻¹⁸, and again with screening, it stalled or failed earlier.

Two simple changes to NucNetPy were tried and rejected. Without the clipping,
hydrogen and helium burning took about six minutes each, but carbon burning
failed with the same collapse of the step size, probably because two small
negative abundances make a positive flow between them, which feeds itself.
With a Jacobian that agrees with the clipping, the integration got as far as
neon burning (its hydrogen-burning result agreed with the unclipped run to
1.5 × 10⁻⁴) and then failed in oxygen burning, after small negative
abundances had built up to 10⁻⁵ of the baryon number.

Implicit Euler does not have this problem. For a decaying abundance it takes
y → y/(1 + λΔt), which can never turn negative, and that is why NucNet Tools
uses it. Its price is first-order accuracy, so the steps must stay small; the
whole calculation takes a few hours in pure Python.

## Two things NucNetPy does not do yet

The calculation needed two small corrections that NucNetPy 1.0.0 does not make
on its own. Both are made by this directory's scripts, not by changes to
NucNetPy.

1. **Electron captures.** ReacLib gives the rates of its three electron
   captures (on ³He, on ⁷Be, and the pep reaction p + p + e⁻) per unit electron
   density, so each must be multiplied by ρY<sub>e</sub>. NucNetPy does not do
   this. `common.RateFactors` supplies the factor through the screening hook
   that NucNetPy evaluates for every rate. The resulting rate coefficients
   agree exactly with pynucastro's.
2. **Which member of a reverse pair gets the partition functions.** NucNetPy
   assumes that ReacLib's derived reverse rate is always the endothermic
   member of a pair. In 168 pairs of this network the forward rate is the
   endothermic one (a capture into a nucleus that is unbound near the drip
   lines), so NucNetPy corrects the capture instead of the
   photodisintegration. Three more pairs are missed because ReacLib rounds
   their large coefficients; one of them is ²⁹Si(γ,p)²⁸Al, whose correction is
   a factor of 1.64 at 3.6 GK. The network file keeps ReacLib's own flag for
   derived rates, and `run_stages.assign_reverse_partition_functions` places
   the factor by that flag.

## Files

| File | What it is |
|---|---|
| `common.py` | The six stages, the species in each figure, the values the paper quotes, and the electron-capture factor |
| `build_inputs.py` | Builds the network and the initial composition in `data/` |
| `stepper.py` | Implicit Euler steps with Newton–Raphson iterations, built from NucNetPy's right-hand side and Jacobian |
| `run_stages.py` | Runs the six stages with NucNetPy |
| `summarise.py` | Compares the results with the paper's values |
| `make_figures.py` | Redraws Figs. 1–8 |
| `data/network_reaclib_default2_20220329_z28.xml.gz` | The network, in the libnucnet XML format NucNetPy reads |
| `data/initial_composition_lodders2003.json` | The initial mass fractions |
| `data/provenance.json` | Versions, checksums and counts for both inputs |

## Sources and licences

* ReacLib: R. H. Cyburt et al., ApJS 189, 240 (2010); snapshot
  `reaclib_default2_20220329` from pynucastro 3.2.0 (BSD-3-Clause).
* AME2020: M. Wang et al., Chinese Physics C 45, 030003 (2021). NUBASE2020:
  F. G. Kondev et al., Chinese Physics C 45, 030001 (2021).
* Partition functions: T. Rauscher and F.-K. Thielemann, Atomic Data and
  Nuclear Data Tables 75, 1 (2000).
* Solar composition: K. Lodders, ApJ 591, 1220 (2003), from `sollo03.dat` in
  StarFit 0.19.1 (GPL-3.0, like NucNetPy).
