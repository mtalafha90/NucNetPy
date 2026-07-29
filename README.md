# NucNetPy

[![CI](https://github.com/mtalafha90/NucNetPy/actions/workflows/ci.yml/badge.svg)](https://github.com/mtalafha90/NucNetPy/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20756798.svg)](https://doi.org/10.5281/zenodo.20756798)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

**NucNetPy** is a pure-Python nuclear reaction-network package that reproduces the
main scientific workflows of the original C/C++ [NucNet Tools](https://sourceforge.net/projects/nucnet-tools/)
ecosystem with a Python-native interface.

It reads JINA/libnucnet-style XML databases, builds nuclear reaction networks,
handles one-zone and multi-zone abundance data, evaluates reaction rates and
flows, solves nuclear statistical equilibrium (NSE), and integrates single-zone
network evolution with SciPy stiff solvers.

NucNetPy does **not** call or wrap the original C++ code — it is a clean
reimplementation in Python on top of NumPy and SciPy.

---

## Features

- **XML I/O** — flexible JINA/libnucnet reader supporting separate nuclide and
  reaction files, combined network files, zone files, partition-function tables,
  `non_smoker_fit`/ReacLib rate fits, tabular rates, and `single_rate`
  weak/decay reactions. Round-trips through a compact `nucnetpy` XML format.
- **Nuclear data** — species with `Z`, `A`, mass excess, spin, and
  temperature-dependent partition functions; mass-fraction/abundance conversion;
  special particles (`gamma`, `electron`, `positron`, `neutrino_e`, …).
- **Reactions** — ReacLib seven-parameter rates, tabular and constant rates,
  statistical factors, stoichiometry, reaction flows, and A/Z conservation
  checks.
- **Evolution** — one-zone and multi-zone integration with SciPy BDF/Radau/LSODA
  (plus fixed-step RK4 / implicit Euler fallbacks), an **analytic Jacobian** of
  the reaction flows (one sweep over the reaction list instead of `N+1`
  right-hand-side evaluations — ~275x faster at 100 species), positivity
  projection, and screening / weak-rate hooks. Choose the Jacobian with
  `evolve_zone(..., jac_mode="analytic"|"numerical"|"sparsity")`.
- **NSE & QSE** — statistical weights use `(2J+1)·G(T)`, matching the libnucnet
  convention in which the tabulated partition function is normalised to one at
  low temperature; species lacking nuclear data are excluded from equilibrium
  solves. Robust nuclear statistical equilibrium solve for `sum(A·Y)=1`
  and `sum(Z·Y)=Ye` using a numerically stable log-sum-exp formulation;
  constrained cluster equilibria (`solve_qse`, the libnuceq cluster workflow);
  optional Bravo & García-Senz **Coulomb corrections** ported from the NucNet
  Tools C++ source.
- **Thermodynamic consistency** — `consistent_reverse_network()` rebuilds every
  reverse rate from detailed balance instead of the library's independently
  fitted value. Rate libraries fit the two directions separately, so a network
  built from them relaxes to a state that is *not* the NSE of the same nuclear
  data: for the alpha chain at T9=5 the fitted reverse rates are off by a median
  factor 1.10, displacing equilibrium mass fractions by a median 4.5%. Rebuilt
  from detailed balance, the network reproduces `solve_nse` to **4.5e-6**.
- **Screening** — `SkyNetScreening` assigns each charge a Coulomb chemical
  potential mu(Z)/kT, blending weak, intermediate and strong regimes. Because it
  is a chemical potential it applies consistently to forward and reverse rates,
  so screening does not break detailed balance. Reproduces the Salpeter weak
  limit to 0.02%.
- **Detailed balance** — reverse reaction rates from the forward rate, masses,
  and partition functions; forward/reverse/net flows that vanish at NSE;
  tabulated photodisintegration partners for (n,γ)-(γ,n) studies.
- **Physics helpers** — electron screening (weak/intermediate), 2-D weak-rate
  tables, decays and fission channels, hydrodynamic trajectories, neutrino
  rates, and thermodynamic utilities.
- **Analysis & validation** — largest mass fractions, element abundances,
  abundance moments, energy and **entropy generation rates**, charge-changing
  flows (dYe/dt), species timescales, s-process neutron exposure, integrated
  currents, separation energies, regression comparisons, and Graphviz DOT
  export.
- **Golden-output regression suite** — numerical-identity tests
  (`tests/test_golden_identity.py`) that pin every rate, flow, and trajectory
  against frozen snapshots, ready to be repointed at outputs of an original
  C++ NucNet Tools build.
- **Tooling** — a `nucnetpy` command-line interface and Jupyter tutorial
  notebooks.

---

## Installation

```bash
git clone https://github.com/mtalafha90/nucnetpy.git
cd nucnetpy
python -m pip install -e .
```

Optional extras:

```bash
python -m pip install -e ".[dev,hdf5,plot]"   # tests, HDF5, matplotlib/networkx
python -m pip install -e ".[notebook]"        # jupyter, pandas, matplotlib
```

Run the test suite:

```bash
pytest -q
```

**Requirements:** Python ≥ 3.9, `numpy ≥ 1.22`, `scipy ≥ 1.9`. Optional:
`h5py` (HDF5), `matplotlib` + `networkx` (plots/notebooks), `pytest` (tests).

---

## Quick start (no input files needed)

Build a small network in code, then evolve a single zone:

```python
from nucnetpy import Network, Species, Reaction, RateFit, Zone
from nucnetpy import evolve_zone, time_grid, constant_thermo

net = Network()
# Give species their measured mass excesses and ground-state spins: equilibrium,
# detailed balance and energy release are all derived from them.
net.add_species(Species("he4", z=2, a=4, mass_excess=2.425, spin=0.0))
net.add_species(Species("c12", z=6, a=12, mass_excess=0.0, spin=0.0))
net.add_species(Species("o16", z=8, a=16, mass_excess=-4.737, spin=0.0))

# triple-alpha and 12C(a,g)16O, with toy ReacLib coefficients
net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                      rate_fits=[RateFit([10, 0, 0, 0, 0, 0, 0])], q_value=7.275))
net.reactions.add(Reaction.from_names(["c12", "he4"], ["o16"],
                                      rate_fits=[RateFit([10, 0, 0, 0, 0, 0, 0])], q_value=7.162))

zone = Zone(abundances={"he4": 0.25})        # Y(4He) = 0.25  ->  X = 1.0
net.add_zone(zone)

result = evolve_zone(net, zone, time_grid(0, 1e-2, 50),
                     thermo=constant_thermo(t9=2.0, rho=1.0e5), method="bdf")

# Always check this.  It is False when the solver failed, and also when
# positivity projection had to invent enough baryon number to break
# conservation -- a stiff network can otherwise return a plausible-looking
# composition that is meaningless.
assert result.success, result.message

for name, y in sorted(result.final_abundances.items()):
    a = net.species[name].a
    print(f"{name:5s} Y={y:.6e}  X={a * y:.6e}")
```

> Note: `constant_thermo`, `evolve_zone`, and `time_grid` are exported at the top
> level of `nucnetpy` (they live in `nucnetpy.solver`).

---

## Working with JINA / libnucnet XML

The typical workflow starts from a nuclide file and a reaction file, with an
optional zone file:

```text
nuclides.xml   reaction_data.xml   zone.xml
```

### Command line

```bash
# Summarize a JINA database (species / reactions / zones / validity)
nucnetpy jina-summary nuclides.xml reaction_data.xml --zones-xml zone.xml

# Combine separate files into one network file most commands consume
nucnetpy jina-combine nuclides.xml reaction_data.xml combined_network.xml --zones-xml zone.xml
```

### Python

```python
from nucnetpy import read_jina_xml

net = read_jina_xml("nuclides.xml", "reaction_data.xml", zones_xml="zone.xml")
print("species:", len(net.species))
print("reactions:", len(net.reactions.reactions))
print("zones:", len(net.zones))

zone = net.zone(0)
for name, y in sorted(zone.abundances.items()):
    sp = net.species.get(name)
    if sp and sp.a * y > 1e-20:
        print(f"{name:8s} X={sp.a * y:.6e}  Y={y:.6e}")
```

---

## Two things that catch people out

Neither announces itself with an error message, so both are worth knowing before
your first production run.

### Keep `gamma` when you cut a network

`limit_network` retains only reactions whose participants *all* survive the cut,
and `select_species` returns nuclides. Drop the photon and every
photodisintegration goes with it — silently — and with it the network's ability
to reach equilibrium at all:

```python
from nucnetpy.network_limiter import select_species, limit_network

limit_network(net, select_species(net, zmax=20, amax=44) + ["gamma"])
```

Photons and leptons are reaction *participants*, not nuclides. They balance
charge and lepton number in a reaction record, but they carry no baryon number,
never enter the abundance vector, and are excluded from the abundance product
and the reactant order of a flow. A photodisintegration is a one-body process.

### A network relaxes to the equilibrium of its rates

Rate libraries fit the forward and reverse directions of a reaction separately,
so their ratio is not exactly the equilibrium constant implied by the masses
shipped alongside them. A network built from them settles somewhere that is not
the NSE of the same nuclear data. On the JINA alpha chain at `T9 = 5` the
reverse fits are off by a median factor of 1.10, displacing equilibrium mass
fractions by a median of 4.5 per cent.

`consistent_reverse_network` rebuilds every reverse rate from detailed balance,
which brings the integrated composition and the independent NSE solve into
agreement at the level of a few parts per million:

```python
from nucnetpy import consistent_reverse_network

consistent = consistent_reverse_network(net)   # tabulate=True to keep it writable to XML
```

It is an option rather than a default: deriving reverse rates from mass
differences enforces consistency, but discards whatever measurement went into
the library's reverse fit. Notebook 09 measures both sides of the trade-off.

---

## Performance: the Jacobian

The reaction flows are monomials in the abundances, so the Jacobian is available
in closed form and costs one sweep over the reaction list, against `N+1`
right-hand-side evaluations for finite differences — about 275 times faster at
one hundred species, and agreeing to the truncation error of the difference
quotient.

```python
evolve_zone(net, zone, times, thermo=..., jac_mode="analytic")   # default
evolve_zone(net, zone, times, thermo=..., jac_mode="numerical")  # finite differences
evolve_zone(net, zone, times, thermo=..., jac_mode="sparsity")   # SciPy estimates it
```

SciPy applies a sparsity pattern only when no Jacobian callable is supplied, so
the pattern and an explicit Jacobian are alternatives, not complements.

---

## Screening and energy release

`SkyNetScreening` assigns each nuclear charge a Coulomb chemical potential
`mu(Z)/kT` and enhances a reaction by the difference between the separated
reactants and the fused compound charge. It extends to more than two charged
reactants and, being a chemical potential, applies consistently to forward and
reverse rates, so screening does not undo detailed balance. It reproduces the
Salpeter weak-screening limit to 0.02 per cent.

```python
from nucnetpy import SkyNetScreening

evolve_zone(net, zone, times, thermo=..., screening=SkyNetScreening(net.species))
```

The nuclear energy generation rate follows from the change in total mass excess,
`eps = -N_A * sum_i (dY_i/dt) * dM_i`, so it needs no Q-values and cannot
disagree with them where a rate library and a nuclide file are inconsistent:

```python
from nucnetpy import nuclear_energy_generation_rate, nuclear_energy_release
```

---

## Nuclear statistical equilibrium (NSE)

NSE solves for the proton/neutron chemical potentials that reproduce
`sum(A·Y) = 1` and `sum(Z·Y) = Ye` at a given temperature and density. It shares
no machinery with the time integration, which is what makes it a genuine check
on a network calculation.

Accurate results require mass excesses on the species, and ground-state spins
where they are known: libnucnet files store the partition function normalised to
one at low temperature, so the statistical weight is `(2J+1) * G(T)` and the
degeneracy has to be supplied separately. Species carrying no nuclear data are
excluded automatically — a placeholder mass excess of zero would make an unbound
nuclide look as tightly bound as the most stable one.

```python
from nucnetpy import Network, Species, solve_nse

net = Network()
for name, mass_excess in [("he4", 2.425), ("si28", -21.49),
                          ("fe56", -60.6), ("ni56", -53.9)]:
    net.add_species(Species.parse(name, mass_excess=mass_excess))

res = solve_nse(net, t9=5.0, rho=1.0e8, ye=0.5)
print("success:", res.success)
print("xsum:", res.xsum, " Ye:", res.computed_ye)   # -> ~1.0 and ~0.5
```

Command-line equivalent:

```bash
nucnetpy nse combined_network.xml --t9 5.0 --rho 1.0e8 --ye 0.5 --min-x 1e-12
```

---

## Command-line interface

After installation the `nucnetpy` command exposes the analysis and conversion
tools. A selection:

| Command | Purpose |
|---|---|
| `summary` | species / reaction / zone counts and validity |
| `print-output` | per-zone properties and abundances |
| `largest-x` | largest mass fractions in a zone |
| `zone-abundances`, `zone-properties` | inspect a single zone |
| `element-abundances` | abundances grouped by element |
| `rates`, `flows`, `ydot` | rate, flow, and derivative evaluation |
| `net-flows` | forward, detailed-balance reverse, and net fluxes |
| `charge-flows` | per-reaction dYe/dt contributions |
| `timescales` | shortest species timescales Y/\|dY/dt\| |
| `conservation`, `validate` | A/Z conservation and network validation |
| `evolve-zone` | integrate one zone in time |
| `nse` | nuclear statistical equilibrium solve (`--coulomb` for plasma corrections) |
| `qse` | constrained cluster equilibrium (libnuceq-style) |
| `energy-generation`, `entropy-generation` | energy and entropy generation rates |
| `net-dot` | export the network as a Graphviz DOT graph |
| `species-history` | track a species across zones |
| `remove-duplicates`, `remove-invalid` | clean a reaction set |
| `export-zone-xml`, `reactions-latex` | export helpers |
| `jina-summary`, `jina-combine` | work with separate JINA XML files |

Run `nucnetpy <command> --help` for the full option list of any command.

---

## Jupyter notebooks

Tutorial notebooks live in `notebooks/` and are best followed in order:

| Notebook | Topic |
|---|---|
| `00_installation_and_first_network` | building a network, evolving it, checking `success` |
| `01_species_zones_and_abundances` | `Species`, `Zone`, `Network`; why mass excess and spin matter |
| `02_xml_read_write_and_cli` | XML round trip and the `nucnetpy` command |
| `03_reaction_rates_flows_and_conservation` | rates, flows, conservation, and how photons enter a reaction |
| `04_one_zone_evolution` | stiff solvers, the analytic Jacobian, tolerances |
| `05_nse_screening_and_weak_rates` | equilibrium, `(2J+1)`, screening, weak rates |
| `06_validation_and_regression_workflow` | conservation, golden files, and their limits |
| `07_using_jina_xml_database` | reading and cutting a JINA database |
| `08_validate_real_jina_files` | a full production database (bring your own) |
| `09_thermodynamic_consistency` | detailed balance vs library reverse rates; energy release |

```bash
python -m pip install -e ".[notebook]"
jupyter lab notebooks/
```

See [`notebooks/README.md`](notebooks/README.md) for details.

---

## Repository structure

```text
nucnetpy/
├── pyproject.toml
├── README.md                    # this file
├── CITATION.cff                 # citation metadata (Zenodo DOI)
├── docs/                        # conversion maps, port status, blog coverage
├── examples/                    # standalone usage scripts
├── notebooks/                   # tutorial notebooks (00 ... 09)
├── tests/                       # pytest suite + golden-output snapshots
├── validation/                  # real-JINA-XML validation & golden generator
└── src/nucnetpy/
    ├── core.py                  # Network, Zone containers
    ├── species.py               # Species parsing and Z/A bookkeeping
    ├── reactions.py             # Reaction, RateFit, ReactionNetwork
    ├── solver.py                # evolution, Jacobians, thermo callables
    ├── nse.py                   # nuclear statistical equilibrium
    ├── qse.py                   # constrained cluster equilibria (libnuceq)
    ├── detailed_balance.py      # reverse rates and net flows
    ├── coulomb.py               # Bravo & Garcia-Senz plasma corrections
    ├── screening.py             # electron-screening factors
    ├── weak.py                  # weak-rate tables
    ├── thermo.py                # thermodynamic helpers
    ├── analysis.py              # flows, timescales, entropy, currents, ...
    ├── validation.py            # validation / regression helpers
    ├── decay.py, hydro.py, neutrino.py, network_limiter.py,
    │   rate_modifiers.py, matrix_solver.py, mathutils.py, graph.py
    ├── cli/                     # the `nucnetpy` command
    └── io/                      # xml.py, jina.py, text.py, hdf5.py
```

Documentation beyond this README lives in [`docs/`](docs/README.md), including
the [blog-workflow coverage map](docs/BLOG_COVERAGE.md) and the
[port status / numerical-identity plan](docs/PURE_PYTHON_PORT_STATUS.md).

---

## Scientific note

NucNetPy is a pure-Python *replacement path* for NucNet Tools workflows: it reads
the same JINA/libnucnet XML data and performs equivalent network analysis and
evolution. Exact bitwise agreement with a specific original C++ build is **not**
guaranteed and requires project-specific regression tests against your own
nuclear data, screening choice, and solver tolerances. A recommended validation
order — with the golden-output framework that implements it — is given in
[`docs/PURE_PYTHON_PORT_STATUS.md`](docs/PURE_PYTHON_PORT_STATUS.md).

The strongest internal evidence the package offers is that two independent
subsystems agree: integrating a silicon-burning trajectory to its stationary
state reproduces an independently solved NSE composition to a few parts per
million, once the reverse rates are made thermodynamically consistent. Golden
files pin behaviour against change, but only agreement between calculations that
share no machinery says anything about correctness. Run `pytest -q` for the
suite, and `validation/demonstration.py` for the full comparison.

---

## License

[GPL-3.0-or-later](LICENSE), matching the original NucNet Tools ecosystem.
Check the license terms of any JINA/libnucnet data files used with this
package before redistributing them.

## Citation

If you use this package in scientific work, please cite it via its archived
release and acknowledge the original NucNet Tools / libnucnet ecosystem and the
JINA reaction-rate database that provides the nuclear data.

**DOI:** [10.5281/zenodo.20756798](https://doi.org/10.5281/zenodo.20756798)

```bibtex
@software{nucnetpy,
  title  = {NucNetPy: a pure-Python nuclear reaction-network package},
  url    = {https://github.com/mtalafha90/nucnetpy},
  doi    = {10.5281/zenodo.20756798},
  year   = {2026}
}
```
