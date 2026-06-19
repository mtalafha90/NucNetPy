# NucNetPy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20756798.svg)](https://doi.org/10.5281/zenodo.20756798)

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
  (plus fixed-step RK4 / implicit Euler fallbacks), analytic sparsity pattern and
  numerical Jacobian, positivity projection, and screening / weak-rate hooks.
- **NSE** — robust nuclear statistical equilibrium solve for `sum(A·Y)=1` and
  `sum(Z·Y)=Ye` using a numerically stable log-sum-exp formulation.
- **Physics helpers** — electron screening (weak/intermediate), 2-D weak-rate
  tables, decays, hydrodynamic trajectories, neutrino rates, and thermodynamic
  utilities.
- **Analysis & validation** — largest mass fractions, element abundances,
  abundance moments, energy generation, regression comparisons, and Graphviz DOT
  export.
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
for name in ["he4", "c12", "o16"]:
    net.add_species(Species.parse(name))

# triple-alpha and 12C(a,g)16O, with toy ReacLib coefficients
net.reactions.add(Reaction.from_names(["he4", "he4", "he4"], ["c12"],
                                      rate_fits=[RateFit([10, 0, 0, 0, 0, 0, 0])], q_value=7.275))
net.reactions.add(Reaction.from_names(["c12", "he4"], ["o16"],
                                      rate_fits=[RateFit([10, 0, 0, 0, 0, 0, 0])], q_value=7.162))

zone = Zone(abundances={"he4": 0.25})        # Y(4He) = 0.25  ->  X = 1.0
net.add_zone(zone)

result = evolve_zone(net, zone, time_grid(0, 1e-2, 50),
                     thermo=constant_thermo(t9=2.0, rho=1.0e5), method="bdf")

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

## Nuclear statistical equilibrium (NSE)

NSE solves for the proton/neutron chemical potentials that reproduce
`sum(A·Y) = 1` and `sum(Z·Y) = Ye` at a given temperature and density. Accurate
results require mass excesses on the species:

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
| `conservation`, `validate` | A/Z conservation and network validation |
| `evolve-zone` | integrate one zone in time |
| `nse` | nuclear statistical equilibrium solve |
| `energy-generation` | nuclear energy generation rate |
| `net-dot` | export the network as a Graphviz DOT graph |
| `species-history` | track a species across zones |
| `remove-duplicates`, `remove-invalid` | clean a reaction set |
| `export-zone-xml`, `reactions-latex` | export helpers |
| `jina-summary`, `jina-combine` | work with separate JINA XML files |

Run `nucnetpy <command> --help` for the full option list of any command.

---

## Jupyter notebooks

Tutorial notebooks live in `notebooks/` and are best followed in order:

```text
00_installation_and_first_network.ipynb
01_species_zones_and_abundances.ipynb
02_xml_read_write_and_cli.ipynb
03_reaction_rates_flows_and_conservation.ipynb
04_one_zone_evolution.ipynb
05_nse_screening_and_weak_rates.ipynb
06_validation_and_regression_workflow.ipynb
07_using_jina_xml_database.ipynb
08_validate_real_jina_files.ipynb
```

```bash
python -m pip install notebook matplotlib
jupyter notebook notebooks
```

---

## Repository structure

```text
nucnetpy/
├── pyproject.toml
├── README.md
├── CONVERSION_MAP.md            # NucNet Tools concept -> nucnetpy mapping
├── JINA_XML_WORKFLOW.md
├── PURE_PYTHON_PORT_STATUS.md
├── examples/                    # standalone usage scripts
├── notebooks/                   # tutorial notebooks
├── tests/                       # pytest suite
├── validation/                  # real-JINA-XML validation script
└── src/nucnetpy/
    ├── core.py                  # Network, Zone containers
    ├── species.py               # Species parsing and Z/A bookkeeping
    ├── reactions.py             # Reaction, RateFit, ReactionNetwork
    ├── solver.py                # evolution, Jacobians, thermo helpers
    ├── nse.py                   # nuclear statistical equilibrium
    ├── screening.py             # electron-screening factors
    ├── weak.py                  # weak-rate tables
    ├── thermo.py                # thermodynamic helpers
    ├── analysis.py              # analysis utilities
    ├── validation.py            # validation / regression helpers
    ├── decay.py, hydro.py, neutrino.py, network_limiter.py,
    │   rate_modifiers.py, matrix_solver.py, mathutils.py, graph.py
    └── io/                      # xml.py, jina.py, text.py, hdf5.py
```

---

## Scientific note

NucNetPy is a pure-Python *replacement path* for NucNet Tools workflows: it reads
the same JINA/libnucnet XML data and performs equivalent network analysis and
evolution. Exact bitwise agreement with a specific original C++ build is **not**
guaranteed and requires project-specific regression tests against your own
nuclear data, screening choice, and solver tolerances. A recommended validation
order is given in [`PURE_PYTHON_PORT_STATUS.md`](PURE_PYTHON_PORT_STATUS.md).

---

## License

Research/educational Python conversion project. Before public redistribution,
check the license terms of the original NucNet Tools project and of any
JINA/libnucnet data files used with this package. See [`LICENSE`](LICENSE).

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
