# NucNetPy

**NucNetPy** is a pure-Python nuclear reaction-network package designed to reproduce the main scientific workflows of the original C/C++ NucNet Tools ecosystem while using a Python-native interface.

The package reads JINA/libnucnet-style XML databases, builds nuclear reaction networks, handles one-zone and multi-zone abundance data, evaluates reaction rates and flows, and provides solver tools for single-zone network evolution.

This version is focused on direct use with real JINA XML files such as:

```text
nuclides.xml
reaction_data.xml
zone.xml
output.xml
```

It does **not** call or wrap the original C++ code.

---

## Main features

- Pure-Python package structure
- JINA/libnucnet XML input support
- Separate nuclide and reaction XML loading
- Zone XML loading for initial abundances and thermodynamic properties
- Combined network XML export
- Nuclide data handling:
  - name
  - charge number `Z`
  - mass number `A`
  - mass excess
  - spin
  - partition-function tables
- Reaction data handling:
  - reactants and products
  - `non_smoker_fit` rates
  - nested `<fit>` blocks
  - `a1`--`a8` coefficients
  - `single_rate` weak/decay reactions
  - special particles such as `gamma`, `electron`, `positron`, `neutrino_e`, and `anti-neutrino_e`
- Mass-fraction and abundance conversion
- Reaction-rate evaluation
- Reaction-flow calculation
- Network conservation checks
- NSE/equilibrium helper routines
- Screening and weak-rate interfaces
- One-zone evolution with SciPy solvers
- Jupyter notebook tutorials
- Command-line tools through `nucnetpy`

---



## Installation

Clone the repository:

```bash
git clone https://github.com/mtalafha90/nucnetpy.git
cd nucnetpy
```

Install in editable mode:

```bash
python -m pip install -e .
```

For development and optional HDF5 support:

```bash
python -m pip install -e .[dev,hdf5]
```

Run the tests:

```bash
pytest -q
```

---

## Required Python packages

Core dependencies:

```text
numpy
scipy
```

Optional dependencies:

```text
pytest       # development/testing
h5py         # HDF5 support
matplotlib   # notebooks/plots
pandas       # notebooks/tables
notebook     # Jupyter tutorials
```

Install optional notebook tools with:

```bash
python -m pip install notebook matplotlib pandas
```

---

## Basic JINA XML workflow

The usual workflow starts from two database files:

```text
nuclides.xml
reaction_data.xml
```

and optionally one zone file:

```text
zone.xml
```

### Summarize the database

```bash
nucnetpy jina-summary nuclides.xml reaction_data.xml
```

With a zone file:

```bash
nucnetpy jina-summary nuclides.xml reaction_data.xml --zones-xml "zone(1).xml"
```

### Combine separate JINA XML files

Many command-line tools operate on one combined XML file. Create it with:

```bash
nucnetpy jina-combine nuclides.xml reaction_data.xml combined_network.xml --zones-xml "zone(1).xml"
```

This produces:

```text
combined_network.xml
```

which contains the nuclides, reactions, and zone data together.

---

## Python example: load JINA XML files

```python
from nucnetpy import read_jina_xml

net = read_jina_xml(
    "nuclides.xml",
    "reaction_data.xml",
    zones_xml="zone(1).xml",
)

print("Number of species:", len(net.species))
print("Number of reactions:", len(net.reactions.reactions))
print("Number of zones:", len(net.zones))

zone = net.zone(0)
print(zone.properties)
```

---

## Python example: inspect the initial zone

```python
from nucnetpy import read_jina_xml

net = read_jina_xml(
    "nuclides.xml",
    "reaction_data.xml",
    zones_xml="zone(1).xml",
)

zone = net.zone(0)

print("Initial properties")
for key, value in zone.properties.items():
    print(f"{key:12s} = {value}")

print("\nInitial mass fractions")
for name, y in sorted(zone.abundances.items()):
    species = net.species.get(name)
    if species is None:
        continue
    x = species.a * y
    if x > 1e-20:
        print(f"{name:8s} X = {x:.8e}   Y = {y:.8e}")
```

---

## Python example: run a single-zone calculation

This example uses the initial temperature and density stored in the zone file. For the uploaded Nova-like zone file, these are typically:

```text
T9  = 0.20
rho = 1.5e4 g cm^-3
```

Create a file called `run_single_zone.py`:

```python
from nucnetpy import read_jina_xml
from nucnetpy.solver import evolve_zone, time_grid
from nucnetpy.thermo import constant_thermo

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
nuclides_xml = "nuclides.xml"
reactions_xml = "reaction_data.xml"
zone_xml = "zone(1).xml"

# ------------------------------------------------------------
# Load network and initial zone
# ------------------------------------------------------------
net = read_jina_xml(
    nuclides_xml,
    reactions_xml,
    zones_xml=zone_xml,
)

zone = net.zone(0)

# ------------------------------------------------------------
# Read thermodynamic properties from zone file
# ------------------------------------------------------------
t9_0 = float(zone.properties.get("t9_0", 0.20))
rho_0 = float(zone.properties.get("rho_0", 1.5e4))
tend = float(zone.properties.get("tend", 100.0))

thermo = constant_thermo(t9=t9_0, rho=rho_0)

# ------------------------------------------------------------
# Time grid
# ------------------------------------------------------------
times = time_grid(0.0, tend, 200)

# ------------------------------------------------------------
# Evolve one zone
# ------------------------------------------------------------
result = evolve_zone(
    network=net,
    zone=zone,
    times=times,
    thermo=thermo,
    method="bdf",
    rtol=1e-8,
    atol=1e-20,
)

# ------------------------------------------------------------
# Print final mass fractions
# ------------------------------------------------------------
print("Final mass fractions")
for name, y in sorted(result.final_abundances.items()):
    species = net.species.get(name)
    if species is None:
        continue
    x = species.a * y
    if x > 1e-20:
        print(f"{name:8s} X = {x:.8e}   Y = {y:.8e}")
```

Run it:

```bash
python run_single_zone.py
```

---

## Command-line example: run a single zone

First combine the JINA files and the zone file:

```bash
nucnetpy jina-combine nuclides.xml reaction_data.xml combined_network.xml --zones-xml "zone(1).xml"
```

Then run the evolution:

```bash
nucnetpy evolve-zone combined_network.xml \
  --zone-index 0 \
  --t0 0.0 \
  --t1 100.0 \
  --steps 200 \
  --t9 0.20 \
  --rho 1.5e4 \
  --method bdf \
  --min-abundance 1e-20
```

---

## Reaction rates

After creating `combined_network.xml`, evaluate rates at a chosen temperature and density:

```bash
nucnetpy rates combined_network.xml --t9 1.0 --rho 1.0e5 --min-rate 1e-30
```

From Python:

```python
from nucnetpy import read_xml

net = read_xml("combined_network.xml")
rates = net.reactions.rates(t9=1.0, rho=1.0e5)

for reaction, rate in list(rates.items())[:10]:
    print(reaction, rate)
```

---

## Reaction flows

Reaction flows require a zone with abundances:

```bash
nucnetpy flows combined_network.xml \
  --zone-index 0 \
  --t9 0.20 \
  --rho 1.5e4 \
  --min-flow 1e-40
```

---

## Network derivative `ydot`

To inspect the abundance derivative for one zone:

```bash
nucnetpy ydot combined_network.xml \
  --zone-index 0 \
  --t9 0.20 \
  --rho 1.5e4 \
  --min-abs 1e-40
```

---

## Conservation checks

Check A/Z conservation for the reaction network:

```bash
nucnetpy conservation combined_network.xml
```

Validate the full XML:

```bash
nucnetpy validate combined_network.xml
```

---

## NSE example

```bash
nucnetpy nse combined_network.xml \
  --t9 5.0 \
  --rho 1.0e7 \
  --ye 0.5 \
  --min-x 1e-12
```

Python version:

```python
from nucnetpy import read_xml
from nucnetpy.nse import solve_nse

net = read_xml("combined_network.xml")
res = solve_nse(net, t9=5.0, rho=1.0e7, ye=0.5)

print(res.success)
print(res.xsum)
print(res.computed_ye)
```

---

## Reading an output XML file

For a full output XML file such as `output_Nova_exp(1).xml`:

```bash
nucnetpy summary "output_Nova_exp(1).xml"
```

Show the largest mass fractions in the first zone:

```bash
nucnetpy largest-x "output_Nova_exp(1).xml" --zone-index 0 -n 20
```

Track a species history across zones:

```bash
nucnetpy species-history "output_Nova_exp(1).xml" h1
nucnetpy species-history "output_Nova_exp(1).xml" he4
nucnetpy species-history "output_Nova_exp(1).xml" c12
```

Show zone properties:

```bash
nucnetpy zone-properties "output_Nova_exp(1).xml" --zone-index 0
```

---

## Jupyter notebooks

Tutorial notebooks are included in the `notebooks/` directory.

Install notebook dependencies:

```bash
python -m pip install notebook matplotlib pandas
```

Start Jupyter:

```bash
jupyter notebook notebooks
```

Suggested order:

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

---

## Repository structure

```text
nucnetpy/
├── pyproject.toml
├── README.md
├── CONVERSION_MAP.md
├── JINA_XML_WORKFLOW.md
├── PURE_PYTHON_PORT_STATUS.md
├── examples/
├── notebooks/
├── tests/
├── validation/
└── src/
    └── nucnetpy/
        ├── core.py
        ├── species.py
        ├── reactions.py
        ├── solver.py
        ├── thermo.py
        ├── nse.py
        ├── screening.py
        ├── weak.py
        ├── analysis.py
        ├── validation.py
        ├── graph.py
        └── io/
            ├── xml.py
            ├── jina.py
            ├── text.py
            └── hdf5.py
```

---

## Important scientific note

This package is a pure-Python replacement path for NucNet Tools workflows. It is designed to read the same JINA/libnucnet XML data and perform equivalent network-analysis and evolution tasks.

However, exact numerical identity with a specific original C++ NucNet Tools build requires dedicated regression tests. The recommended validation sequence is:

1. XML parsing and round-trip tests.
2. ReacLib rate comparison at fixed `T9` values.
3. Reaction-flow comparison at fixed `T9`, `rho`, and abundance state.
4. `ydot` comparison for one-zone states.
5. Screening-factor comparison.
6. Weak-rate table interpolation comparison.
7. NSE comparison.
8. Full single-zone evolution comparison with identical timestep controls and tolerances.

Until those tests are complete, the package should be described as a **pure-Python replacement under validation**, not a guaranteed bitwise numerical clone.

---

## Typical development workflow

```bash
python -m pip install -e .[dev,hdf5]
pytest -q
nucnetpy jina-summary nuclides.xml reaction_data.xml --zones-xml "zone(1).xml"
nucnetpy jina-combine nuclides.xml reaction_data.xml combined_network.xml --zones-xml "zone(1).xml"
nucnetpy validate combined_network.xml
nucnetpy evolve-zone combined_network.xml --t9 0.20 --rho 1.5e4 --t1 100 --steps 200 --method bdf
```

---

## License

This is a research/educational Python conversion project. Before public redistribution, check the license terms of the original NucNet Tools project and any JINA/libnucnet data files used with this package.

---

## Citation and acknowledgement

If this package is used in a scientific project, acknowledge the original NucNet Tools/libnucnet ecosystem and the JINA reaction-rate database used to provide the nuclear data.
