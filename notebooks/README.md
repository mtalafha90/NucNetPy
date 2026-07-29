# NucNetPy tutorial notebooks

A guided tour of the package, from building a three-species network by hand to
reading a production JINA/libnucnet database and checking a calculation against
an independently solved equilibrium.

## Running them

```bash
python -m pip install -e ".[notebook]"   # jupyter, pandas, matplotlib
jupyter lab notebooks/
```

The `notebook` extra exists because these tutorials import `pandas` and
`matplotlib`, which the package itself does not. NucNetPy needs only NumPy and
SciPy.

Every notebook is self-contained and runs in seconds, with one exception:
notebook 08 needs a production database, which is not redistributed here. Its
cells skip cleanly when the files are absent, so it can still be read and run.

## Contents

| Notebook | What it covers |
|---|---|
| `00_installation_and_first_network` | Installing, building a network, evolving it, checking `success` |
| `01_species_zones_and_abundances` | `Species`, `Zone`, `Network`; why mass excess and spin matter |
| `02_xml_read_write_and_cli` | Compact XML, the round-trip fixed point, the `nucnetpy` command |
| `03_reaction_rates_flows_and_conservation` | Rate forms, flows, `dY/dt`, conservation — **and how photons enter a reaction** |
| `04_one_zone_evolution` | Stiff solvers, the analytic Jacobian, tolerances, prescribed trajectories |
| `05_nse_screening_and_weak_rates` | Equilibrium, the `(2J+1)` statistical weight, screening models, weak-rate tables |
| `06_validation_and_regression_workflow` | Conservation along a trajectory, golden files, and what they do *not* prove |
| `07_using_jina_xml_database` | Reading JINA files, validating them, cutting a network down to size |
| `08_validate_real_jina_files` | Working with a full production database (bring your own) |
| `09_thermodynamic_consistency` | Detailed balance against library reverse rates; energy release |

## Using notebook 08 with your own data

Copy the files into `notebooks/data/` with these names:

```text
nuclides.xml         nuclide records: mass excesses, spins, partition functions
reaction_data.xml    reaction records with ReacLib / non-smoker fits
zone.xml             an initial composition (optional)
```

The notebook then reports parser counts, mass-fraction normalisation, reaction
conservation, species lacking nuclear data, and sample rates.

## Two things that catch people out

Both have their own section in the notebooks, but they are worth stating here
because neither announces itself with an error message.

**Keep `gamma` when you cut a network.** `limit_network` retains only reactions
whose participants all survive the cut, and `select_species` returns nuclides.
Drop the photon and every photodisintegration goes with it — silently, and with
it the network's ability to reach equilibrium:

```python
limit_network(net, select_species(net, zmax=20, amax=44) + ["gamma"])
```

**A network relaxes to the equilibrium of the rates it is given**, which is not
necessarily the equilibrium of the nuclear data, because rate libraries fit the
forward and reverse directions separately. Notebook 09 measures the gap and
shows how to close it with `consistent_reverse_network`.

## Data

`data/` holds the small fixtures the notebooks create or read. Notebooks 02, 06
and 07 write files there as they run.
