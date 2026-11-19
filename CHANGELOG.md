# Changelog

The changes in each release of NucNetPy. Every release is archived on Zenodo
with its own DOI; [10.5281/zenodo.20756797](https://doi.org/10.5281/zenodo.20756797)
covers all versions and resolves to the latest.

## 1.2.0 (2026-11-19)

### Changed

- The right-hand side and the analytic Jacobian evaluate every reaction at
  once, with NumPy arrays and a sparse stoichiometric matrix, instead of in a
  Python loop over the reactions. On a network of 900 species and 10,800
  reactions one right-hand side takes 0.4 ms instead of 110 ms, and the
  analytic Jacobian 2 ms instead of 60 ms. The results do not change: on the
  test machine every number of the demonstration record and of the archived
  reference case is the same to the last digit.
- `jacobian_sparsity` builds the pattern as a sparse matrix from the start.
- The adaptive implicit Euler integrator factorises its Newton matrix with
  diagonal pivots, which keep the fill-reducing ordering: on a complete
  REACLIB network each factorisation is 2 to 6 times faster, and a hydrogen-rich
  zone at T9 = 2 and rho = 1e6 g/cm^3 integrates for a full second in about
  five minutes. The results are those of threshold pivoting to rounding error.
- `evolve-zone` reports the cost of one right-hand side above 1000 species,
  and suggests settings that get through explosive burning when two
  successive progress reports find an integration crawling. It does not
  repeat settings already in use.
- `evolve-zone` takes the temperature, density, end time and first step from
  the zone's `t9` (or `t9_0`), `rho` (or `rho_0`), `tend` and `dt` properties
  whenever `--t9`, `--rho`, `--t1` or `--dt` is left out, as NucNet Tools
  does, and reports on standard error which values it uses and where they
  came from. It refuses to run without a temperature or a density, which it
  used to take as 1 without a word; the end time still defaults to 1 s.
- `evolve-zone` ends with a short summary instead of the final abundance of
  every species: the verdict, the steps and evaluations, the change in the
  sum of the mass fractions, `Ye` and the ten largest mass fractions. Its
  progress reports give the number of steps and the size of the last one.
- `evolve_zone` calls `progress` after every accepted step, whatever the
  method, instead of on every right-hand-side evaluation for the SciPy
  methods, and `EvolutionResult.steps` counts the steps of every method. The
  results do not change.
- Every command has a one-line description in `nucnetpy --help`.
- The documentation describes the package as it is; the history of each
  release is kept here.

### Added

- `evolve-zone --dt` and `evolve_zone(..., first_step=...)`: the first time
  step of the SciPy and adaptive implicit Euler integrators.
- `Zone.end_time()`, `Zone.time_step()` and `Zone.find_property()`.
- `evolve-zone --max-change` and `--min-mass-fraction`, the step control of
  the adaptive implicit Euler integrator.
- `evolve-zone` saves the whole history of every run, in a new folder named
  after the date and time: `output.xml` holds one zone per output time, with
  its time, `t9` and `rho`, in the libnucnet format of NucNet Tools, which
  wnutils and the `nucnetpy` zone commands read. `--output` names the folder,
  `--no-output` saves nothing, and `--min-abundance` leaves the smallest
  abundances out. A run never overwrites an earlier `output.xml`.
- `write_evolution_xml(network, result, path, thermo=...)` writes that file
  from any evolution.
- `validation/benchmark_performance.py` integrates every network size, with
  BDF and with adaptive implicit Euler, and stops each integration after
  `--evolve-timeout` seconds; `--fd-limit` bounds the size up to which it
  times the finite-difference Jacobian.

## 1.1.1 (2026-11-17)

DOI [10.5281/zenodo.24323777](https://doi.org/10.5281/zenodo.24323777)

### Fixed

- Combined files keep the record of species without nuclear data. A
  placeholder is a nuclide that the reactions name but the nuclide file leaves
  out. `nucnetpy jina-combine` and `write_xml` used to write each one as an
  ordinary nuclide with a mass excess of zero, so an equilibrium solve on the
  combined file included unbound nuclides again. Placeholders are now left out
  of the nuclide records, and both readers recreate and record them.
  **Regenerate any combined file written by an earlier version.**
- `nucnetpy validate` no longer reports photons and leptons as missing
  species. It reports the species without nuclear data in one warning.
- HDF5 files no longer list photons and leptons as species.

### Added

- Options of `nucnetpy evolve-zone` for large networks: `--zmax` and `--amax`
  cut the network, `--rtol` and `--atol` set the tolerances, and `--progress`
  reports progress on standard error (every 30 s by default). Above 400 species
  the command warns before it starts.
- `evolve_zone(..., progress=callable)`.

## 1.1.0 (2026-11-12)

DOI [10.5281/zenodo.24152490](https://doi.org/10.5281/zenodo.24152490)

### Fixed

- Electron captures. ReacLib fits its electron captures (on ³He and ⁷Be, and
  the pep reaction p + p + e⁻) per unit electron density, ρYe. Version 1.0.0
  left that factor out, so these captures ran too slowly: by a factor of about
  3 in the hydrogen-burning core of a massive star, and about 100 in the Sun.
- Reverse rates. The partition-function factor of a ReacLib reverse rate goes
  on the member that ReacLib derived: the one its `v` flag marks when the fit
  labels carry it, otherwise the photodisintegration of a capture pair,
  otherwise the endothermic member. Version 1.0.0 put the factor on the wrong
  member of 1,862 of the 32,676 pairs in the ReacLib snapshot of 29 March 2022.
- The PDF manual shows its equations again.

### Added

- `evolve_zone(..., method="adaptive_implicit_euler")`, the integrator of
  NucNet Tools. It keeps working where BDF stalls as a fuel runs out in a large
  network.
- `electron_fraction()`. `ReactionNetwork.rates()` accepts `ye`.
  `EvolutionResult` gains `steps` and `rejected_steps`.

### Changed

- Python 3.10 or later is required.

## 1.0.0 (2026-06-19)

DOI [10.5281/zenodo.20756798](https://doi.org/10.5281/zenodo.20756798)

First public release.
