# Changelog

The changes in each release of NucNetPy. Every release is archived on Zenodo
with its own DOI; [10.5281/zenodo.20756797](https://doi.org/10.5281/zenodo.20756797)
covers all versions and resolves to the latest.

## Unreleased

### Changed

- Every command now has a one-line description in `nucnetpy --help`.
- The documentation describes the package as it is; the history of each
  release is kept here.

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
