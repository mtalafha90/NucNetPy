# Reference case: silicon burning

This directory is the sample input and output that ships with the program. It
makes the principal scientific result of the accompanying article reproducible
**without the third-party JINA/libnucnet database**, which cannot be
redistributed.

```bash
python validation/reproduce_si_burning.py
```

Exit status 0 means every archived value was reproduced within the tolerances
recorded beside them. Nothing but the installed package is needed.

## What is here

| File | Contents |
|---|---|
| `si_burning_network.xml` | The extracted network: 15 nuclides, 32 reactions, and the initial zone. Complete input — rates, mass excesses, spins and partition functions included |
| `si_burning_expected.json` | Conditions, expected final abundances, the NSE solution, and the comparison tolerances |
| `si_burning_provenance.json` | Checksums, sizes and data-release labels of the database this was extracted from |

The network is 58 kB, so it costs nothing to archive, and it is the whole
input: the 15-nuclide alpha chain from `he4` to `ni56` with free nucleons and
the photon, which must be kept or every photodisintegration record is dropped
with it.

## The conditions

Pure <sup>28</sup>Si, `Y = 1/28`, giving `Ye = 0.5` exactly, burned at a
constant `T9 = 5` and `rho = 1e8 g/cm^3` for 10 s on 120 output points, with
BDF at `rtol = 1e-8` and `atol = 1e-14`. The stationary composition is then
compared with `solve_nse` at the same conditions.

## Provenance of the source database

The database itself is not redistributed. `si_burning_provenance.json`
identifies it: SHA-256 and byte size of each file, the libnucnet schema version
(2019-01-15), and a histogram of the `<source>` labels that identifies the data
release — `frdm`, `ame11` and `reac1` for the masses, `rath`, `ths8`, `mo03`
and `wc12` among the rates.

It also records the species counts, which are easy to misread. The nuclide file
holds 7852 entries; the parsed network holds 8131 species, because 279 are
placeholders synthesised for species that reaction records name and the nuclide
file does not describe.

## Regenerating this directory

Only necessary if the reference case itself should change. It requires the
database and overwrites the expectations, so the diff should be read carefully:

```bash
python validation/extract_reference_case.py \
    --nuclides nuclides.xml --reactions reaction_data.xml
```

`tests/test_reference_case.py` checks the archive on every test run, including
that the network still matches the checksum its expectations were built from.
A change in the code that moves the reference result fails there first.
