# Installation

## Requirements

Python 3.9 or later, with NumPy 1.22+ and SciPy 1.9+. Nothing else is required
to use the package.

## From a checkout

```bash
git clone https://github.com/mtalafha90/NucNetPy.git
cd NucNetPy
python -m pip install -e .
pytest -q
```

The test suite runs in about a second and should report all tests passing
before you rely on any result.

## Optional extras

| Extra | Adds | For |
|---|---|---|
| `dev` | `pytest` | running the test suite |
| `hdf5` | `h5py` | HDF5 persistence |
| `plot` | `matplotlib`, `networkx` | plotting and graph export |
| `notebook` | `jupyter`, `pandas`, `matplotlib` | the tutorial notebooks |
| `docs` | `sphinx`, `furo`, `myst-parser` | building this manual |

```bash
python -m pip install -e ".[dev,hdf5,plot,notebook]"
```

`pandas` is used only by the notebooks to display tables; the package itself
never imports it.

## Verifying the installation

```python
import nucnetpy as nn
print(nn.__version__)
```

The command-line interface is installed alongside:

```bash
nucnetpy --help
```

## Building this manual

```bash
python -m pip install -e ".[docs]"
sphinx-build -b html docs/manual docs/manual/_build/html
```

For a single self-contained file, suitable for a submission archive:

```bash
sphinx-build -b singlehtml docs/manual docs/manual/_build/singlehtml
```
