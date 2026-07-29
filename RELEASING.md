# Releasing NucNetPy

A PyPI upload is permanent: a version number can never be reused, even after
the file is deleted. Everything below is designed so that mistakes are caught
before that point.

## 1. Before building

```bash
pytest -q                                   # must be all green
make -C docs/manual all                     # rebuild the manual if code changed
cp docs/manual/_build/pdf/nucnetpy.pdf docs/manual/nucnetpy-manual.pdf
```

Check that `version` in `pyproject.toml`, `CITATION.cff`, and the Zenodo record
agree.

## 2. Build

```bash
python -m pip install --upgrade build twine
rm -rf dist build src/*.egg-info
python -m build
```

This produces `dist/nucnetpy-<version>-py3-none-any.whl` and
`dist/nucnetpy-<version>.tar.gz`.

## 3. Validate

```bash
python -m twine check dist/*
```

Then verify a clean install actually works, which `twine check` does not do:

```bash
python -m venv /tmp/relcheck
/tmp/relcheck/bin/pip install dist/nucnetpy-*.whl
/tmp/relcheck/bin/python -c "import nucnetpy; print(nucnetpy.__version__)"
/tmp/relcheck/bin/nucnetpy --help
```

And confirm the sdist is self-contained, since that is what an archive or a
reviewer receives:

```bash
tar -xzf dist/nucnetpy-*.tar.gz -C /tmp
python -m venv /tmp/sdistcheck
/tmp/sdistcheck/bin/pip install /tmp/nucnetpy-<version> pytest
cd /tmp/nucnetpy-<version> && /tmp/sdistcheck/bin/python -m pytest -q
```

The sdist must carry the tests, `validation/`, `notebooks/`, `examples/`, the
manual, `LICENSE` and `CITATION.cff`. `MANIFEST.in` controls this.

## 4. TestPyPI first

TestPyPI is a separate index with separate accounts and tokens. Register at
<https://test.pypi.org/account/register/> and create an API token.

```bash
python -m twine upload --repository testpypi dist/*
#   username: __token__
#   password: pypi-...        (the TestPyPI token)
```

Install from it in a fresh environment. TestPyPI does not mirror the real
index, so dependencies must still come from PyPI:

```bash
python -m venv /tmp/testpypi
/tmp/testpypi/bin/pip install \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    nucnetpy
/tmp/testpypi/bin/python -c "import nucnetpy; print(nucnetpy.__version__)"
```

Look at the rendered project page on TestPyPI: the README, the author names,
the licence, the classifiers and the links should all be right. This is the
last cheap opportunity to fix them.

## 5. PyPI

```bash
python -m twine upload dist/*
#   username: __token__
#   password: pypi-...        (the PyPI token)
```

Use a project-scoped API token once the project exists. Never put a token in a
file that git tracks; `~/.pypirc` should be mode 600.

## 6. After

```bash
git tag -a v<version> -m "NucNetPy <version>"
git push origin v<version>
```

Then update the Zenodo deposit so the archived record matches the released
source, and confirm the DOI in `CITATION.cff`, `pyproject.toml` and the
manuscript still resolve to it.

## Version numbers

`version` lives in `pyproject.toml` and is the single source of truth;
`nucnetpy.__version__` is set from the same string in `src/nucnetpy/__init__.py`
and must be updated with it. A released version is never rebuilt or reuploaded
— fix forward with a new patch version.
