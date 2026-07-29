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

## 4. Credentials

PyPI and TestPyPI are separate services with separate accounts; a login on one
does not work on the other. Both authenticate uploads with an API token rather
than a password.

| | TestPyPI | PyPI |
|---|---|---|
| Register | <https://test.pypi.org/account/register/> | <https://pypi.org/account/register/> |
| Create a token | <https://test.pypi.org/manage/account/token/> | <https://pypi.org/manage/account/token/> |

Two-factor authentication must be enabled before a token can be created. The
token is displayed **once**; copy it then.

Scope the first token to the whole account, because a project-scoped token can
only be created after the project exists. Once `nucnetpy` is on the index,
delete it and create a project-scoped token in its place.

When `twine` prompts, the username is the literal string `__token__` and the
password is the entire token including its `pypi-` prefix. To avoid retyping,
put it in `~/.pypirc`:

```ini
[distutils]
index-servers = pypi testpypi

[pypi]
username = __token__
password = pypi-<your-PyPI-token>

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-<your-TestPyPI-token>
```

```bash
chmod 600 ~/.pypirc
```

That file holds credentials in clear text. Never commit it, and never put a
token in anything git tracks.

### Trusted Publishing

PyPI can instead accept a short-lived OIDC identity from GitHub Actions, so no
long-lived token exists to leak or rotate. Register the repository once under
*Publishing* in the PyPI project settings, and a tagged release publishes
itself. For a package accompanying a paper, where releases are rare and the
archived artefact should match the tag exactly, this is the more robust
arrangement.

## 5. TestPyPI first

Upload to the test index and prove an install from it before touching PyPI.

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

## 6. PyPI

```bash
python -m twine upload dist/*
#   username: __token__
#   password: pypi-...        (the PyPI token)
```

Use a project-scoped API token once the project exists. Never put a token in a
file that git tracks; `~/.pypirc` should be mode 600.

## 7. After

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
