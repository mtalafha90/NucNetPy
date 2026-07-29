# NucNetPy documentation

## User manual

The user manual is the place to start. It covers installation, the object
model, the common workflows, the pitfalls that produce plausible but wrong
answers, how to verify a calculation, and a full API and command-line
reference generated from the source.

| Form | Location |
|---|---|
| PDF | [`manual/nucnetpy-manual.pdf`](manual/nucnetpy-manual.pdf) |
| Sources | [`manual/`](manual/) |

Rebuild it after changing the code or the docstrings:

```bash
python -m pip install -e ".[docs]"
make -C docs/manual all      # html, singlehtml and pdf
```

The PDF is produced by rinohtype and needs no LaTeX installation.

## Developer notes


| Document | Contents |
|---|---|
| [BLOG_COVERAGE.md](BLOG_COVERAGE.md) | Map of every NucNet Tools blog workflow to its nucnetpy feature, including the items ported directly from the C++ source (r647) and what remains open |
| [CONVERSION_MAP.md](CONVERSION_MAP.md) | Original NucNet Tools / libnucnet concept → nucnetpy API mapping |
| [JINA_XML_WORKFLOW.md](JINA_XML_WORKFLOW.md) | Working with separate JINA nuclide/reaction/zone XML files |
| [PURE_PYTHON_PORT_STATUS.md](PURE_PYTHON_PORT_STATUS.md) | Port status per C++ area and the golden-output numerical-identity test framework |
| [EXACT_NUMERICAL_IDENTITY.md](EXACT_NUMERICAL_IDENTITY.md) | Notes on achieving strict numerical identity with a specific C++ build |
| [REAL_JINA_FILES_UPDATE.md](REAL_JINA_FILES_UPDATE.md), [REAL_XML_VALIDATION_REPORT.md](REAL_XML_VALIDATION_REPORT.md) | Historical validation reports against real JINA XML databases |

The main [README](../README.md) covers installation, quick start, and the CLI.
Tutorial notebooks live in [`notebooks/`](../notebooks/).
