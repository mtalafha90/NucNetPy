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

## Notes

These notes are also part of the manual.

| Document | Contents |
|---|---|
| [CONVERSION_MAP.md](CONVERSION_MAP.md) | Where each NucNet Tools / libnucnet concept and source area is found in NucNetPy |
| [BLOG_COVERAGE.md](BLOG_COVERAGE.md) | Each workflow of the NucNet Tools blog and the NucNetPy feature that covers it, including the routines written from the C++ source |
| [EXACT_NUMERICAL_IDENTITY.md](EXACT_NUMERICAL_IDENTITY.md) | How to compare NucNetPy with a NucNet Tools build, and the golden-output tests that make the comparison |
| [JINA_XML_WORKFLOW.md](JINA_XML_WORKFLOW.md) | Working with separate JINA nuclide, reaction and zone XML files |

To check a production JINA database of your own, run
`validation/validate_real_jina_xml.py`; it prints a summary of the database and
writes it to a JSON file. Notebook
[`08_validate_real_jina_files.ipynb`](../notebooks/08_validate_real_jina_files.ipynb)
walks through the same check interactively.

The main [README](../README.md) covers installation, quick start, and the CLI.
Tutorial notebooks live in [`notebooks/`](../notebooks/), and the history of
each release is in [`CHANGELOG.md`](../CHANGELOG.md).
