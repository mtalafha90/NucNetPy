"""Sphinx configuration for the NucNetPy user manual.

Build with:

    python -m pip install -e ".[docs]"
    sphinx-build -b html docs/manual docs/manual/_build/html
    sphinx-build -b singlehtml docs/manual docs/manual/_build/singlehtml

The single-file build is the one to place in a submission archive.
"""
from __future__ import annotations

import importlib.metadata

project = "NucNetPy"
copyright = "2026, M. H. Talafha and N. M. Ershaidat"
author = "M. H. Talafha and N. M. Ershaidat"
try:
    release = importlib.metadata.version("nucnetpy")
except importlib.metadata.PackageNotFoundError:  # building without an install
    release = "1.0.0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",     # NumPy-style docstrings
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",             # so the existing .md notes can be included
]

templates_path = []
exclude_patterns = ["_build"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"

# Order members as they appear in the source rather than alphabetically: the
# modules are written so that reading top to bottom makes sense.
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}

myst_enable_extensions = ["dollarmath", "colon_fence"]

html_theme = "furo"
html_title = f"NucNetPy {release}"
html_static_path = []
