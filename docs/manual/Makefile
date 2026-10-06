# Build the NucNetPy user manual.
#
#   python -m pip install -e ".[docs]"
#   make -C docs/manual html      # browsable HTML
#   make -C docs/manual pdf       # single PDF, no LaTeX required
#   make -C docs/manual singlehtml
#
# The PDF is the artefact to place in a submission archive.

SPHINXBUILD ?= sphinx-build
SOURCEDIR    = .
BUILDDIR     = _build

PYTHON      ?= python

.PHONY: help cli html singlehtml pdf clean all

help:
	@echo "targets: cli  html  singlehtml  pdf  all  clean"

# The command-line reference is the parser's own help text; rebuild it first
# so the manual cannot drift from the code.
cli:
	$(PYTHON) generate_cli_reference.py

html: cli
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html

singlehtml: cli
	$(SPHINXBUILD) -b singlehtml $(SOURCEDIR) $(BUILDDIR)/singlehtml

pdf: cli
	$(SPHINXBUILD) -b rinoh $(SOURCEDIR) $(BUILDDIR)/pdf

all: html singlehtml pdf

clean:
	rm -rf $(BUILDDIR)
