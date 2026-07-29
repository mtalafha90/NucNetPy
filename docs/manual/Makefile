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

.PHONY: help html singlehtml pdf clean all

help:
	@echo "targets: html  singlehtml  pdf  all  clean"

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(SOURCEDIR) $(BUILDDIR)/singlehtml

pdf:
	$(SPHINXBUILD) -b rinoh $(SOURCEDIR) $(BUILDDIR)/pdf

all: html singlehtml pdf

clean:
	rm -rf $(BUILDDIR)
