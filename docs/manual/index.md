# NucNetPy user manual

NucNetPy is a pure-Python nuclear reaction-network package. It reads
JINA/libnucnet XML data, evaluates reaction rates, builds and integrates the
abundance equations with stiff solvers, solves nuclear statistical and
quasi-statistical equilibrium, and provides the analysis quantities used in the
NucNet Tools workflows.

It is a clean reimplementation: it does not call or wrap the original C/C++
code, and NumPy and SciPy are its only mandatory dependencies.

```{toctree}
:maxdepth: 2
:caption: User guide

guide/installation
guide/concepts
guide/workflows
guide/pitfalls
guide/verification
```

```{toctree}
:maxdepth: 2
:caption: Reference

reference/api
reference/cli
```

```{toctree}
:maxdepth: 1
:caption: Developer notes

notes/conversion
notes/blog_coverage
notes/port_status
notes/numerical_identity
notes/jina_workflow
```

## Where to start

If you have never used the package, read {doc}`guide/installation` and
{doc}`guide/concepts`, then work through the tutorial notebooks in
`notebooks/`, which are the fastest way in.

If you are migrating an existing NucNet Tools workflow, {doc}`notes/conversion`
maps the original concepts onto this API and {doc}`notes/blog_coverage` maps
each blog workflow onto the feature that replaces it.

If you are about to run a production calculation, read {doc}`guide/pitfalls`
first. It describes two behaviours that produce plausible but wrong answers
without raising an error.

## Citing

Please cite the archived release
([10.5281/zenodo.20756798](https://doi.org/10.5281/zenodo.20756798)) and
acknowledge the NucNet Tools and libnucnet ecosystem, whose data model and
workflows this package reimplements, together with the JINA reaction-rate
database that supplies the nuclear data.
