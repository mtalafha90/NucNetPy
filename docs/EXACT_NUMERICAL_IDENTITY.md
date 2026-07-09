# Exact numerical identity policy

A pure Python rewrite cannot honestly promise exact numerical identity with the
original C++ NucNet Tools by itself.  Identity depends on the original C/C++
source, libnucnet/libnuceq/statmech versions, GSL/libxml2 behavior, compiler,
optimization flags, reaction ordering, matrix solver details, and floating-point
rounding.

Therefore this package uses two backends:

1. **Exact backend**: calls the original C++ NucNet Tools executables or compiled
   wrappers.  This is the backend to use when you need exact legacy identity.
2. **Python backend**: a Python implementation/scaffold for analysis, testing,
   and future line-by-line porting.  This is useful, but it is not the source of
   truth for exact identity until validated against golden C++ outputs.

## Build the original backend

The original C++ source tree is bundled under:

```text
external/original_cpp
```

On your machine, install the original dependencies, usually including:

```bash
sudo apt install build-essential g++ make gsl-bin libgsl-dev libxml2-dev xsltproc wget
```

Then try:

```bash
nucnetpy-exact source
nucnetpy-exact make --jobs 4
```

The legacy Makefiles may download specific vendor packages such as libnucnet,
libnuceq, wn_matrix, and statmech.  For strict reproducibility, keep the same
versions and compile flags as your original C++ installation.

## Run an original executable through Python

Example:

```bash
nucnetpy-exact run print_output my_output.xml
```

or in Python:

```python
from nucnetpy.exact import CppBackend

cpp = CppBackend(bin_dir="/path/to/original/build/bin")
result = cpp.run("print_output", ["my_output.xml"])
print(result.stdout)
```

## Compare outputs

Byte-for-byte text comparison:

```bash
nucnetpy-exact compare cpp_output.txt python_output.txt
```

Floating-value comparison:

```bash
nucnetpy-exact compare cpp_output.txt python_output.txt --float --rtol 1e-14 --atol 0
```

For exact identity, use `--rtol 0 --atol 0`.  In practice, different compilers or
CPUs may produce last-bit differences even for the same C++ code.

## Development rule

For every C++ example converted to Python, create a golden-output test:

1. Run the original C++ executable on the same input XML.
2. Save stdout as the golden file.
3. Run the Python equivalent.
4. Compare with zero tolerance first; if platform differences appear, record the
   exact tolerance and reason.

Only modules passing these golden tests should be advertised as exact replacements.
