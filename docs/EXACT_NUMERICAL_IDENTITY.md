# Exact numerical identity policy

NucNetPy is pure Python. It contains no C++ code, does not call the original
NucNet Tools programs, and has no "exact backend". It therefore cannot promise
bit-for-bit agreement with a particular C++ build on its own.

## Why identity cannot simply be assumed

Two codes that implement the same equations can still disagree in the last
digits, or by more, because results depend on:

- the nuclear data: mass excesses, partition-function tables and their
  interpolation, and the rate fits;
- modelling choices: the screening prescription, how reverse rates are
  formed, and which species enter an equilibrium solve;
- numerical choices: the integrator, its tolerances and time grid, the order
  in which reactions are summed, and the linear-algebra library;
- the platform: compiler, optimisation flags and floating-point rounding.

Agreement has to be shown by comparison, not asserted.

## How to compare with a C++ build you run yourself

The repository carries the machinery for this; the port status notes
(`docs/PURE_PYTHON_PORT_STATUS.md`) describe it in detail.

1. Run the original NucNet Tools on the inputs in `tests/golden/` (or replace
   those inputs with your own network and run both codes on them).
2. Write the C++ outputs into the `data` blocks of the JSON files in
   `tests/golden/`, and set each file's `source` to a label for your build.
3. Set each file's `rtol` and `atol` to the agreement you need. The tests read
   the tolerances from the files, so no test code has to change.
4. Run `pytest tests/test_golden_identity.py`. A failure names the quantity,
   the value from each code and the tolerance.

Out of the box those files hold snapshots of NucNetPy's own output, so the
same tests also catch any unintended change to the Python numerics.
`validation/generate_golden.py` regenerates the snapshots after an intended
change.

## Development rule

When a C++ example is converted to Python, add a golden-output test for it:

1. Run the original C++ program on a fixed input file and save its output.
2. Run the Python equivalent on the same input.
3. Compare with zero tolerance first. If the platforms differ in the last
   bits, record the tolerance you need and why.

Only modules that pass such tests should be described as exact replacements.
