# Comparing with a NucNet Tools build

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

## The golden-output tests

`tests/test_golden_identity.py` compares NucNetPy with reference values kept in
`tests/golden/`, in this order:

1. the XML round trip: `test_xml_round_trip` and
   `test_round_trip_is_a_fixed_point`, which need no reference values;
2. ReacLib rates on a $T_9$ grid: `test_reaclib_rates_on_t9_grid` against
   `rates_reaclib.json`;
3. screening factors: `test_screening_factors` against `screening.json`;
4. weak-rate table interpolation: `test_weak_rate_table_interpolation`
   against `weak_rates.json`;
5. the one-zone right-hand side at fixed states: `test_ydot_at_fixed_states`
   against `ydot.json`, which also checks that the right-hand side conserves
   nucleon number;
6. a full trajectory on a fixed time grid: `test_full_trajectory` against
   `trajectory.json`, integrated with fixed-step RK4, and
   `test_trajectory_solver_cross_check`, which checks that adaptive BDF agrees
   with it.

As shipped, the reference values are snapshots of NucNetPy's own output, so the
same tests catch any unintended change to the numerics.
`validation/generate_golden.py` regenerates them after an intended change.

## How to compare with a C++ build you run yourself

1. Run the original NucNet Tools on the inputs in `tests/golden/`
   (`golden_network.xml` and `golden_weak_rates.txt`), or replace those inputs
   with your own network and run both codes on them.
2. Write the C++ outputs into the `data` blocks of the JSON files in
   `tests/golden/`, and set each file's `source` to a label for your build.
3. Set each file's `rtol` and `atol` to the agreement you need. The tests read
   the tolerances from the files, so no test code has to change.
4. Run `pytest tests/test_golden_identity.py`. A failure names the quantity,
   the value from each code and the tolerance.

## Adding a routine written from a C++ example

When a routine is written from a C++ example, add a golden-output test for it:

1. Run the original C++ program on a fixed input file and save its output.
2. Run the Python routine on the same input.
3. Compare with zero tolerance first. If the platforms differ in the last
   bits, record the tolerance you need and why.

Only routines that pass such a test should be described as exact replacements.
