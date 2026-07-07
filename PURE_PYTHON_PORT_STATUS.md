# Pure Python port status

The C++ source contains helper files under `nnt/` and `user/`. This alpha adds Python equivalents for the practical computational layers rather than a backend wrapper.

| Original area | Python module | Status |
|---|---|---|
| `nnt/math.*` | `nucnetpy.mathutils` | Ported: linear, bilinear, 2-D interpolation behavior |
| wrappers/containers | `core.py`, `species.py` | Ported core containers |
| network/reaction utilities | `reactions.py`, `validation.py` | Ported main rate/flow/conservation APIs |
| `screen.*` | `screening.py` | Pure-Python weak and intermediate screening; exact formula variants still need golden validation |
| `weak_utilities.*`, `two_d_weak_rates.*` | `weak.py` | Ported 2-D weak-rate table interpolation and weak source terms |
| `thermo.*` | `thermo.py`, `nse.py` | Partial: common thermodynamic helpers and NSE; detailed libstatmech parity requires validation |
| `matrix_solver.*`, `ilu_solvers.*` | `matrix_solver.py`, `solver.py` | Ported to SciPy sparse solvers and stiff integrators |
| `hydro.*`, `hydro_helper.*` | `hydro.py` | Ported trajectory interpolation and exponential expansion |
| `network_limiter.*` | `network_limiter.py` | Ported species/reaction subnetwork limiting |
| `nuclear_decay_utilities.*` | `decay.py` | Ported half-life/decay record tools |
| `neutrino_rate_functions.*` | `neutrino.py` | Ported neutrino luminosity/cross-section helpers |
| `hdf5_routines.*` | `io/hdf5.py` | Ported compact HDF5 save/load |
| `graph.h` | `graph.py` | Ported DOT export |

## What remains for strict numerical identity

No C++ is used in this package. To make the Python code exactly match a chosen original C++ build, add golden-output tests generated once from the original tools, then tune the pure-Python formulas and tolerances until those tests pass.

The framework for this now exists in-repo. `tests/test_golden_identity.py` implements the full recommended order below against golden files in `tests/golden/`, and `validation/generate_golden.py` regenerates those files:

1. XML round-trip tests for your real networks — `test_xml_round_trip`, `test_round_trip_is_a_fixed_point` (identity tests, no golden data needed).
2. Rate-by-rate ReacLib comparisons on a T9 grid — `test_reaclib_rates_on_t9_grid` vs `rates_reaclib.json`.
3. Screening-factor comparisons for your chosen `user/screen.cpp` model — `test_screening_factors` vs `screening.json`.
4. Weak-rate table comparisons — `test_weak_rate_table_interpolation` vs `weak_rates.json`.
5. One-zone RHS `ydot` comparisons at fixed state — `test_ydot_at_fixed_states` vs `ydot.json` (also asserts nucleon-number conservation of the RHS).
6. Full trajectory comparisons with identical time grids and tolerances — `test_full_trajectory` vs `trajectory.json` (deterministic fixed-step rk4), plus a solver cross-check that adaptive BDF agrees with the rk4 path.

Out of the box the goldens are self-consistent snapshots of the current Python numerics, so any future change to a formula fails the suite immediately. To validate against a specific C++ NucNet Tools build:

1. Replace the fixture inputs (`tests/golden/golden_network.xml`, `golden_weak_rates.txt`) with your production data, or keep them and reproduce the same inputs on the C++ side.
2. Run the same grids through the original tools and write their outputs into the `data` blocks of the JSON files (set `source` to your build id).
3. Loosen each file's `rtol`/`atol` to the agreement you require — the tolerances are read from the golden files, not hard-coded in the tests.
