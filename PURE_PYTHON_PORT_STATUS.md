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

Recommended order:

1. XML round-trip tests for your real networks.
2. Rate-by-rate ReacLib comparisons on a T9 grid.
3. Screening-factor comparisons for your chosen `user/screen.cpp` model.
4. Weak-rate table comparisons.
5. One-zone RHS `ydot` comparisons at fixed state.
6. Full trajectory comparisons with identical time grids and tolerances.
