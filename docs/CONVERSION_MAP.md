# From NucNet Tools to NucNetPy

NucNetPy is written in Python from the start: it does not call or wrap the
C++ code of NucNet Tools. This note shows where each part of NucNet Tools and
libnucnet is found in NucNetPy, for readers who know the original.

## Concepts

| NucNet Tools / libnucnet | NucNetPy |
|---|---|
| `Libnucnet__Nuc` / species handling | `nucnetpy.species.Species`, `Network.species` |
| `Libnucnet__Reac` / reactions | `nucnetpy.reactions.Reaction`, `ReactionNetwork`, `RateFit`, `TabularRate` |
| ReacLib seven-parameter rates | `RateFit.rate(t9)` |
| Reaction duplicate/invalid utilities | `ReactionNetwork.remove_duplicates`, `filter_valid`, CLI `remove-duplicates`, `remove-invalid` |
| Zone XML and abundance output | `Zone`, `Network.zones`, `nucnetpy.io.xml` |
| Analysis examples (`print_output`, largest X, flows, ydot) | `nucnetpy.analysis` and CLI commands |
| Graph utilities | `nucnetpy.graph.reaction_network_dot`, CLI `net-dot` |
| Screening user callbacks | `nucnetpy.screening` (`SkyNetScreening`, weak and intermediate screening) and the solver's `screening=` hook |
| Weak-rate utilities | `nucnetpy.weak.WeakRateTable`, `read_weak_table`, and the solver's weak-rate hook |
| `libnuceq` / NSE workflows | `nucnetpy.nse.solve_nse`, `equilibrium_ratio`; `nucnetpy.qse.solve_qse` for constrained clusters |
| Matrix/Jacobian solver utilities | `nucnetpy.solver.analytic_jacobian`, `jacobian`, `jacobian_sparsity`; SciPy BDF/Radau/LSODA and the adaptive implicit Euler scheme |
| Multi-zone evolution | `evolve_network_zones` |
| HDF5 helper workflow | `nucnetpy.io.hdf5` |

## Source areas

| NucNet Tools area | NucNetPy module | What it provides |
|---|---|---|
| `nnt/math.*` | `nucnetpy.mathutils` | Linear, bilinear and two-dimensional interpolation |
| wrappers and containers | `core.py`, `species.py` | `Network`, `Zone` and `Species` |
| network and reaction utilities | `reactions.py`, `validation.py` | Rates, flows and conservation checks |
| `screen.*` | `screening.py` | Weak, intermediate and SkyNet-style screening |
| `weak_utilities.*`, `two_d_weak_rates.*` | `weak.py` | Two-dimensional weak-rate tables and their source terms |
| `thermo.*` | `thermo.py`, `electrons.py`, `nse.py` | Entropy of ions, electrons and positrons, and photons; density and temperature from the entropy; NSE. The electron gas is checked against closed-form limits and thermodynamic identities |
| `matrix_solver.*`, `ilu_solvers.*` | `matrix_solver.py`, `solver.py` | SciPy sparse solvers and stiff integrators, and the adaptive implicit Euler scheme |
| `hydro.*`, `hydro_helper.*` | `hydro.py` | Trajectory interpolation and exponential expansion |
| `network_limiter.*` | `network_limiter.py` | A network limited to chosen species and the reactions among them |
| `nuclear_decay_utilities.*` | `decay.py` | Half-life and decay records, and fission channels |
| `neutrino_rate_functions.*` | `neutrino.py` | Neutrino luminosity and cross-section helpers |
| `hdf5_routines.*` | `io/hdf5.py` | Compact HDF5 output and input |
| `graph.h` | `graph.py` | Graphviz DOT export |

## Matching a particular C++ run

To reproduce a particular C++ run, both codes need the same inputs and the same
choices. The points to match are:

1. the XML schema variants of your input files;
2. the mass excesses and partition-function data;
3. the screening callback chosen in the C++ build;
4. the source and interpolation convention of the weak-rate tables;
5. the solver, its tolerances and time steps, and any rate modifiers.

`docs/EXACT_NUMERICAL_IDENTITY.md` describes how to compare the two codes.
