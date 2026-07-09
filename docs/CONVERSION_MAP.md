# Conversion map from original NucNet Tools concepts to `nucnetpy`

| Original area | Python replacement |
|---|---|
| `Libnucnet__Nuc` / species handling | `nucnetpy.species.Species`, `Network.species` |
| `Libnucnet__Reac` / reactions | `nucnetpy.reactions.Reaction`, `ReactionNetwork`, `RateFit`, `TabularRate` |
| ReacLib seven-parameter rates | `RateFit.rate(t9)` |
| Reaction duplicate/invalid utilities | `ReactionNetwork.remove_duplicates`, `filter_valid`, CLI `remove-duplicates`, `remove-invalid` |
| Zone XML and abundance output | `Zone`, `Network.zones`, `nucnetpy.io.xml` |
| Analysis examples (`print_output`, largest X, flows, ydot) | `nucnetpy.analysis` and CLI commands |
| Graph utilities | `nucnetpy.graph.reaction_network_dot`, CLI `net-dot` |
| Screening user callbacks | `nucnetpy.screening` with weak/intermediate models and solver callback hook |
| Weak-rate utilities | `nucnetpy.weak.WeakRateTable`, `read_weak_table`, solver weak-rate hook |
| `libnuceq`/NSE workflows | `nucnetpy.nse.solve_nse`, `equilibrium_ratio` |
| Matrix/Jacobian solver utilities | `nucnetpy.solver.jacobian`, `jacobian_sparsity`, SciPy BDF/Radau/LSODA |
| Multi-zone evolution | `evolve_network_zones` |
| HDF5 helper workflow | `nucnetpy.io.hdf5` |

## What is implemented in alpha 2

- Flexible XML parsing for common legacy structures and compact Python-written XML.
- Partition functions and optional properties are preserved.
- ReacLib, tabular, and constant rates are supported.
- NSE solves are available for network-level equilibrium initialization.
- Screening and weak-rate hooks are now solver-compatible.
- A sparse numerical Jacobian and sparsity pattern are available for stiff integration.
- Validation/regression helpers are included.

## What still needs project-specific validation

A full replacement can now be used and extended, but exact reproduction of a historical C++ run requires comparing against your original input/output files.  The main tunable points are:

1. exact XML schema variants used by your input files,
2. original mass excess and partition-function data,
3. exact screening callback selected in the original C++ build,
4. weak-rate table source and interpolation convention,
5. solver tolerance, timestep, and rate-modifier settings.
