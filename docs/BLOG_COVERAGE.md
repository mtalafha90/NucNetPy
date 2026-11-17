# Blog workflow coverage

This note maps the science workflows of Bradley Meyer's NucNet Tools blog
(<https://sourceforge.net/u/mbradle/blog/>, 94 posts, 2012–2019) onto
NucNetPy. The 24 posts on installation and infrastructure (Cygwin, Boost,
compilers, AWS, SourceForge housekeeping) do not apply to a pure-Python package
and are left out.

## Covered

| Blog workflow | NucNetPy feature |
|---|---|
| Running a first network calculation; changing input conditions; dynamic H burning | `evolve_zone`, `constant_thermo`, CLI `evolve-zone` |
| Analysing a network calculation; printing properties and mass fractions | `analysis`, CLI `print-output`, `largest-x`, `zone-*` |
| Understanding, updating and validating nuclear and reaction data | `io.xml`, `io.jina`, `validate_network`, CLI `validate` |
| Merging nuclear and reaction data | `combine_jina_xml`, CLI `jina-combine` |
| Calculating reaction rates | `RateFit`, `TabularRate`, CLI `rates` |
| Comparing forward and reverse reaction rates | `reverse_rate`, `reverse_reaction` (detailed balance) |
| Computing reaction flows (forward and reverse) | `net_flows`, CLI `net-flows`; forward only: `flows` |
| Constructing the rate equations | `solver.rhs`, `analytic_jacobian`, `jacobian`, `jacobian_sparsity` |
| Finding valid reactions; conservation | `invalid_reactions`, `conserves_a_z`, CLI `conservation` |
| Modifying reaction rates during a calculation; user-defined rate modifiers | `rate_modifiers` |
| Specifying an abundance in a network calculation | `Zone.set_abundance` |
| Selecting an input network; limiting the reaction network | `select_species`, `limit_network` |
| Nuclide XPath / XSLT node-set extraction | `select_species` (selection by Z, A or element) |
| Sorting species in a network calculation | `Network.species_names` (Z, A order) |
| Calculating nuclear statistical equilibrium | `solve_nse`, CLI `nse` |
| Comparing network calculations to equilibrium | `solve_nse` and `regression_summary` |
| Studying the nuclear partition function | `Species.partition`, partition-aware NSE |
| Understanding abundance moments | `abundance_moment` |
| Computing the number of heavy nuclei | `heavy_nuclei_abundance` |
| Computing system timescales | `system_timescales`, CLI `timescales` |
| Computing charge-changing flows | `charge_changing_flows`, CLI `charge-flows` |
| Computing the s-process neutron exposure | `neutron_exposure` |
| Plotting abundances versus nucleon number | `abundances_vs_nucleon_number` |
| Plotting isotopic abundances versus time | `species_history`, `EvolutionResult.mass_fraction_history` |
| Computing separation energies | `separation_energy` |
| Computing and understanding the entropy generation rate (including weak decay) | `entropy_generation_rate`, CLI `entropy-generation` |
| Computing thermodynamic quantities | `thermo`, `electrons` (electron-positron gas: chemical potential, pressure, energy, entropy) |
| Inverting the entropy to find the density or temperature | `entropy_per_nucleon` (ions with masses, spins and partition functions; electrons and positrons; photons; optional Coulomb term), `density_for_entropy`, `t9_for_entropy` |
| Constant-entropy evolution | `constant_entropy_thermo`, which sets T9 from the current composition at each step; the reactions' heat is not fed back |
| Defining your own thermodynamic (trajectory) function | `ThermoFunction` callables, `hydro.Trajectory` |
| Using a thermodynamic trajectory file | `read_trajectory` |
| Running a network calculation with simple hydrodynamics | `exponential_expansion` |
| Including electron screening; supplying your own screening function | `screening`, the solver's `screening=` hook |
| Electron screening in non-degenerate matter | `SkyNetScreening`, which takes the electron and positron term from the electron gas (`ElectronGas.screening_term`); `pair_term` overrides it |
| Adding fission to an r-process calculation | `fission_reaction` |
| Running an r-process calculation | `evolve_zone` with a JINA XML network |
| (n,γ)-(γ,n) equilibrium studies | `reverse_reaction` to build photodisintegration partners; `solve_nse` |
| Creating flow diagrams (colours, scalings, currents) | `graph.reaction_network_dot`, CLI `net-dot` (basic) |
| Creating webnucleo reaction XML from text with Python | `read_reaclib_text`, `io.text`, `write_xml` |
| Using zone properties; global strings | `Zone.properties`, `Zone.optional_properties`, `Network.metadata` |
| Selecting a zone as input for a new calculation | CLI `export-zone-xml` |
| Controlling the frequency of writing output | the `times` grid passed to `evolve_zone` |
| Using wnutils to analyse output | the `analysis` module and `io.hdf5` |
| Weak decays; entropy for weak decay | `weak` tables and `charge_changing_flows` |
| Updating mass excesses in a nuclide XML | edit `Species.mass_excess`, then `write_xml` |
| Counting matrix elements; arrow matrix storage | `jacobian_sparsity` (SciPy sparse matrices replace the arrow solver) |

## Written from the C++ source

These routines were written by working from the named files of the NucNet
Tools SVN trunk (`nucnet-tools-code`, revision 647, from SourceForge):

| Blog workflow | C++ origin | NucNetPy feature |
|---|---|---|
| Studying Coulomb corrections to NSE; including NSE corrections in a network calculation | `user/nse_corr.cpp` (Bravo & García-Senz 1999, Ogata & Ichimaru constants) | `coulomb` module: `gamma_e`, `species_coulomb_chemical_potential` / `_energy` / `_entropy`, `nse_correction`; `solve_nse(..., nse_correction=...)` hook |
| Computing NSE corrections to the entropy | `user/nse_corr.cpp` (`species_coulomb_entropy`) | `coulomb_entropy_per_nucleon` |
| Comparing network calculations to equilibrium; (n,γ)-(γ,n) equilibrium; constraining the equilibrium | libnuceq clusters via `examples/analysis/compare_equil.cpp` | `solve_qse` with `QSECluster` constraints (one multiplier per cluster on top of μ_p and μ_n) |
| Calculating cluster flows | `examples/analysis/compute_Ycdot.cpp` | `cluster_ydot`, `cluster_abundance` |
| Creating and analysing integrated currents diagrams | `user/flow_utilities.cpp` (`update_flow_currents`: current += (f−r)·dt) | `integrated_currents` over an `EvolutionResult` |
| Computing the entropy generation rate (per reaction, with reverse flows) | `user/flow_utilities.cpp` (`compute_entropy_generation_rate`) | `entropy_generation_rate(use_reverse=True)` — Σ(f−r)·ln(f/r), non-negative, vanishes at NSE; `reaction_entropy_changes` |

## Not covered, or left to other tools

| Blog workflow | In NucNetPy |
|---|---|
| Electron and positron chemical-potential terms in the entropy generation (libstatmech) | Not covered: weak reactions carry tabulated rates instead |
| Fission cycling studies with automatic fragment distributions | Partly covered: `fission_reaction` builds the channels; the fragment distributions are user input |
| Downloading webnucleo XML files | Outside the package: download them from the addresses in the blog post; NucNetPy reads them directly |
| XML inclusion (XInclude) and XSLT transforms | Outside the package: apply them with `lxml`, then call `read_xml` |
