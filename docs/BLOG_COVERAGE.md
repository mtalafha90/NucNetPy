# Blog workflow coverage

Mapping of the science workflows from Bradley Meyer's NucNet Tools blog
(<https://sourceforge.net/u/mbradle/blog/>, 94 posts, 2012–2019) to nucnetpy
features. Installation/infrastructure posts (Cygwin, Boost, compilers, AWS,
SourceForge housekeeping — 24 posts) do not apply to a pure-Python package and
are omitted.

## Covered

| Blog workflow | nucnetpy feature |
|---|---|
| Running a first network calculation; changing input conditions; dynamic H burning | `evolve_zone`, `constant_thermo`, CLI `evolve-zone` |
| Analyzing a network calculation; printing properties and mass fractions | `analysis`, CLI `print-output`, `largest-x`, `zone-*` |
| Understanding/updating/validating nuclear and reaction data | `io.xml`, `io.jina`, `validate_network`, CLI `validate` |
| Merging nuclear and reaction data | `combine_jina_xml`, CLI `jina-combine` |
| Calculating reaction rates | `RateFit`, `TabularRate`, CLI `rates` |
| **Comparing forward and reverse reaction rates** | **`reverse_rate`, `reverse_reaction` (detailed balance)** |
| **Computing reaction flows (forward + reverse)** | **`net_flows`, CLI `net-flows`**; forward-only: `flows` |
| Constructing the rate equations | `solver.rhs`, `jacobian`, `jacobian_sparsity` |
| Finding valid reactions; conservation | `invalid_reactions`, `conserves_a_z`, CLI `conservation` |
| Modifying reaction rates during a calculation; user-defined rate modifiers | `rate_modifiers` |
| Specifying an abundance in a network calculation | `Zone.set_abundance` |
| Selecting an input network; limiting the reaction network | **`select_species` (new)**, `limit_network` |
| Nuclide XPath / XSLT node-set extraction | **`select_species` (Z/A/element selection)** |
| Sorting species in a network calculation | `Network.species_names` (Z, A order) |
| Calculating nuclear statistical equilibrium | `solve_nse`, CLI `nse` |
| Comparing network calculations to equilibrium | `solve_nse` + `regression_summary` |
| Studying the nuclear partition function | `Species.partition`, partition-aware NSE |
| Understanding abundance moments | `abundance_moment` |
| **Computing the number of heavy nuclei** | **`heavy_nuclei_abundance` (new)** |
| **Computing system timescales** | **`system_timescales` (new), CLI `timescales`** |
| **Computing charge-changing flows** | **`charge_changing_flows` (new), CLI `charge-flows`** |
| **Computing the s-process neutron exposure** | **`neutron_exposure` (new)** |
| Plotting abundances versus nucleon number | `abundances_vs_nucleon_number` |
| Plotting isotopic abundances versus time | `species_history`, `EvolutionResult.mass_fraction_history` |
| **Computing separation energies** | **`separation_energy` (proton branch fixed)** |
| **Computing/understanding the entropy generation rate (incl. weak decay)** | **`entropy_generation_rate` (new), CLI `entropy-generation`** |
| Inverting the entropy to find the density | `thermo.density_from_entropy`, `temperature_from_entropy` |
| Computing thermodynamic quantities | `thermo` |
| Defining your own thermodynamic (trajectory) function | `ThermoFunction` callables, `hydro.Trajectory` |
| Using a thermodynamic trajectory file | `read_trajectory` |
| Running a network calculation with simple hydrodynamics | `exponential_expansion` |
| Including electron screening; supplying your own screening function | `screening`, solver `screening=` hook |
| **Adding fission to an r-process calculation** | **`fission_reaction` (new)** |
| Running an r-process calculation | `evolve_zone` + JINA XML network |
| (n,γ)-(γ,n) equilibrium studies | **`reverse_reaction`** to build photodisintegration partners; `solve_nse` |
| Creating flow diagrams (colors, scalings, currents) | `graph.reaction_network_dot`, CLI `net-dot` (basic) |
| Creating webnucleo reaction XML from text with Python | `read_reaclib_text`, `io.text`, `write_xml` |
| Using zone properties; global strings | `Zone.properties`, `Network.metadata` |
| Selecting a zone as input for a new calculation | CLI `export-zone-xml` |
| Controlling the frequency of writing output | `times` grid passed to `evolve_zone` |
| Using wnutils to analyze output | `analysis` module + `io.hdf5` |
| Weak decays / entropy for weak decay | `weak` tables + `charge_changing_flows` |
| Updating mass excesses in a nuclide XML | edit `Species.mass_excess`, re-`write_xml` |
| Counting matrix elements; arrow matrix storage | `jacobian_sparsity` (SciPy sparse replaces arrow solver) |

## Ported from the C++ source (r647)

These items were verified against the actual NucNet Tools SVN trunk
(`nucnet-tools-code` r647, downloaded from SourceForge) and ported from the
named C++ files:

| Blog workflow | C++ origin | nucnetpy feature |
|---|---|---|
| Studying Coulomb corrections to NSE; including NSE corrections in a network calculation | `user/nse_corr.cpp` (Bravo & García-Senz 1999, Ogata & Ichimaru constants) | `coulomb` module: `gamma_e`, `species_coulomb_chemical_potential` / `_energy` / `_entropy`, `coulomb_nse_correction`; `solve_nse(..., nse_correction=...)` hook |
| Computing NSE corrections to the entropy | `user/nse_corr.cpp` (`species_coulomb_entropy`) | `coulomb_entropy_per_nucleon` |
| Comparing network calculations to equilibrium; (n,γ)-(γ,n) equilibrium; constraining the equilibrium | libnuceq clusters via `examples/analysis/compare_equil.cpp` | `solve_qse` with `QSECluster` constraints (one multiplier per cluster on top of μ_p/μ_n) |
| Calculating cluster flows | `examples/analysis/compute_Ycdot.cpp` | `cluster_ydot`, `cluster_abundance` |
| Creating/analyzing integrated currents diagrams | `user/flow_utilities.cpp` (`update_flow_currents`: current += (f−r)·dt) | `integrated_currents` over an `EvolutionResult` |
| Computing the entropy generation rate (per reaction, with reverse flows) | `user/flow_utilities.cpp` (`compute_entropy_generation_rate`) | `entropy_generation_rate(use_reverse=True)` — Σ(f−r)·ln(f/r), non-negative, vanishes at NSE; `reaction_entropy_changes` |

## Not ported (by design or still open)

| Blog workflow | Status |
|---|---|
| Constant-entropy evolution with self-consistent entropy inverter | partial — `thermo` inverters exist; no closed-loop driver |
| Electron/positron chemical-potential terms in entropy generation (libstatmech) | open — weak reactions carry tabulated rates instead |
| Fission *cycling* studies (automatic fragment distributions) | partial — `fission_reaction` builds channels; distributions are user input |
| Downloading webnucleo XML files | use the URLs in the blog post; nucnetpy reads them directly |
| XML inclusion (XInclude) / XSLT transforms | use `lxml` externally, then `read_xml` |
