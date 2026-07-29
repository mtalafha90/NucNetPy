# Python API

Everything below is exported at the top level of `nucnetpy` unless noted.
Pages are generated from the docstrings in the source.

## Core objects

```{eval-rst}
.. automodule:: nucnetpy.species
   :members: Species, species_from_za, normalize_species_name, is_massless,
             mass_fraction, abundance_from_mass_fraction

.. automodule:: nucnetpy.core
   :members: Zone, Network

.. automodule:: nucnetpy.reactions
   :members: Reaction, ReactionParticipant, ReactionNetwork, RateFit,
             TabularRate, read_reaclib_text
```

## Evolution

```{eval-rst}
.. automodule:: nucnetpy.solver
   :members: evolve_zone, evolve_network_zones, EvolutionResult, rhs,
             analytic_jacobian, jacobian, jacobian_sparsity,
             constant_thermo, zone_thermo, time_grid
```

## Equilibrium

```{eval-rst}
.. automodule:: nucnetpy.nse
   :members: solve_nse, NSEResult, nse_prefactor, equilibrium_ratio

.. automodule:: nucnetpy.qse
   :members: solve_qse, QSECluster, QSEResult, cluster_abundance, cluster_ydot

.. automodule:: nucnetpy.coulomb
   :members: nse_correction, gamma_e, species_coulomb_chemical_potential,
             species_coulomb_energy, species_coulomb_entropy,
             coulomb_entropy_per_nucleon
```

## Detailed balance

```{eval-rst}
.. automodule:: nucnetpy.detailed_balance
   :members: consistent_reverse_network, reverse_rate, reverse_reaction,
             log_equilibrium_constant, net_flows
```

## Screening, weak rates, decays

```{eval-rst}
.. automodule:: nucnetpy.screening
   :members: SkyNetScreening, ScreeningContext, weak_screening_factor,
             graboske_intermediate_factor, reaction_screening_factor,
             debye_radius, ion_strength

.. automodule:: nucnetpy.weak
   :members: WeakRateTable, read_weak_table, compute_yedot

.. automodule:: nucnetpy.decay
   :members: DecayRecord, decay_constant_from_half_life, add_decay_records,
             fission_reaction

.. automodule:: nucnetpy.neutrino
   :members: NeutrinoLuminosity, NeutrinoQuantity, geometric_flux_rate
```

## Analysis

```{eval-rst}
.. automodule:: nucnetpy.analysis
   :members: largest_mass_fractions, element_abundances, abundance_moment,
             abundances_vs_nucleon_number, flows, ydot,
             nuclear_energy_generation_rate, nuclear_energy_release,
             energy_generation_rate, entropy_generation_rate,
             reaction_entropy_changes, charge_changing_flows,
             system_timescales, heavy_nuclei_abundance, neutron_exposure,
             integrated_currents, separation_energy, species_history,
             compare_rates
```

## Input and output

```{eval-rst}
.. automodule:: nucnetpy.io.xml
   :members: read_xml, write_xml, read_xml_string, read_network_xml,
             write_zone_xml

.. automodule:: nucnetpy.io.jina
   :members: read_jina_xml, combine_jina_xml, jina_database_summary
```

## Utilities

```{eval-rst}
.. automodule:: nucnetpy.network_limiter
   :members: select_species, limit_network

.. automodule:: nucnetpy.hydro
   :members: Trajectory, read_trajectory, exponential_expansion

.. automodule:: nucnetpy.thermo
   :members:

.. automodule:: nucnetpy.validation
   :members: validate_network, validate_zone, regression_summary,
             ValidationIssue

.. automodule:: nucnetpy.rate_modifiers
   :members: RateModifierRegistry, constant_factor, exp_temperature_factor

.. automodule:: nucnetpy.graph
   :members:

.. automodule:: nucnetpy.constants
```
