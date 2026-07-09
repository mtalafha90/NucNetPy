"""nucnetpy: pure-Python replacement layer for NucNet Tools workflows."""
from .species import Species, species_from_za, normalize_species_name
from .core import Network, Zone
from .reactions import Reaction, ReactionParticipant, ReactionNetwork, RateFit, TabularRate, read_reaclib_text
from .solver import evolve_zone, evolve_network_zones, constant_thermo, time_grid, EvolutionResult, jacobian, jacobian_sparsity
from .io.xml import read_xml, write_xml, read_network_xml, write_zone_xml, read_xml_string
from .io.jina import read_jina_xml, combine_jina_xml, jina_database_summary

__version__ = "1.0.0"

from .nse import solve_nse, NSEResult, equilibrium_ratio
from .screening import ScreeningContext, weak_screening_factor, reaction_screening_factor, graboske_intermediate_factor
from .weak import WeakRateTable, read_weak_table, compute_yedot
from .validation import validate_network, validate_zone, regression_summary

from .mathutils import linear_interpolation, bilinear_interpolation, two_d_interpolation
from .decay import DecayRecord, decay_constant_from_half_life, add_decay_records, fission_reaction
from .hydro import Trajectory, read_trajectory, exponential_expansion
from .matrix_solver import solve_linear
from .network_limiter import limit_network, select_species
from .neutrino import NeutrinoLuminosity, NeutrinoQuantity, geometric_flux_rate
from .rate_modifiers import RateModifierRegistry, constant_factor, exp_temperature_factor
from .detailed_balance import log_equilibrium_constant, reverse_rate, reverse_reaction, net_flows
from .analysis import (charge_changing_flows, system_timescales, heavy_nuclei_abundance,
                       neutron_exposure, entropy_generation_rate, separation_energy,
                       integrated_currents, reaction_entropy_changes,
                       flows, ydot, energy_generation_rate, largest_mass_fractions,
                       element_abundances, abundance_moment, abundances_vs_nucleon_number)
from .coulomb import (nse_correction as coulomb_nse_correction, gamma_e,
                      species_coulomb_chemical_potential, species_coulomb_energy,
                      species_coulomb_entropy, coulomb_entropy_per_nucleon)
from .qse import solve_qse, QSECluster, QSEResult, cluster_abundance, cluster_ydot
