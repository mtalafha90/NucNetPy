"""Analysis helpers equivalent to common NucNet Tools example programs."""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np

from .core import Network, Zone
from .species import Species, normalize_species_name


def largest_mass_fractions(zone: Zone, species_map: Mapping[str, Species], n: int = 10, min_x: float = 0.0):
    items = [(name, x) for name, x in zone.mass_fractions(species_map).items() if x >= min_x]
    return sorted(items, key=lambda kv: kv[1], reverse=True)[:n]


def abundances_for_element(zone: Zone, element: str, species_map: Mapping[str, Species]) -> Dict[str, float]:
    element = element.lower()
    return {name: y for name, y in zone.abundances.items() if species_map.get(name, Species.parse(name)).element == element}


def abundances_vs_nucleon_number(zone: Zone, species_map: Mapping[str, Species]) -> Dict[int, float]:
    out = defaultdict(float)
    for name, y in zone.abundances.items():
        out[species_map.get(name, Species.parse(name)).a] += y
    return dict(sorted(out.items()))


def element_abundances(zone: Zone, species_map: Mapping[str, Species]) -> Dict[str, float]:
    out = defaultdict(float)
    for name, y in zone.abundances.items():
        out[species_map.get(name, Species.parse(name)).element] += y
    return dict(sorted(out.items()))


def abundance_moment(zone: Zone, species_map: Mapping[str, Species], moment: int = 1, kind: str = "a") -> float:
    total = 0.0
    for name, y in zone.abundances.items():
        sp = species_map.get(name, Species.parse(name))
        q = sp.a if kind.lower() == "a" else sp.z if kind.lower() == "z" else sp.n
        total += (q ** moment) * y
    return float(total)


def species_history(network: Network, species: str) -> List[Tuple[int, Tuple[str, str, str], float]]:
    key = normalize_species_name(species)
    return [(i, z.label, z.get_abundance(key)) for i, z in enumerate(network.zones)]


def flows(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None):
    z = network.zone(zone_index)
    return network.reactions.flows(z.abundances, t9=t9 or z.temperature9(), rho=rho or z.density())


def ydot(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None):
    z = network.zone(zone_index)
    return network.reactions.ydot(z.abundances, t9=t9 or z.temperature9(), rho=rho or z.density())


def energy_generation_rate(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> float:
    # MeV per nucleon-ish proxy: sum Q * reaction flux. Users should calibrate
    # to their original libnucnet build for publication-grade energetics.
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    total = 0.0
    for r in network.reactions.reactions:
        total += r.q_value * r.flux(z.abundances, t9=t9, rho=rho)
    return float(total)


def compare_rates(net_a: Network, net_b: Network, t9: float, rho: float = 1.0) -> List[Tuple[str, float, float, float]]:
    ra = net_a.reactions.rates(t9, rho)
    rb = net_b.reactions.rates(t9, rho)
    keys = sorted(set(ra) | set(rb))
    return [(k, ra.get(k, 0.0), rb.get(k, 0.0), rb.get(k, 0.0) - ra.get(k, 0.0)) for k in keys]


def separation_energy(species_name: str, species_map: Mapping[str, Species], particle: str = "n") -> Optional[float]:
    """Return S_n or S_p in MeV, or None if the needed nuclides are absent.

    ``S_n(Z, A) = ME(Z, A-1) + ME(n) - ME(Z, A)`` and
    ``S_p(Z, A) = ME(Z-1, A-1) + ME(p) - ME(Z, A)``.
    """
    from .species import species_from_za
    sp = species_map.get(normalize_species_name(species_name), Species.parse(species_name))
    if sp.a <= 1:
        return None
    try:
        if particle.lower() == "n":
            daughter_name = species_from_za(sp.z, sp.a - 1).name
            particle_name = "n"
        else:
            if sp.z < 1:
                return None
            daughter_name = species_from_za(sp.z - 1, sp.a - 1).name
            particle_name = "h1"
    except ValueError:
        return None
    if daughter_name not in species_map or particle_name not in species_map:
        return None
    return species_map[daughter_name].mass_excess + species_map[particle_name].mass_excess - sp.mass_excess


def charge_changing_flows(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> Dict[str, float]:
    """Return each reaction's contribution to dYe/dt (flux times net ΔZ).

    Only reactions with a nonzero nuclear charge change (weak reactions such as
    beta decays and electron captures) appear in the result; strong reactions
    conserve Z and contribute nothing.  ``sum(result.values())`` is dYe/dt.
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    out: Dict[str, float] = {}
    for r in network.reactions.reactions:
        dz = 0
        for name, coeff in r.stoichiometry().items():
            sp = network.species.get(name)
            if sp is None:
                try:
                    sp = Species.parse(name)
                except Exception:
                    continue
            if sp.a > 0 and sp.z >= 0:
                dz += coeff * sp.z
        if dz != 0:
            out[r.string] = dz * r.flux(z.abundances, t9=t9, rho=rho)
    return out


def system_timescales(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> Dict[str, float]:
    """Return per-species timescales ``Y / |dY/dt|`` in seconds.

    Species with zero derivative get ``inf``.  The shortest timescales identify
    the stiffest components of the system (blog: "Computing system timescales").
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    dy = network.reactions.ydot(z.abundances, t9=t9, rho=rho)
    out: Dict[str, float] = {}
    for name in set(z.abundances) | set(dy):
        y = z.get_abundance(name)
        rate = abs(dy.get(name, 0.0))
        out[name] = float(y / rate) if rate > 0.0 else float("inf")
    return out


def heavy_nuclei_abundance(zone: Zone, species_map: Mapping[str, Species], zmin: int = 3) -> float:
    """Return Y_h, the total abundance of nuclei with Z >= ``zmin``.

    The photon-to-heavy-nucleus ratio and the r-process neutron-to-seed ratio
    are built from this quantity (blog: "Computing the number of heavy nuclei").
    """
    total = 0.0
    for name, y in zone.abundances.items():
        sp = species_map.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        if sp.z >= zmin:
            total += y
    return float(total)


def neutron_exposure(result, thermo) -> float:
    """Return the s-process neutron exposure tau = ∫ n_n v_T dt in mb^-1.

    ``result`` is an :class:`~nucnetpy.solver.EvolutionResult` whose species
    include the neutron; ``thermo`` is the same ``(t, abundances) -> (t9, rho)``
    function used for the evolution.  ``n_n = rho N_A Y_n`` and the thermal
    velocity is ``v_T = sqrt(2 k T / m_n)`` (blog: "Computing the s-process
    neutron exposure").
    """
    from .constants import AVOGADRO, KB_CGS, MN_G
    if "n" not in result.species:
        return 0.0
    j = result.species.index("n")
    integrand = np.zeros(len(result.time))
    for i, t in enumerate(result.time):
        abund = {s: float(v) for s, v in zip(result.species, result.y[i])}
        t9, rho = thermo(float(t), abund)
        n_n = max(float(rho), 0.0) * AVOGADRO * max(float(result.y[i, j]), 0.0)
        v_t = np.sqrt(2.0 * KB_CGS * max(float(t9), 1e-30) * 1.0e9 / MN_G)
        integrand[i] = n_n * v_t
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 compat
    tau_cm2 = float(trapezoid(integrand, np.asarray(result.time, dtype=float)))
    return tau_cm2 * 1.0e-27  # cm^-2 -> mb^-1


def reaction_entropy_changes(network: Network, abundances: Mapping[str, float], t9: float, rho: float) -> Dict[str, float]:
    """Return per-reaction entropy change ΔS in units of k_B per reaction.

    ``ΔS = -sum_i nu_i mu_i/kT`` with Maxwell–Boltzmann chemical potentials
    ``mu_i/kT = ln Y_i - ln pref_i``; this equals the C++ per-reaction
    ``Q/kT + sum_r ln(Y/Y_Q) - sum_p ln(Y/Y_Q)`` of NucNet Tools
    ``flow_utilities.cpp`` and, under detailed balance, ``ln(f/r)``.
    """
    from .nse import _log_prefactor
    import math
    out: Dict[str, float] = {}
    for r in network.reactions.reactions:
        ds = 0.0
        for name, nu in r.stoichiometry().items():
            sp = network.species.get(name)
            if sp is None:
                try:
                    sp = Species.parse(name)
                except Exception:
                    continue
            if sp.a <= 0:
                continue
            y = max(float(abundances.get(name, 0.0)), 1e-300)
            ds -= nu * (math.log(y) - _log_prefactor(sp, t9, rho))
        out[r.string] = float(ds)
    return out


def entropy_generation_rate(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None, use_reverse: bool = True) -> float:
    """Return dS/dt per nucleon in units of k_B per second.

    Ports NucNet Tools ``compute_entropy_generation_rate``: the sum over
    reactions of ``(f - r) * ΔS`` with ``r`` the detailed-balance reverse flux
    and ``ΔS = ln(f/r)`` per reaction, so every term is non-negative and the
    total vanishes at NSE (blog series "Computing the entropy generation
    rate").  With ``use_reverse=False`` only forward fluxes are used
    (``dS/dt = -sum_i (mu_i/kT) dY_i/dt``).  Electron/neutrino chemical
    potential terms of the C++ version are not included; nucnetpy networks
    carry weak reactions with their own tabulated rates instead.
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    ds = reaction_entropy_changes(network, z.abundances, t9, rho)
    if use_reverse:
        from .detailed_balance import net_flows
        flows_frn = net_flows(network, z.abundances, t9=t9, rho=rho)
        return float(sum(flows_frn[key][2] * ds[key] for key in ds))
    total = 0.0
    for r in network.reactions.reactions:
        total += r.flux(z.abundances, t9=t9, rho=rho) * ds[r.string]
    return float(total)


def integrated_currents(network: Network, result, thermo, use_reverse: bool = True) -> Dict[str, float]:
    """Return per-reaction time-integrated net currents over an evolution.

    Ports NucNet Tools ``update_flow_currents``: for each reaction the current
    accumulates ``(forward - reverse) * dt`` across the time grid of
    ``result`` (an :class:`~nucnetpy.solver.EvolutionResult`), integrated here
    with the trapezoidal rule.  The integrated current of a reaction is the
    net number of transitions per nucleon it produced over the calculation
    (blog: "Creating integrated currents diagrams", "Analyzing integrated
    currents quantitatively").
    """
    from .detailed_balance import net_flows
    times = np.asarray(result.time, dtype=float)
    strings = [r.string for r in network.reactions.reactions]
    rates = np.zeros((len(times), len(strings)))
    for i, t in enumerate(times):
        abund = {s: float(v) for s, v in zip(result.species, result.y[i])}
        t9, rho = thermo(float(t), abund)
        if use_reverse:
            frn = net_flows(network, abund, t9=t9, rho=rho)
            rates[i] = [frn[s][2] for s in strings]
        else:
            rates[i] = [r.flux(abund, t9=t9, rho=rho) for r in network.reactions.reactions]
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 compat
    return {s: float(trapezoid(rates[:, j], times)) for j, s in enumerate(strings)}
