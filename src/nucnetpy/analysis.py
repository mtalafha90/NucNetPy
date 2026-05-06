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
    sp = species_map.get(normalize_species_name(species_name), Species.parse(species_name))
    if particle.lower() == "n":
        daughter_name = Species.parse(f"{sp.element}{sp.a-1}").name if sp.a > 1 else None
        particle_name = "n"
    else:
        daughter_name = Species.parse(f"{sp.element}{sp.a-1}").name if sp.a > 1 else None
        particle_name = "h1"
    if not daughter_name or daughter_name not in species_map or particle_name not in species_map:
        return None
    return species_map[daughter_name].mass_excess + species_map[particle_name].mass_excess - sp.mass_excess
