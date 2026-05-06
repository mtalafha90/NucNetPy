"""Core containers: Network, Zone, and multi-zone collections."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np

from .species import Species, normalize_species_name, mass_fraction
from .reactions import ReactionNetwork

@dataclass
class Zone:
    label: Tuple[str, str, str] = ("0", "0", "0")
    abundances: Dict[str, float] = field(default_factory=dict)
    properties: Dict[str, str] = field(default_factory=dict)
    optional_properties: Dict[str, str] = field(default_factory=dict)

    def set_abundance(self, species: str, y: float) -> None:
        self.abundances[normalize_species_name(species)] = float(y)

    def get_abundance(self, species: str, default: float = 0.0) -> float:
        return float(self.abundances.get(normalize_species_name(species), default))

    def mass_fractions(self, species_map: Optional[Mapping[str, Species]] = None) -> Dict[str, float]:
        out = {}
        for name, y in self.abundances.items():
            a = species_map[name].a if species_map and name in species_map else Species.parse(name).a
            out[name] = mass_fraction(y, a)
        return out

    def ye(self, species_map: Optional[Mapping[str, Species]] = None) -> float:
        total = 0.0
        for name, y in self.abundances.items():
            sp = species_map[name] if species_map and name in species_map else Species.parse(name)
            total += sp.z * y
        return float(total)

    def ysum(self) -> float:
        return float(sum(self.abundances.values()))

    def xsum(self, species_map: Optional[Mapping[str, Species]] = None) -> float:
        return float(sum(self.mass_fractions(species_map).values()))

    def temperature9(self, default: float = 1.0) -> float:
        props = {**self.optional_properties, **self.properties}
        for key in ("t9", "T9", "temperature9", "temperature", "t9_0"):
            if key in props:
                val = float(props[key])
                return val / 1.0e9 if key == "temperature" and val > 1e6 else val
        return default

    def density(self, default: float = 1.0) -> float:
        props = {**self.optional_properties, **self.properties}
        for key in ("rho", "density", "rho_0"):
            if key in props:
                return float(props[key])
        return default

@dataclass
class Network:
    species: Dict[str, Species] = field(default_factory=dict)
    reactions: ReactionNetwork = field(default_factory=ReactionNetwork)
    zones: List[Zone] = field(default_factory=list)
    metadata: Dict[str, str] = field(default_factory=dict)

    def add_species(self, species: Species) -> None:
        self.species[normalize_species_name(species.name)] = species

    def ensure_species(self, name: str) -> Species:
        key = normalize_species_name(name)
        if key not in self.species:
            self.species[key] = Species.parse(key)
        return self.species[key]

    def add_zone(self, zone: Zone) -> None:
        self.zones.append(zone)
        for name in zone.abundances:
            try:
                self.ensure_species(name)
            except Exception:
                pass

    def species_names(self) -> List[str]:
        names = set(self.species)
        names.update(self.reactions.species_names())
        for z in self.zones:
            names.update(z.abundances)
        return sorted(names, key=lambda n: (self.species.get(n, Species.parse(n) if n != '' else Species('n',0,1)).z, self.species.get(n, Species.parse(n) if n != '' else Species('n',0,1)).a))

    def zone(self, index: int = 0) -> Zone:
        return self.zones[index]

    def validate(self) -> Dict[str, object]:
        missing = set()
        for r in self.reactions.reactions:
            for p in r.reactants + r.products:
                if p.species not in self.species:
                    missing.add(p.species)
        invalid = self.reactions.invalid_reactions(self.species)
        return {"missing_species": sorted(missing), "invalid_reactions": invalid, "n_species": len(self.species), "n_reactions": len(self.reactions.reactions), "n_zones": len(self.zones)}

    def abundance_matrix(self, species: Optional[Sequence[str]] = None) -> Tuple[List[str], np.ndarray]:
        names = [normalize_species_name(s) for s in (species or self.species_names())]
        arr = np.zeros((len(self.zones), len(names)))
        for i, zone in enumerate(self.zones):
            for j, name in enumerate(names):
                arr[i, j] = zone.get_abundance(name)
        return names, arr
