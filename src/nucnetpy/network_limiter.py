"""Network limiting utilities for reduced-network calculations."""
from __future__ import annotations
from typing import Iterable, Set
from .species import normalize_species_name


def limit_network(network, species: Iterable[str], include_linked: bool = False):
    keep: Set[str] = {normalize_species_name(s) for s in species}
    if include_linked:
        changed = True
        while changed:
            changed = False
            for r in network.reactions.reactions:
                names = {p.species for p in (r.reactants + r.products)}
                if names & keep and not names <= keep:
                    keep |= names; changed = True
    network.species = {k:v for k,v in network.species.items() if k in keep}
    network.reactions.reactions = [r for r in network.reactions.reactions if all(p.species in keep for p in r.reactants + r.products)]
    for z in network.zones:
        z.abundances = {k:v for k,v in z.abundances.items() if k in keep}
    return network
