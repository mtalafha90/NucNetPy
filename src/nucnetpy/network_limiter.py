"""Network limiting utilities for reduced-network calculations."""
from __future__ import annotations
from typing import Iterable, List, Optional, Set
from .species import normalize_species_name


def select_species(network, zmin: Optional[int] = None, zmax: Optional[int] = None,
                   amin: Optional[int] = None, amax: Optional[int] = None,
                   elements: Optional[Iterable[str]] = None) -> List[str]:
    """Return species names matching Z/A ranges and/or element symbols.

    A pure-Python replacement for the libnucnet nuclide XPath selections used
    in the blog workflows "Selecting an input network" and "Applying a nuclide
    XPath expression".  Combine with :func:`limit_network`::

        limit_network(net, select_species(net, zmin=1, zmax=28, amax=60))
    """
    wanted = {e.lower() for e in elements} if elements else None
    out = []
    for name, sp in network.species.items():
        if zmin is not None and sp.z < zmin: continue
        if zmax is not None and sp.z > zmax: continue
        if amin is not None and sp.a < amin: continue
        if amax is not None and sp.a > amax: continue
        if wanted is not None and sp.element not in wanted: continue
        out.append(name)
    return sorted(out, key=lambda n: (network.species[n].z, network.species[n].a))


def limit_network(network, species: Iterable[str], include_linked: bool = False):
    """Restrict ``network`` in place to the given nuclides and the reactions among them.

    A reaction survives if every nuclide it involves is kept.  Photons and
    leptons are not nuclides, so they never decide the outcome: a
    photodisintegration among kept nuclides survives whether or not ``"gamma"``
    is listed, and listing it is harmless.  With ``include_linked`` the kept set
    grows to every nuclide that shares a reaction with it; the photon does not
    count as a link, or every photon reaction would join the network.
    """
    from .species import is_massless

    def nuclides(reaction) -> Set[str]:
        return {p.species for p in reaction.reactants + reaction.products
                if not is_massless(p.species)}

    keep: Set[str] = {normalize_species_name(s) for s in species if not is_massless(s)}
    if include_linked:
        changed = True
        while changed:
            changed = False
            for r in network.reactions.reactions:
                names = nuclides(r)
                if names & keep and not names <= keep:
                    keep |= names; changed = True
    network.species = {k: v for k, v in network.species.items() if k in keep or is_massless(k)}
    network.reactions.reactions = [r for r in network.reactions.reactions if nuclides(r) <= keep]
    for z in network.zones:
        z.abundances = {k: v for k, v in z.abundances.items() if k in keep or is_massless(k)}
    return network
