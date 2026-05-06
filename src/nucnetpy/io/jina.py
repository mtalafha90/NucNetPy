"""JINA/libnucnet XML database helpers.

The original NucNet Tools workflow usually reads nuclide data and reaction data
from JINA ReacLib/libnucnet XML files.  These helpers keep that workflow in the
pure-Python port: one file may contain nuclides, another may contain reactions,
and both can be merged into one :class:`nucnetpy.core.Network` object.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional

from ..core import Network
from ..reactions import ReactionNetwork
from .xml import _parse_species, _parse_reactions, _parse_zones, write_xml


def read_jina_xml(
    nuclides_xml: str | Path,
    reactions_xml: str | Path,
    zones_xml: Optional[str | Path] = None,
) -> Network:
    """Read separate JINA/libnucnet nuclide and reaction XML files.

    Parameters
    ----------
    nuclides_xml:
        XML file containing ``<nuclide>`` records with Z, A, mass excess,
        spin, and optional partition-function data.
    reactions_xml:
        XML file containing ``<reaction>`` records and ReacLib rate fits.
    zones_xml:
        Optional XML file containing zones/abundances.  If omitted, the returned
        network contains species and reactions but no zones.

    Returns
    -------
    Network
        A single Python network object with species, reactions, and optionally
        zones.
    """
    net = Network(metadata={
        "nuclides_xml": str(nuclides_xml),
        "reactions_xml": str(reactions_xml),
        "zones_xml": str(zones_xml) if zones_xml is not None else "",
        "database": "JINA/libnucnet XML",
    })

    nuc_root = ET.parse(nuclides_xml).getroot()
    reac_root = ET.parse(reactions_xml).getroot()

    _parse_species(nuc_root, net)
    _parse_reactions(reac_root, net)

    # Ensure every species referenced by reactions has a nuclide entry.  If the
    # nuclide XML is complete, this will not create anything new; if a light
    # species alias is missing, it gives the user a useful placeholder.
    for name in net.reactions.species_names():
        try:
            net.ensure_species(name)
        except Exception:
            pass

    if zones_xml is not None:
        zone_root = ET.parse(zones_xml).getroot()
        _parse_zones(zone_root, net)

    return net


def combine_jina_xml(
    nuclides_xml: str | Path,
    reactions_xml: str | Path,
    output_xml: str | Path,
    zones_xml: Optional[str | Path] = None,
) -> Network:
    """Read separate JINA XML files and write one combined nucnetpy XML file."""
    net = read_jina_xml(nuclides_xml, reactions_xml, zones_xml=zones_xml)
    write_xml(net, output_xml)
    return net


def jina_database_summary(net_or_nuclides, reactions_xml=None, zones_xml=None) -> dict[str, int]:
    """Return a compact summary of a loaded JINA database.

    Accepts either an already loaded Network, or file paths
    ``(nuclides_xml, reactions_xml[, zones_xml])``.
    """
    if isinstance(net_or_nuclides, Network):
        net = net_or_nuclides
    else:
        net = read_jina_xml(net_or_nuclides, reactions_xml, zones_xml)
    return {
        "species": len(net.species),
        "reactions": len(net.reactions.reactions),
        "zones": len(net.zones),
        "invalid_reactions": len(net.reactions.invalid_reactions(net.species)),
        "missing_species": len(net.validate()["missing_species"]),
    }
