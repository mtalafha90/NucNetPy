"""Nuclear decay helpers ported into Python."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Optional
import math
from .reactions import Reaction, ReactionNetwork, ReactionParticipant
from .species import Species, normalize_species_name

LN2 = math.log(2.0)

@dataclass
class DecayRecord:
    parent: str
    daughter: str
    half_life: float
    branch: float = 1.0
    mode: str = "decay"

    @property
    def rate(self) -> float:
        if self.half_life <= 0:
            return 0.0
        return self.branch * LN2 / self.half_life

    def reaction(self) -> Reaction:
        return Reaction([ReactionParticipant(self.parent)], [ReactionParticipant(self.daughter)], constant_rate=self.rate, source="decay", label=self.mode)


def add_decay_records(network, records: Iterable[DecayRecord]) -> None:
    for rec in records:
        network.reactions.add(rec.reaction())


def decay_constant_from_half_life(half_life: float) -> float:
    return 0.0 if half_life <= 0 else LN2 / float(half_life)


def fission_reaction(parent: str, fragments: Iterable[str], neutrons: int = 0, rate: float = 0.0, mode: str = "sf") -> Reaction:
    """Return a fission reaction ``parent -> fragments + neutrons * n``.

    ``rate`` is the fission decay constant in 1/s (use
    :func:`decay_constant_from_half_life` for spontaneous-fission half-lives).
    ``mode`` labels the channel, e.g. ``"sf"`` (spontaneous) or ``"bdf"``
    (beta-delayed); for neutron-induced fission build a normal two-reactant
    Reaction instead.  Blog workflow: "Adding fission to an r-process
    calculation".
    """
    products = [ReactionParticipant(f) for f in fragments]
    if int(neutrons) > 0:
        products.append(ReactionParticipant("n", int(neutrons)))
    return Reaction([ReactionParticipant(parent)], products,
                    constant_rate=float(rate), source="fission", label=str(mode))
