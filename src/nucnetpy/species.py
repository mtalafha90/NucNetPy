"""Nuclide/species utilities.

This module provides the most frequently needed pieces of libnucnet species
handling in pure Python: name parsing, Z/A/N bookkeeping, mass fractions,
charge fractions, and common element symbols.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Dict, Optional

_SYMBOLS = [
    "n", "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al",
    "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co",
    "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr", "nb",
    "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe", "cs",
    "ba", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm",
    "yb", "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi",
    "po", "at", "rn", "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk",
    "cf", "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds", "rg",
    "cn", "nh", "fl", "mc", "lv", "ts", "og",
]
SYMBOL_TO_Z = {s: i for i, s in enumerate(_SYMBOLS)}
Z_TO_SYMBOL = {i: s for i, s in enumerate(_SYMBOLS)}
ALIASES = {"p": "h1", "d": "h2", "t": "h3", "alpha": "he4", "neutron": "n", "photon": "gamma"}
SPECIAL_SPECIES = {
    "gamma": (0, 0),
    "electron": (-1, 0),
    "positron": (1, 0),
    "neutrinoe": (0, 0),
    "antineutrinoe": (0, 0),
    "neutrino_e": (0, 0),
    "anti-neutrino_e": (0, 0),
}

@dataclass(frozen=True)
class Species:
    """A nuclear species.

    Parameters are intentionally similar to libnucnet's Nuc objects.
    `mass_excess` is in MeV when available. `partition` may store temperature
    dependent partition-function data.
    """
    name: str
    z: int
    a: int
    mass_excess: float = 0.0
    spin: Optional[float] = None
    source: Optional[str] = None
    partition: Dict[float, float] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return self.a - self.z

    @property
    def element(self) -> str:
        return Z_TO_SYMBOL.get(self.z, f"z{self.z}")

    @property
    def latex(self) -> str:
        if self.name == "n":
            return r"n"
        sym = self.element.capitalize()
        return rf"^{{{self.a}}}\mathrm{{{sym}}}"

    @classmethod
    def parse(cls, value: str, **kwargs) -> "Species":
        name = normalize_species_name(value)
        if name == "n":
            return cls("n", 0, 1, **kwargs)
        if name in SPECIAL_SPECIES:
            z, a = SPECIAL_SPECIES[name]
            return cls(name, z, a, **kwargs)
        m = re.match(r"^([a-z]+)(\d+)$", name)
        if not m:
            raise ValueError(f"Cannot parse species name {value!r}; expected names like h1, he4, ni56, or n")
        sym, a_s = m.group(1), m.group(2)
        if sym not in SYMBOL_TO_Z:
            raise ValueError(f"Unknown element symbol in species name {value!r}")
        return cls(name, SYMBOL_TO_Z[sym], int(a_s), **kwargs)

def normalize_species_name(value: str) -> str:
    v = str(value).strip().lower().replace("-", "").replace("_", "")
    return ALIASES.get(v, v)

def species_from_za(z: int, a: int, **kwargs) -> Species:
    if z == 0 and a == 1:
        return Species("n", z, a, **kwargs)
    sym = Z_TO_SYMBOL.get(int(z))
    if sym is None:
        raise ValueError(f"Unknown Z={z}")
    return Species(f"{sym}{int(a)}", int(z), int(a), **kwargs)

def mass_fraction(y: float, a: int) -> float:
    return float(y) * int(a)

def abundance_from_mass_fraction(x: float, a: int) -> float:
    a = int(a)
    return 0.0 if a == 0 else float(x) / a
