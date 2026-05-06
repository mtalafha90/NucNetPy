"""Electron-screening corrections for thermonuclear reaction rates.

The original NucNet Tools exposes user-defined screening callbacks.  This
module provides a Python callback-compatible implementation with common weak
and intermediate-screening approximations.  The API returns multiplicative
factors applied to bare nuclear rates.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional
import math

from .constants import KB_CGS
from .species import Species, normalize_species_name

# cgs electrostatic constants
_E_CHARGE_ESU = 4.803204712570263e-10


def electron_number_density(rho: float, ye: float) -> float:
    from .constants import AVOGADRO
    return max(float(rho) * float(ye) * AVOGADRO, 0.0)


def debye_radius(t9: float, rho: float, ye: float, ion_strength: float = 0.0) -> float:
    """Return Debye radius in cm.

    ``ion_strength`` may be supplied as ``sum Z_i^2 Y_i``.  If omitted, the
    electron term alone is used, which is the usual safe weak-screening fallback.
    """
    t = max(float(t9) * 1.0e9, 1.0)
    ne_factor = max(float(rho) * max(float(ye), 0.0), 0.0)
    plasma = ne_factor + max(float(ion_strength), 0.0) * max(float(rho), 0.0)
    if plasma <= 0:
        return float("inf")
    from .constants import AVOGADRO
    denom = 4.0 * math.pi * (_E_CHARGE_ESU ** 2) * AVOGADRO * plasma
    return math.sqrt(KB_CGS * t / denom)


def ion_strength(abundances: Mapping[str, float], species_map: Optional[Mapping[str, Species]] = None) -> float:
    total = 0.0
    for name, y in abundances.items():
        try:
            sp = species_map.get(normalize_species_name(name)) if species_map else Species.parse(name)
            if sp is None:
                sp = Species.parse(name)
            total += sp.z * sp.z * float(y)
        except Exception:
            continue
    return float(total)


@dataclass
class ScreeningContext:
    t9: float
    rho: float
    ye: float
    abundances: Mapping[str, float]
    species_map: Optional[Mapping[str, Species]] = None

    @property
    def ion_strength(self) -> float:
        return ion_strength(self.abundances, self.species_map)


def weak_screening_factor(z1: int, z2: int, t9: float, rho: float, ye: float, ion_strength_value: float = 0.0) -> float:
    """Salpeter weak-screening enhancement factor.

    The factor is ``exp(Z1 Z2 e^2 / (kT R_D))``.  A cap prevents numerical
    overflow while retaining monotonic enhancement.
    """
    if z1 == 0 or z2 == 0:
        return 1.0
    rd = debye_radius(t9, rho, ye, ion_strength=ion_strength_value)
    if not math.isfinite(rd) or rd <= 0:
        return 1.0
    kt = KB_CGS * max(float(t9) * 1e9, 1.0)
    h = int(z1) * int(z2) * (_E_CHARGE_ESU ** 2) / (kt * rd)
    return float(math.exp(max(min(h, 200.0), -200.0)))


def graboske_intermediate_factor(z1: int, z2: int, t9: float, rho: float, ye: float, ion_strength_value: float = 0.0) -> float:
    """Smooth weak-to-intermediate approximation inspired by Graboske et al.

    This is not tied to one compiled NucNet user file; it gives a stable,
    callback-compatible enhancement and can be replaced by a user function when
    exact project-specific screening is required.
    """
    weak = math.log(weak_screening_factor(z1, z2, t9, rho, ye, ion_strength_value))
    gamma = 0.188 * (max(float(rho) * max(float(ye), 1e-30), 1e-99) ** (1.0/3.0)) / max(float(t9), 1e-30)
    inter = weak * (1.0 + 0.25 * min(gamma, 10.0))
    return float(math.exp(max(min(inter, 200.0), -200.0)))


def reaction_screening_factor(reaction, context: ScreeningContext, model: str = "weak") -> float:
    charged = []
    for p in reaction.reactants:
        try:
            sp = context.species_map.get(p.species) if context.species_map else Species.parse(p.species)
            if sp is None:
                sp = Species.parse(p.species)
            for _ in range(p.count):
                if sp.z > 0:
                    charged.append(sp.z)
        except Exception:
            continue
    if len(charged) < 2:
        return 1.0
    fac = 1.0
    strength = context.ion_strength
    for i in range(len(charged)):
        for j in range(i + 1, len(charged)):
            if model.lower() in {"graboske", "intermediate"}:
                fac *= graboske_intermediate_factor(charged[i], charged[j], context.t9, context.rho, context.ye, strength)
            else:
                fac *= weak_screening_factor(charged[i], charged[j], context.t9, context.rho, context.ye, strength)
    return float(fac)
