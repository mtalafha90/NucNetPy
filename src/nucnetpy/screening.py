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


# --------------------------------------------------------------------------
# Composition-dependent screening after SkyNet (Lippuner & Roberts 2017)
# --------------------------------------------------------------------------

_LN25 = math.log(25.0)
#: Exponent of the intermediate-screening interpolation of Graboske et al.
_B_INTERMEDIATE = 0.860
#: Prefactor of the plasma parameter lambda_0 in the units used here (baryon
#: number density in cm^-3, temperature in MeV).
_LAMBDA0_FACTOR = 1.9370131349470739099e-19
_KB_MEV_PER_GK = 8.6173324e-2


class SkyNetScreening:
    """Screening as a Coulomb chemical potential per nuclear charge.

    Instead of a pairwise enhancement factor, each charge ``Z`` is assigned a
    Coulomb chemical potential ``mu(Z)/kT`` in the plasma, and a reaction is
    enhanced by

        exp( sum_reactants mu(Z_i)/kT  -  mu(sum_i Z_i)/kT ),

    the difference between the potentials of the separated reactants and of the
    fused compound charge.  This has two advantages over a pairwise factor: it
    extends naturally to reactions with more than two charged reactants, and
    because it is a chemical potential it can be applied to the forward and
    reverse directions consistently, so screening does not break detailed
    balance.

    ``mu(Z)`` blends three regimes with smooth ``tanh`` switches: weak
    screening, ``mu_w = -Z^2 zeta lambda_0 / 2``; intermediate screening,
    ``mu_i = -0.380 lambda_0^b eta_b Z^(b+1)`` with ``b = 0.860``; and strong
    screening from a one-component-plasma fit.

    The electron-positron pair contribution to ``zeta``, which requires the
    electron degeneracy parameter from an equation of state, is not included.
    Pass ``pair_term`` if it is available.  Neglecting it is accurate while
    pairs are unimportant, which covers ordinary stellar and explosive burning,
    but it should not be relied on where positrons are abundant.

    The object is a callable suitable for the ``screening`` argument of
    :func:`nucnetpy.evolve_zone`.  It caches ``mu(Z)`` and recomputes it in
    :meth:`update`, which the network calls once per right-hand-side
    evaluation, so the cost does not scale with the number of reactions.
    """

    def __init__(self, species_map: Optional[Mapping[str, Species]] = None,
                 max_z: int = 110, pair_term: float = 0.0):
        self.species_map = species_map
        self.max_z = int(max_z)
        self.pair_term = float(pair_term)
        self._mu = [0.0] * (self.max_z + 1)

    def _species(self, name: str) -> Optional[Species]:
        sp = self.species_map.get(normalize_species_name(name)) if self.species_map else None
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                return None
        return sp

    def update(self, abundances: Mapping[str, float], t9: float, rho: float,
               ye: Optional[float] = None) -> None:
        """Recompute ``mu(Z)`` for this composition and thermodynamic state."""
        b = _B_INTERMEDIATE
        sum_y = sum_zy = sum_z2y = sum_z3bm1y = 0.0
        for name, y in abundances.items():
            sp = self._species(name)
            if sp is None or sp.a <= 0:
                continue
            y = max(float(y), 0.0)
            z = float(sp.z)
            sum_y += y
            sum_zy += z * y
            sum_z2y += z * z * y
            if z > 0.0:
                sum_z3bm1y += (z ** (3.0 * b - 1.0)) * y
        if sum_y <= 0.0:
            self._mu = [0.0] * (self.max_z + 1)
            return

        from .constants import AVOGADRO
        zbar = sum_zy / sum_y
        nb = max(float(rho), 1e-300) * AVOGADRO
        t_mev = max(float(t9), 1e-30) * _KB_MEV_PER_GK
        lambda0 = _LAMBDA0_FACTOR * math.sqrt(nb * sum_y) / (t_mev ** 1.5)
        zeta = math.sqrt(max(sum_z2y + self.pair_term, 0.0) / sum_y)
        if zeta <= 0.0 or zbar <= 0.0:
            self._mu = [0.0] * (self.max_z + 1)
            return
        etab = (sum_z3bm1y / sum_y) / ((zeta ** (3.0 * b - 2.0))
                                       * (zbar ** (2.0 - 2.0 * b)))

        mu = [0.0] * (self.max_z + 1)
        for z_int in range(1, self.max_z + 1):
            z = float(z_int)
            p = (zeta + z) * zeta * zeta * lambda0
            two_ln_p = 2.0 * math.log(max(p, 1e-300))
            f_weak = 0.5 * (math.tanh(-two_ln_p - _LN25) + 1.0)
            f_strong = 0.5 * (math.tanh(two_ln_p - _LN25) + 1.0)
            f_inter = 1.0 - f_weak - f_strong
            mu_weak = -0.5 * z * z * zeta * lambda0
            mu_inter = -0.380 * (lambda0 ** b) * etab * (z ** (b + 1.0))
            mu_strong = (-(lambda0 ** (2.0 / 3.0))
                         * (0.6240 * (z ** (5.0 / 3.0)) * (zbar ** (1.0 / 3.0))
                            + 0.1971 * (z ** (4.0 / 3.0)) * (zbar ** (2.0 / 3.0))
                            - 0.0374 * z * zbar)
                         + 9.0 / 16.0 * z / zbar
                         - 0.4600 * ((z / zbar) ** (2.0 / 3.0)))
            mu[z_int] = f_weak * mu_weak + f_inter * mu_inter + f_strong * mu_strong
        self._mu = mu

    def chemical_potential(self, z: int) -> float:
        """Return ``mu(Z)/kT`` for charge ``z``."""
        z = int(z)
        if z <= 0:
            return 0.0
        return self._mu[min(z, self.max_z)]

    def factor(self, reaction) -> float:
        """Return the multiplicative screening enhancement of ``reaction``."""
        total_z = 0
        total_mu = 0.0
        for p in reaction.nuclear_reactants:
            sp = self._species(p.species)
            if sp is None or sp.z <= 0:
                continue
            total_mu += p.count * self.chemical_potential(sp.z)
            total_z += p.count * sp.z
        if total_z <= 0:
            return 1.0
        exponent = total_mu - self.chemical_potential(total_z)
        return float(math.exp(max(min(exponent, 200.0), -200.0)))

    def __call__(self, reaction, t9: float = 0.0, rho: float = 0.0,
                 ye: Optional[float] = None) -> float:
        return self.factor(reaction)
