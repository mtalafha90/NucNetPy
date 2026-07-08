"""Coulomb (plasma) corrections to NSE, ported from NucNet Tools.

This is a direct port of ``user/nse_corr.cpp`` (r647): the Coulomb correction
of Bravo & Garcia-Senz (MNRAS 307, 984, 1999) with the fit constants of
Ogata & Ichimaru (Phys. Rev. A 36, 5451, 1987) and the one-component-plasma
value ``f(Gamma=1)/kT = -0.420`` of Slattery, Doolen & DeWitt (1980).

The species Coulomb chemical potential per kT is

    Gamma_i > 1:  mu_C/kT = a G + 4b G^(1/4) - 4c G^(-1/4) + d ln G - o
    Gamma_i <= 1: mu_C/kT = (beta/gamma) G^gamma - G^(3/2)/sqrt(3)

with ``G = Gamma_i = Z^(5/3) Gamma_e`` and the two branches matched in value
and first derivative at ``Gamma_i = 1``.  ``Gamma_e = e^2 / (a_e k T)`` uses
the electron-cloud radius ``a_e = (4 pi rho N_A Ye / 3)^(-1/3)``.

The NSE correction factor is ``f_corr = -mu_C/kT`` (libnucnet convention): it
is added to the exponent of the NSE abundance of each species and to the
equilibrium constants used for detailed-balance reverse rates.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .constants import AVOGADRO, KB_CGS

# electron charge in esu (matches GSL's CGSM electron charge * speed of light)
_E_ESU = 4.80320471257e-10

# Ogata & Ichimaru (1987) fit constants and Slattery, Doolen & DeWitt f(1)/kT.
_A = -0.898004
_B = 0.9678
_C = 0.22070
_D = -0.86097
_F1KT = -0.420
_O = _A + 4.0 * (_B - _C) - _F1KT
_BETA = _A + _B + _C + _D + math.sqrt(3.0) / 2.0
_GAMMA = _BETA / ((1.0 / math.sqrt(3.0)) + _F1KT)


def gamma_e(t9: float, rho: float, ye: float) -> float:
    """Return the electron Coulomb coupling parameter ``Gamma_e``."""
    ne = max(float(rho), 1e-300) * AVOGADRO * max(float(ye), 1e-300)
    a_e = (4.0 * math.pi * ne / 3.0) ** (-1.0 / 3.0)
    return (_E_ESU ** 2) / (a_e * KB_CGS * max(float(t9), 1e-30) * 1.0e9)


def species_coulomb_chemical_potential(z: int, t9: float, rho: float, ye: float) -> float:
    """Return ``mu_C/kT`` for a species of charge ``z`` (0 for neutrons)."""
    z = int(z)
    if z <= 0:
        return 0.0
    g = gamma_e(t9, rho, ye) * (float(z) ** (5.0 / 3.0))
    if g > 1.0:
        return (_A * g + 4.0 * _B * g ** 0.25 - 4.0 * _C * g ** -0.25
                + _D * math.log(g) - _O)
    return (_BETA / _GAMMA) * g ** _GAMMA - g ** 1.5 / math.sqrt(3.0)


def species_coulomb_energy(z: int, t9: float, rho: float, ye: float) -> float:
    """Return the species Coulomb energy per particle in units of kT."""
    z = int(z)
    if z <= 0:
        return 0.0
    g = gamma_e(t9, rho, ye) * (float(z) ** (5.0 / 3.0))
    if g > 1.0:
        return _A * g + _B * g ** 0.25 + _C * g ** -0.25 + _D
    return _BETA * g ** _GAMMA - math.sqrt(3.0) * g ** 1.5 / 2.0


def species_coulomb_entropy(z: int, t9: float, rho: float, ye: float) -> float:
    """Return the species Coulomb entropy per particle in units of k_B."""
    z = int(z)
    if z <= 0:
        return 0.0
    g = gamma_e(t9, rho, ye) * (float(z) ** (5.0 / 3.0))
    if g > 1.0:
        return (-3.0 * _B * g ** 0.25 + 5.0 * _C * g ** -0.25
                + _D * (1.0 - math.log(g)) + _O)
    return _BETA * g ** _GAMMA * (1.0 - 1.0 / _GAMMA) - g ** 1.5 / (2.0 * math.sqrt(3.0))


def nse_correction(species, t9: float, rho: float, ye: float) -> float:
    """Return the NSE correction factor ``f_corr = -mu_C/kT`` for a species.

    Signature-compatible with the ``nse_correction`` hooks accepted by
    :func:`nucnetpy.nse.solve_nse` and the detailed-balance helpers: pass this
    function (or your own with the same signature) to enable Coulomb
    corrections.
    """
    z = species.z if hasattr(species, "z") else int(species)
    return -species_coulomb_chemical_potential(z, t9, rho, ye)


def coulomb_entropy_per_nucleon(abundances, species_map, t9: float, rho: float, ye: float) -> float:
    """Return the total Coulomb entropy per nucleon in units of k_B.

    ``sum_i Y_i * s_C(Z_i)`` — the plasma part of the entropy used in the blog
    workflow "Computing NSE corrections to the entropy".
    """
    total = 0.0
    for name, y in abundances.items():
        sp = species_map.get(name)
        if sp is None:
            continue
        total += float(y) * species_coulomb_entropy(sp.z, t9, rho, ye)
    return float(total)
