"""Thermodynamic helper functions used by NucNet Tools examples.

These are practical approximations for pure-Python workflows.  For exact
publication comparisons, validate constants and EOS choices against the C++
libnucnet/libnuceq configuration used in the original calculation.
"""
from __future__ import annotations

import math
from typing import Mapping, Optional
from .constants import KB_CGS, AMU_G, MEV_TO_ERG
from .species import Species


def ideal_gas_pressure(rho: float, temperature: float, ytot: float) -> float:
    return float(rho) * KB_CGS * float(temperature) * float(ytot) / AMU_G


def radiation_pressure(temperature: float) -> float:
    a_rad = 7.5657e-15
    return a_rad * float(temperature) ** 4 / 3.0


def entropy_ideal_ions(rho: float, temperature: float, ytot: float) -> float:
    """Return ``Ytot * (5/2 + ln(T^1.5 / rho))``, a scale-free entropy proxy.

    This is not the physical ion entropy.  The Sackur-Tetrode entropy per
    nucleon is ``sum_i Y_i (5/2 + ln(Y_Q,i / Y_i))``, and the quantum abundance
    ``Y_Q,i`` carries the nuclear mass, the partition function and physical
    constants that are absent here.  The proxy therefore has the right
    dependence on ``T`` and ``rho`` at fixed composition, so differences in it
    are meaningful, but its absolute value is not: for pure helium-4 at
    ``T = 5e9 K`` and ``rho = 1e8 g/cm^3`` it gives 4.39 against the true 3.24.
    ``temperature`` is in kelvin.
    """
    rho = max(float(rho), 1e-99); temperature = max(float(temperature), 1e-99)
    return float(ytot) * (2.5 + math.log((temperature ** 1.5) / rho))


def density_from_entropy(entropy: float, temperature: float, ytot: float) -> float:
    """Invert :func:`entropy_ideal_ions` for the density.

    Because the proxy omits the quantum abundance, a physical entropy does not
    give a physical density: inverting the true entropy of pure helium-4 at
    ``rho = 1e8 g/cm^3`` returns about ``1e10``.  Use it only with values from
    :func:`entropy_ideal_ions` itself, for example to find the density change
    that keeps the proxy constant as the temperature changes.
    """
    return float((temperature ** 1.5) / math.exp(float(entropy) / max(ytot, 1e-99) - 2.5))


def temperature_from_entropy(entropy: float, rho: float, ytot: float) -> float:
    """Invert :func:`entropy_ideal_ions` for the temperature, in kelvin.

    The same caveat as :func:`density_from_entropy` applies: only values from
    :func:`entropy_ideal_ions` give meaningful results.
    """
    return float((rho * math.exp(float(entropy) / max(ytot, 1e-99) - 2.5)) ** (2.0/3.0))


def sound_speed_ideal_gamma(rho: float, pressure: float, gamma: float = 5.0/3.0) -> float:
    return math.sqrt(max(gamma * pressure / max(rho, 1e-99), 0.0))


def q_energy_erg_per_g(q_mev: float) -> float:
    from .constants import AVOGADRO
    return float(q_mev) * MEV_TO_ERG * AVOGADRO
