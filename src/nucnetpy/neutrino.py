"""Neutrino reaction-rate helper functions."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import math
import numpy as np
from .mathutils import linear_interpolation

@dataclass
class NeutrinoLuminosity:
    l0: float
    tau: float | None = None

    def luminosity(self, t: float) -> float:
        if self.tau is None or self.tau <= 0:
            return float(self.l0)
        return float(self.l0) * math.exp(-float(t) / float(self.tau))

@dataclass
class NeutrinoQuantity:
    temperature: Sequence[float]
    log10_xsec: Sequence[float]

    def cross_section(self, tnu: float) -> float:
        return 10.0 ** linear_interpolation(self.temperature, self.log10_xsec, tnu)


def geometric_flux_rate(luminosity: float, mean_energy_mev: float, radius_cm: float, cross_section_cm2: float) -> float:
    mev_to_erg = 1.602176634e-6
    number_luminosity = float(luminosity) / max(float(mean_energy_mev) * mev_to_erg, 1e-300)
    flux = number_luminosity / (4.0 * math.pi * max(float(radius_cm), 1e-99)**2)
    return float(flux * float(cross_section_cm2))
