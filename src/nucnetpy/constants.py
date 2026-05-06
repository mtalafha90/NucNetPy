"""Physical constants used by nucnetpy.

The values are SI/cgs mixed in the conventional astrophysical style used by
network codes.  They are intentionally centralized so advanced users can swap
or audit constants when comparing with libnucnet/NucNet Tools builds.
"""
from __future__ import annotations

AVOGADRO = 6.02214076e23  # mol^-1
KB_MEV = 8.617333262145e-11  # MeV K^-1
KB_CGS = 1.380649e-16  # erg K^-1
MEV_TO_ERG = 1.602176634e-6
AMU_G = 1.66053906660e-24
CLIGHT = 2.99792458e10  # cm s^-1
HBAR = 1.054571817e-27  # erg s
NA = AVOGADRO

HBAR_C_MEV_FM = 197.3269804  # MeV fm
E2_MEV_FM = 1.43996448  # e^2 in MeV fm
