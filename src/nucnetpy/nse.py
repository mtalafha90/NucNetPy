"""Nuclear statistical equilibrium (NSE) and equilibrium helpers.

This module is a pure-Python counterpart to the most commonly used libnuceq
workflows.  It solves for chemical potentials that reproduce a requested
``rho``, ``T9`` and ``Ye`` and then returns NSE abundances for the species that
are present in a :class:`nucnetpy.Network`.

The formulas are intentionally transparent and unit-aware.  They are suitable
for replacement workflows and regression testing, but exact equality with a
specific compiled libnuceq build should be checked with the same nuclear masses,
partition functions, Coulomb corrections, and numerical tolerances.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple
import math
import numpy as np

from .constants import AVOGADRO, KB_MEV, HBAR_C_MEV_FM, AMU_G
from .core import Network, Zone
from .species import Species, normalize_species_name

# hbar*c in MeV cm
_HBAR_C_MEV_CM = HBAR_C_MEV_FM * 1.0e-13
_MU_BOUNDS = (-500.0, 500.0)


def _partition(sp: Species, t9: float) -> float:
    if not sp.partition:
        return 1.0
    keys = np.array(sorted(sp.partition), dtype=float)
    vals = np.array([sp.partition[float(k)] for k in keys], dtype=float)
    if len(keys) == 1:
        return float(vals[0])
    return float(np.interp(float(t9), keys, vals))


def nse_prefactor(sp: Species, t9: float, rho: float, include_partition: bool = True) -> float:
    """Return the NSE abundance prefactor excluding proton/neutron potentials.

    The abundance is represented as ``Y_i = prefactor_i * exp(Z_i*mu_p/kT +
    N_i*mu_n/kT)`` where ``mu_p`` and ``mu_n`` are solved dimensionless chemical
    potential parameters in the internal solver.  The mass excess term is
    included with the sign convention common in NSE abundance formulas.
    """
    t9 = max(float(t9), 1e-30)
    rho = max(float(rho), 1e-300)
    kt = KB_MEV * t9 * 1.0e9
    a = max(int(sp.a), 1)
    g = _partition(sp, t9) if include_partition else 1.0
    # Translational factor in cgs.  This compact form is used for numerical
    # replacement workflows; use identical constants/masses for strict regression.
    theta = ((AMU_G * kt) / (2.0 * math.pi * (_HBAR_C_MEV_CM ** 2))) ** 1.5
    pref = g * (a ** 1.5) / (rho * AVOGADRO) * theta
    # More tightly bound nuclei have lower mass excess and are enhanced.
    pref *= math.exp(-float(sp.mass_excess) / kt) if sp.mass_excess else 1.0
    return max(float(pref), 1e-300)


@dataclass
class NSEResult:
    t9: float
    rho: float
    ye: float
    mu_p: float
    mu_n: float
    abundances: Dict[str, float]
    success: bool
    message: str = ""

    def zone(self) -> Zone:
        return Zone(abundances=dict(self.abundances), properties={"t9": str(self.t9), "rho": str(self.rho), "ye": str(self.ye)})

    @property
    def xsum(self) -> float:
        return float(sum(Species.parse(k).a * v for k, v in self.abundances.items()))

    @property
    def computed_ye(self) -> float:
        return float(sum(Species.parse(k).z * v for k, v in self.abundances.items()))


def _safe_exp(x: float) -> float:
    if x > 700:
        return math.exp(700)
    if x < -745:
        return 0.0
    return math.exp(x)


def _moments(species: Sequence[Species], pref: np.ndarray, mu_p: float, mu_n: float) -> Tuple[float, float, Dict[str, float]]:
    xsum = 0.0
    ye = 0.0
    abund: Dict[str, float] = {}
    for sp, pr in zip(species, pref):
        expo = sp.z * mu_p + sp.n * mu_n
        y = float(pr) * _safe_exp(expo)
        abund[sp.name] = y
        xsum += sp.a * y
        ye += sp.z * y
    return float(xsum), float(ye), abund


def solve_nse(network: Network, t9: float, rho: float, ye: float, species: Optional[Sequence[str]] = None, include_partition: bool = True, tol: float = 1e-10, max_iter: int = 100) -> NSEResult:
    """Solve an NSE composition for a network.

    The constraints are ``sum(A_i Y_i)=1`` and ``sum(Z_i Y_i)=Ye``.  A robust
    SciPy root solve is used when available; otherwise a damped Newton method is
    used.
    """
    names = [normalize_species_name(s) for s in (species or network.species_names())]
    sps = []
    for name in names:
        sp = network.species.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        if sp.a > 0:
            sps.append(sp)
    if not sps:
        raise ValueError("No valid species available for NSE solve")
    pref = np.array([nse_prefactor(sp, t9, rho, include_partition=include_partition) for sp in sps], dtype=float)

    def residual(mus):
        xs, ys, _ = _moments(sps, pref, float(mus[0]), float(mus[1]))
        return np.array([xs - 1.0, ys - float(ye)], dtype=float)

    guess = np.array([0.0, 0.0], dtype=float)
    try:
        from scipy.optimize import root
        sol = root(residual, guess, method="hybr", tol=tol)
        mu = np.asarray(sol.x, dtype=float)
        xs, ys, abund = _moments(sps, pref, mu[0], mu[1])
        return NSEResult(t9, rho, ye, float(mu[0]), float(mu[1]), abund, bool(sol.success), str(sol.message))
    except Exception as exc:
        mu = guess.copy()
        msg = f"used internal damped Newton fallback ({exc})"
        ok = False
        for _ in range(max_iter):
            r = residual(mu)
            if float(np.linalg.norm(r, ord=2)) < tol:
                ok = True
                break
            eps = 1e-5
            j = np.column_stack([(residual(mu + [eps, 0]) - r) / eps, (residual(mu + [0, eps]) - r) / eps])
            try:
                step = np.linalg.solve(j, -r)
            except np.linalg.LinAlgError:
                step = -0.1 * r
            damp = 1.0
            old = np.linalg.norm(r)
            while damp > 1e-4:
                trial = np.clip(mu + damp * step, *_MU_BOUNDS)
                if np.linalg.norm(residual(trial)) < old:
                    mu = trial
                    break
                damp *= 0.5
        xs, ys, abund = _moments(sps, pref, mu[0], mu[1])
        return NSEResult(t9, rho, ye, float(mu[0]), float(mu[1]), abund, ok, msg)


def equilibrium_ratio(reaction, network: Network, t9: float, rho: float, ye: float) -> float:
    """Return product/reactant NSE abundance ratio for a reaction."""
    nse = solve_nse(network, t9=t9, rho=rho, ye=ye)
    num = 1.0
    den = 1.0
    for p in reaction.products:
        num *= max(nse.abundances.get(p.species, 0.0), 1e-300) ** p.count
    for p in reaction.reactants:
        den *= max(nse.abundances.get(p.species, 0.0), 1e-300) ** p.count
    return float(num / max(den, 1e-300))
