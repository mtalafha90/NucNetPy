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

from .constants import AVOGADRO, KB_MEV, HBAR_C_MEV_FM, AMU_MEV
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
    return math.exp(min(_log_prefactor(sp, t9, rho, include_partition), 700.0))


def _log_prefactor(sp: Species, t9: float, rho: float, include_partition: bool = True) -> float:
    """Natural log of :func:`nse_prefactor`, computed without overflow.

    Working in log space keeps the iron-peak prefactors (which span tens of
    orders of magnitude) representable and lets the solver use a numerically
    stable log-sum-exp for the abundance moments.
    """
    t9 = max(float(t9), 1e-30)
    rho = max(float(rho), 1e-300)
    kt = KB_MEV * t9 * 1.0e9  # kT in MeV
    a = max(int(sp.a), 1)
    g = _partition(sp, t9) if include_partition else 1.0
    # Quantum concentration (cm^-3) for a nucleus of mass A*m_u.  All energies
    # are in MeV and hbar*c is in MeV*cm, so the units are consistent; the
    # nucleon mass enters as its rest energy m_u c^2 = AMU_MEV, not its mass in
    # grams.  The species mass is folded in through the a**1.5 term.
    log_theta = 1.5 * math.log((AMU_MEV * kt) / (2.0 * math.pi * (_HBAR_C_MEV_CM ** 2)))
    log_pref = math.log(max(g, 1e-300)) + 1.5 * math.log(a) + log_theta - math.log(rho * AVOGADRO)
    # More tightly bound nuclei have lower mass excess and are enhanced.  The
    # Z*ME(p)+N*ME(n) piece of the binding energy is linear in Z and N and is
    # absorbed by the proton/neutron chemical potentials, so only -ME_i remains.
    log_pref += -float(sp.mass_excess) / kt
    return float(log_pref)


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


def _log_moments(log_pref: np.ndarray, z: np.ndarray, a: np.ndarray, mu_p: float, mu_n: float) -> Tuple[float, np.ndarray]:
    """Return ``(log(sum A_i Y_i), weights)`` via a stable log-sum-exp.

    ``weights[i] = A_i Y_i / sum_j A_j Y_j`` and sums to one, so charge moments
    can be formed as ``sum (Z_i/A_i) weights_i`` without overflow even when the
    prefactors span tens of orders of magnitude.
    """
    n = a - z
    log_w = log_pref + np.log(a) + z * mu_p + n * mu_n
    m = float(np.max(log_w))
    log_sum_a = m + math.log(float(np.sum(np.exp(log_w - m))))
    weights = np.exp(log_w - log_sum_a)
    return log_sum_a, weights


def solve_nse(network: Network, t9: float, rho: float, ye: float, species: Optional[Sequence[str]] = None, include_partition: bool = True, tol: float = 1e-8, max_iter: int = 200, nse_correction=None) -> NSEResult:
    """Solve an NSE composition for a network.

    The constraints are ``sum(A_i Y_i)=1`` and ``sum(Z_i Y_i)=Ye``.  The solve is
    performed on the well-scaled residuals ``log(sum A_i Y_i)`` and
    ``sum(Z_i Y_i)/sum(A_i Y_i) - Ye`` using a numerically stable log-sum-exp,
    which keeps the iron-peak prefactors representable.  A Levenberg--Marquardt
    least-squares solve is used when SciPy is available; otherwise an internal
    damped LM iteration is used.  Several starting points are tried so the solve
    is robust across temperature/density regimes.

    ``nse_correction`` is an optional callable ``(species, t9, rho, ye) ->
    f_corr`` added to each species' NSE exponent, the libnucnet NSE correction
    factor hook.  Pass :func:`nucnetpy.coulomb.nse_correction` for the Bravo &
    Garcia-Senz Coulomb correction.
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
    log_pref = np.array([_log_prefactor(sp, t9, rho, include_partition=include_partition) for sp in sps], dtype=float)
    if nse_correction is not None:
        log_pref = log_pref + np.array([float(nse_correction(sp, t9, rho, float(ye))) for sp in sps], dtype=float)
    z = np.array([sp.z for sp in sps], dtype=float)
    a = np.array([sp.a for sp in sps], dtype=float)
    ye = float(ye)

    def residual(mus):
        log_sum_a, weights = _log_moments(log_pref, z, a, float(mus[0]), float(mus[1]))
        charge = float(np.sum((z / a) * weights))
        return np.array([log_sum_a, charge - ye], dtype=float)

    guesses = [(0.0, 0.0), (-1.0, -1.0), (-5.0, -5.0), (-10.0, -10.0), (-20.0, -20.0)]
    best_mu = np.array(guesses[0], dtype=float)
    best_norm = float("inf")
    message = ""
    for guess in guesses:
        mu, norm, msg = _solve_mu(residual, np.array(guess, dtype=float), tol, max_iter)
        if norm < best_norm:
            best_mu, best_norm, message = mu, norm, msg
        if best_norm < tol:
            break

    mu_p, mu_n = float(best_mu[0]), float(best_mu[1])
    abund = {sp.name: float(np.exp(min(lp + sp.z * mu_p + sp.n * mu_n, 700.0)))
             for sp, lp in zip(sps, log_pref)}
    success = best_norm < max(tol, 1e-6)
    return NSEResult(t9, rho, ye, mu_p, mu_n, abund, success, message)


def _solve_mu(residual, guess: np.ndarray, tol: float, max_iter: int) -> Tuple[np.ndarray, float, str]:
    """Solve ``residual(mu)=0`` returning ``(mu, residual_norm, message)``."""
    try:
        from scipy.optimize import least_squares
        sol = least_squares(residual, guess, method="lm", xtol=1e-14, ftol=1e-14, max_nfev=max_iter * 10)
        return np.asarray(sol.x, dtype=float), float(np.linalg.norm(sol.fun)), str(sol.message)
    except Exception as exc:
        return _levenberg_marquardt(residual, guess, tol, max_iter) + (f"internal LM ({exc})",)


def _levenberg_marquardt(residual, mu: np.ndarray, tol: float, max_iter: int) -> Tuple[np.ndarray, float]:
    mu = np.clip(mu.astype(float), *_MU_BOUNDS)
    lam = 1e-3
    r = residual(mu)
    norm = float(np.linalg.norm(r))
    eps = 1e-6
    for _ in range(max_iter):
        if norm < tol:
            break
        j = np.column_stack([(residual(mu + [eps, 0]) - r) / eps, (residual(mu + [0, eps]) - r) / eps])
        jtj = j.T @ j
        grad = j.T @ r
        improved = False
        for _ in range(30):
            try:
                step = np.linalg.solve(jtj + lam * np.eye(2), -grad)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            trial = np.clip(mu + step, *_MU_BOUNDS)
            trial_norm = float(np.linalg.norm(residual(trial)))
            if trial_norm < norm:
                mu, r, norm = trial, residual(trial), trial_norm
                lam = max(lam / 3.0, 1e-12)
                improved = True
                break
            lam *= 3.0
            if lam > 1e12:
                break
        if not improved:
            break
    return mu, norm


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
