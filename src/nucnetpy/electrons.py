"""The electron-positron gas: chemical potential, pressure, energy, entropy.

NucNet Tools takes the properties of the electrons from libstatmech.  This
module is the pure-Python counterpart.  Given the temperature, the density and
the electron fraction ``Ye`` it finds the electron chemical potential and,
from it, everything else the rest of the package needs:

* the pressure, energy and entropy of the electrons and positrons, used by
  :func:`nucnetpy.thermo.entropy_per_nucleon`;
* the electrons' contribution to Debye screening, used by
  :class:`nucnetpy.screening.SkyNetScreening`.

The gas is treated exactly: the electrons may be relativistic and degenerate
to any degree, and the positrons that appear in hot, thin matter are included.
Their chemical potential is the negative of the electrons' (the reaction
``e- + e+ <-> photons`` is in equilibrium), and charge neutrality fixes the
net number of electrons per nucleon to ``Ye``.

Notation.  ``beta = kT / (m_e c^2)`` and ``eta`` is the electron chemical
potential without the rest mass, in units of kT.  The positrons then have
``eta+ = -eta - 2/beta``.  All quantities come from the integrals

    I = integral over x of  (powers of x and (1 + beta x)) / (exp(x - eta) + 1)

where ``x`` is the kinetic energy in units of kT.  They are evaluated by
Gauss-Legendre quadrature in ``t = sqrt(x)``, which removes the square-root
behaviour at ``x = 0``, with the intervals packed around the Fermi edge at
``x = eta``.  The quadrature agrees with SciPy's adaptive integration to about
one part in 10^12 (see the tests).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .constants import AVOGADRO, ELECTRON_MASS_MEV, HBAR_C_MEV_FM, KB_MEV, MEV_TO_ERG

# (m_e c / hbar)^3 in cm^-3 and m_e c^2 in erg.
_COMPTON3 = (ELECTRON_MASS_MEV / (HBAR_C_MEV_FM * 1.0e-13)) ** 3
_ME_ERG = ELECTRON_MASS_MEV * MEV_TO_ERG
# Prefactors of the number density, pressure and energy density (see below).
_N0 = math.sqrt(2.0) / math.pi ** 2 * _COMPTON3
_P0 = 2.0 * math.sqrt(2.0) / (3.0 * math.pi ** 2) * _COMPTON3 * _ME_ERG
_U0 = math.sqrt(2.0) / math.pi ** 2 * _COMPTON3 * _ME_ERG

# Fermi-Dirac occupancy is below exp(-_TAIL) beyond x = eta + _TAIL.
_TAIL = 80.0
_EDGE_OFFSETS = (-60.0, -30.0, -15.0, -8.0, -4.0, -2.0, 0.0, 2.0, 4.0, 8.0, 15.0, 30.0)
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(32)


def _grid(eta: float) -> Tuple[np.ndarray, np.ndarray]:
    """Quadrature points ``t`` and weights covering ``0 <= t^2 <= eta + _TAIL``."""
    t_max = math.sqrt(max(eta, 0.0) + _TAIL)
    edges = {0.0, t_max}
    for d in _EDGE_OFFSETS:
        if 0.0 < eta + d < t_max * t_max:
            edges.add(math.sqrt(eta + d))
    # Doubling intervals below the edge keep each one short compared with its
    # distance from the branch point of sqrt(1 + beta x / 2), so a strongly
    # degenerate, relativistic gas is integrated as accurately as a hot one.
    t = 0.125
    while t < t_max:
        edges.add(t)
        t *= 2.0
    edges = sorted(edges)
    lo = np.array(edges[:-1])[:, None]
    hi = np.array(edges[1:])[:, None]
    half = 0.5 * (hi - lo)
    t = (lo + half * (_NODES + 1.0)).ravel()
    w = (half * _WEIGHTS).ravel()
    return t, w


def _integrals(eta: float, beta: float) -> Tuple[float, np.ndarray]:
    """Return ``(log_scale, I)`` for one species of the gas.

    ``I * exp(log_scale)`` holds, in order, the integrals for the number
    density, its derivative with respect to ``eta``, the pressure and the
    kinetic energy density, each without its prefactor.  The scale keeps the
    numbers representable when the gas is far from degenerate (``eta`` very
    negative), where the occupancy is close to ``exp(eta - x)``.
    """
    t, w = _grid(eta)
    x = t * t
    rel = np.sqrt(1.0 + 0.5 * beta * x)              # sqrt(1 + beta x / 2)
    one_bx = 1.0 + beta * x
    log_scale = min(eta, 0.0)
    # log of the occupancy f and of f (1 - f), divided by exp(log_scale)
    log_f = -np.logaddexp(0.0, x - eta) - log_scale
    log_df = log_f - np.logaddexp(0.0, eta - x)
    f = np.exp(log_f)
    df = np.exp(log_df)
    # dx = 2 t dt, so x^(1/2) dx = 2 t^2 dt and x^(3/2) dx = 2 t^4 dt.
    jac_n = 2.0 * t * t * rel * one_bx
    jac_p = 2.0 * t ** 4 * rel ** 3
    jac_u = 2.0 * t ** 4 * rel * one_bx
    out = np.array([np.dot(w, jac_n * f), np.dot(w, jac_n * df),
                    np.dot(w, jac_p * f), np.dot(w, jac_u * f)])
    return log_scale, out


def _species(eta: float, beta: float) -> Tuple[float, float, float, float]:
    """Number density, d(n)/d(eta), pressure and kinetic energy density (cgs)."""
    log_scale, (i_n, i_dn, i_p, i_u) = _integrals(eta, beta)
    if log_scale < -700.0:
        return 0.0, 0.0, 0.0, 0.0
    scale = math.exp(log_scale)
    b32 = beta ** 1.5
    b52 = beta ** 2.5
    return (float(_N0 * b32 * i_n * scale), float(_N0 * b32 * i_dn * scale),
            float(_P0 * b52 * i_p * scale), float(_U0 * b52 * i_u * scale))


def _net_density(xi: float, beta: float) -> Tuple[float, float]:
    """Return the net electron density ``n- - n+`` and its derivative.

    The argument is ``xi = eta + 1/beta``, the electron chemical potential
    measured from the point where electrons and positrons cancel, so
    ``eta = xi - 1/beta`` and ``eta+ = -xi - 1/beta``.  The derivative is with
    respect to ``xi`` (equal to the one with respect to ``eta``).

    In hot, thin matter the pairs can outnumber the net electrons by many
    orders of magnitude, and subtracting the two densities would lose that
    many digits.  Instead the difference of the two occupancies is formed
    analytically,

        f(eta) - f(eta+) = exp(x - eta) expm1(2 xi) / ((exp(x - eta) + 1)(exp(x - eta+) + 1)),

    and ``xi`` itself is kept as the unknown, so it keeps its full precision
    however small it is.
    """
    eta = xi - 1.0 / beta
    eta_p = -xi - 1.0 / beta
    t, w = _grid(eta)
    x = t * t
    log_scale = min(eta, 0.0)
    jac_n = 2.0 * t * t * np.sqrt(1.0 + 0.5 * beta * x) * (1.0 + beta * x)
    log_f = -np.logaddexp(0.0, x - eta)
    # d(n-)/d(eta): occupancy times (1 - occupancy), as in _integrals.
    log_df = log_f - np.logaddexp(0.0, eta - x) - log_scale
    deriv = float(np.dot(w, jac_n * np.exp(log_df)))
    if xi > 0.0:
        two_xi = 2.0 * xi
        log_expm1 = two_xi + math.log(-math.expm1(-two_xi))
        log_diff = (x - eta) + log_expm1 + log_f - np.logaddexp(0.0, x - eta_p) - log_scale
        net = float(np.dot(w, jac_n * np.exp(log_diff)))
    else:
        net = 0.0
    if log_scale < -700.0:
        net = deriv = 0.0
    else:
        factor = _N0 * beta ** 1.5 * math.exp(log_scale)
        net *= factor
        deriv *= factor
    return net, deriv + _species(eta_p, beta)[1]


@dataclass(frozen=True)
class ElectronGas:
    """The state of the electron-positron gas at one ``(T9, rho, Ye)``.

    Densities are in cm^-3, the pressure and energy densities in erg cm^-3
    and the entropy density in units of k_B per cm^3.  The energies are
    kinetic energies; the rest mass of the positron pairs is in
    ``pair_rest_energy``.
    """
    t9: float
    rho: float
    ye: float
    eta: float                 # electron chemical potential / kT, without rest mass
    beta: float                # kT / (m_e c^2)
    n_electrons: float
    n_positrons: float
    pressure: float
    energy: float
    entropy_density: float
    dn_deta: float             # d(n- - n+)/d(eta) at fixed T, both species

    @property
    def mu_mev(self) -> float:
        """Electron chemical potential in MeV, including the rest mass."""
        return (self.eta * self.beta + 1.0) * ELECTRON_MASS_MEV

    @property
    def pair_rest_energy(self) -> float:
        """Rest energy of the positrons and the electrons paired with them, erg cm^-3."""
        return 2.0 * self.n_positrons * _ME_ERG

    @property
    def entropy_per_nucleon(self) -> float:
        """Entropy of the electrons and positrons per nucleon, in units of k_B."""
        return self.entropy_density / (self.rho * AVOGADRO)

    @property
    def degeneracy_factor(self) -> float:
        """How strongly the electrons screen compared with a classical gas.

        ``(d n / d eta) / (n- - n+)``: one for a non-degenerate electron gas,
        falling towards zero as the electrons become degenerate, and above one
        when positron pairs add to the screening.
        """
        net = self.n_electrons - self.n_positrons
        return self.dn_deta / net if net > 0.0 else float("inf")

    @property
    def screening_term(self) -> float:
        """The electrons' share of the Debye sum, per nucleon.

        Debye screening by the ions is ``sum Z^2 Y``; the electrons and
        positrons add ``kT dn/dmu / (rho N_A)``, which is this value.  It
        equals ``Ye`` for non-degenerate electrons and tends to zero for
        degenerate ones.
        """
        return self.dn_deta / (self.rho * AVOGADRO)


def electron_chemical_potential(t9: float, rho: float, ye: float) -> float:
    """Return ``eta``, the electron chemical potential over kT without rest mass.

    It is the value for which the net electron density ``n- - n+`` equals
    ``rho N_A Ye``.
    """
    xi, beta = _solve_xi(t9, rho, ye)
    return xi - 1.0 / beta


def _solve_xi(t9: float, rho: float, ye: float) -> Tuple[float, float]:
    """Return ``(xi, beta)`` with ``xi = eta + 1/beta`` for the state ``(t9, rho, ye)``.

    The equation ``n- - n+ = rho N_A Ye`` is solved for ``ln(n- - n+)`` by
    Newton's method, kept inside a shrinking bracket, in the variable ``xi``
    (see :func:`_net_density`).
    """
    t9 = float(t9)
    rho = float(rho)
    ye = float(ye)
    if t9 <= 0.0 or rho <= 0.0 or ye <= 0.0:
        raise ValueError("t9, rho and ye must all be positive")
    beta = KB_MEV * t9 * 1.0e9 / ELECTRON_MASS_MEV
    target = rho * AVOGADRO * ye
    log_target = math.log(target)

    # First guess: the larger of the non-degenerate, cold degenerate and
    # pair-dominated estimates.  At xi = 0 the net density is zero.
    n_quantum = 2.0 * _COMPTON3 * (beta / (2.0 * math.pi)) ** 1.5
    if target < n_quantum:
        eta = math.log(target / n_quantum)
    else:
        p_fermi = (3.0 * math.pi ** 2 * target / _COMPTON3) ** (1.0 / 3.0)   # in m_e c
        eta = (math.sqrt(1.0 + p_fermi ** 2) - 1.0) / beta
    xi = eta + 1.0 / beta
    if xi <= 0.0:
        # Pairs dominate: the net density grows linearly from xi = 0.
        _, slope0 = _net_density(0.0, beta)
        xi = target / slope0 if slope0 > 0.0 else 1.0

    lo, hi = 0.0, math.inf
    for _ in range(200):
        net, deriv = _net_density(xi, beta)
        if abs(net - target) <= 1.0e-13 * target:
            return xi, beta
        if not net < target:            # includes an overflow to inf or nan
            hi = xi
            if not math.isfinite(net):
                xi = 0.5 * (lo + hi) if lo > 0.0 else 0.5 * hi
                continue
        else:
            lo = xi
        if net > 0.0 and deriv > 0.0:
            new = xi + (log_target - math.log(net)) * net / deriv
        else:
            new = 2.0 * xi
        if not (lo < new < hi):
            if math.isinf(hi):
                new = max(2.0 * xi, xi + 1.0)
            elif lo > 0.0:
                new = math.sqrt(lo * hi) if hi / lo > 4.0 else 0.5 * (lo + hi)
            else:
                new = 0.1 * hi
        if math.isfinite(hi) and hi - lo <= 1.0e-15 * hi:
            return 0.5 * (lo + hi), beta
        xi = new
    raise RuntimeError(f"electron chemical potential did not converge at t9={t9}, rho={rho}, ye={ye}")


def electron_gas(t9: float, rho: float, ye: float) -> ElectronGas:
    """Return the full :class:`ElectronGas` state at ``(t9, rho, ye)``.

    ``t9`` is the temperature in 10^9 K, ``rho`` the density in g cm^-3 and
    ``ye`` the net number of electrons per nucleon.
    """
    xi, beta = _solve_xi(t9, rho, ye)
    eta = xi - 1.0 / beta
    eta_p = -xi - 1.0 / beta
    n_m, dn_m, p_m, u_m = _species(eta, beta)
    n_p, dn_p, p_p, u_p = _species(eta_p, beta)
    kt_erg = beta * _ME_ERG
    # T s = u + P - mu n for each species, with kinetic u and mu.
    s_m = (u_m + p_m) / kt_erg - eta * n_m
    s_p = (u_p + p_p) / kt_erg - eta_p * n_p
    return ElectronGas(t9=float(t9), rho=float(rho), ye=float(ye), eta=eta, beta=beta,
                       n_electrons=n_m, n_positrons=n_p, pressure=p_m + p_p,
                       energy=u_m + u_p, entropy_density=s_m + s_p,
                       dn_deta=dn_m + dn_p)
