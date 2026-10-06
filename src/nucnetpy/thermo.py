"""Thermodynamics of the burning matter: entropy, and finding T or rho from it.

:func:`entropy_per_nucleon` adds up the entropy of the three parts of the
matter, each treated in full:

* **ions** - the Sackur-Tetrode entropy of each nuclear species, with its
  real mass, ground-state spin and temperature-dependent partition function;
* **electrons and positrons** - from :mod:`nucnetpy.electrons`, at any
  degeneracy and including pairs;
* **photons** - black-body radiation, ``4 a T^3 / (3 rho N_A k)``;

and, optionally, the Coulomb (plasma) correction of :mod:`nucnetpy.coulomb`.

:func:`density_for_entropy` and :func:`t9_for_entropy` search for the density
or temperature that gives a target entropy at fixed composition, and
:func:`constant_entropy_thermo` uses the second of these to drive a network
calculation along a constant-entropy (adiabatic) expansion.

All entropies are per nucleon in units of Boltzmann's constant.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from .constants import A_RAD, AMU_G, AVOGADRO, KB_CGS, MEV_TO_ERG
from .species import Species, is_massless, normalize_species_name


def ideal_gas_pressure(rho: float, temperature: float, ytot: float) -> float:
    return float(rho) * KB_CGS * float(temperature) * float(ytot) / AMU_G


def radiation_pressure(temperature: float) -> float:
    return A_RAD * float(temperature) ** 4 / 3.0


# ---------------------------------------------------------------------------
# Entropy of the matter
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EntropyPerNucleon:
    """The entropy per nucleon, in units of k_B, split into its parts."""
    ions: float
    electrons: float
    photons: float
    coulomb: float = 0.0

    @property
    def total(self) -> float:
        return self.ions + self.electrons + self.photons + self.coulomb


def _species_map(species) -> Mapping[str, Species]:
    """Accept a network, a mapping of species, or ``None``."""
    if species is None:
        return {}
    return getattr(species, "species", species)


def _lookup(name: str, species_map: Mapping[str, Species]) -> Optional[Species]:
    key = normalize_species_name(name)
    if is_massless(key):
        return None
    sp = species_map.get(key)
    if sp is None:
        try:
            sp = Species.parse(key)
        except Exception:
            return None
    return sp if sp.a > 0 else None


class _IonTable:
    """The parts of the ion entropy that depend on the species, as arrays.

    Each partition-function table is resampled onto the union of all the
    tables' temperature points.  Because :func:`nucnetpy.nse._partition`
    interpolates linearly and holds the end values outside the table, the
    resampling is exact, and one interpolation then serves every species.
    """

    def __init__(self, species: Sequence[Species], include_partition: bool = True):
        from .nse import _log_quantum_abundance

        self.species = list(species)
        self.z = np.array([sp.z for sp in self.species], dtype=float)
        # ln Y_Q without the partition function, at T9 = 1 and rho = 1; the
        # temperature and density enter as 1.5 ln T9 - ln rho.
        self.log_q1 = np.array([_log_quantum_abundance(sp, 1.0, 1.0, include_partition=False)
                                for sp in self.species], dtype=float)
        tables = [sp.partition for sp in self.species if include_partition and sp.partition]
        knots = sorted({float(k) for table in tables for k in table})
        self.grid = np.array(knots, dtype=float)
        self.table = np.ones((len(self.species), len(knots)))
        if include_partition:
            for i, sp in enumerate(self.species):
                if sp.partition:
                    keys = np.array(sorted(sp.partition), dtype=float)
                    vals = np.array([sp.partition[float(k)] for k in keys], dtype=float)
                    self.table[i] = np.interp(self.grid, keys, vals)

    def log_partition(self, t9: float) -> Tuple[np.ndarray, np.ndarray]:
        """Return ``ln G(T)`` and ``d ln G / d ln T`` for every species."""
        n = len(self.species)
        if self.grid.size == 0:
            return np.zeros(n), np.zeros(n)
        if t9 <= self.grid[0] or t9 >= self.grid[-1] or self.grid.size == 1:
            column = 0 if t9 <= self.grid[0] else -1
            return np.log(self.table[:, column]), np.zeros(n)
        k = int(np.searchsorted(self.grid, t9, side="right")) - 1
        x0, x1 = self.grid[k], self.grid[k + 1]
        slope = (self.table[:, k + 1] - self.table[:, k]) / (x1 - x0)
        g = self.table[:, k] + (t9 - x0) * slope
        return np.log(g), t9 * slope / g

    def entropy(self, y: np.ndarray, t9: float, rho: float) -> float:
        """Sackur-Tetrode entropy per nucleon for abundances ``y``."""
        present = y > 0.0
        if not np.any(present):
            return 0.0
        log_g, dlog_g = self.log_partition(t9)
        log_q = self.log_q1 + log_g + 1.5 * math.log(t9) - math.log(rho)
        yp = y[present]
        return float(np.dot(yp, 2.5 + log_q[present] - np.log(yp) + dlog_g[present]))


def photon_entropy_per_nucleon(t9: float, rho: float) -> float:
    """Entropy of black-body radiation per nucleon, ``4 a T^3 / (3 rho N_A k)``."""
    temperature = float(t9) * 1.0e9
    return 4.0 * A_RAD * temperature ** 3 / (3.0 * float(rho) * AVOGADRO * KB_CGS)


def ion_entropy_per_nucleon(abundances: Mapping[str, float], species, t9: float, rho: float,
                            include_partition: bool = True) -> float:
    """Sackur-Tetrode entropy of the nuclei per nucleon, in units of k_B.

    ``sum_i Y_i (5/2 + ln(Y_Q,i / Y_i) + d ln G_i / d ln T)``, where
    ``Y_Q,i = (2J_i+1) G_i(T) n_Q,i / (rho N_A)`` uses the nuclear mass through
    the quantum concentration ``n_Q,i = (A_i m_u kT / 2 pi hbar^2)^(3/2)``.
    The last term is the entropy of the nuclei's excited states.  Species with
    zero abundance and massless particles contribute nothing.
    """
    names, sps = _present_species(abundances, _species_map(species))
    if not sps:
        return 0.0
    y = np.array([float(abundances[n]) for n in names])
    return _IonTable(sps, include_partition).entropy(y, float(t9), float(rho))


def _present_species(abundances: Mapping[str, float], species_map: Mapping[str, Species]):
    names, sps = [], []
    for name, y in abundances.items():
        if float(y) <= 0.0:
            continue
        sp = _lookup(name, species_map)
        if sp is not None:
            names.append(name)
            sps.append(sp)
    return names, sps


def _ye(abundances: Mapping[str, float], species_map: Mapping[str, Species]) -> float:
    total = 0.0
    for name, y in abundances.items():
        sp = _lookup(name, species_map)
        if sp is not None:
            total += sp.z * float(y)
    return total


def entropy_per_nucleon(abundances: Mapping[str, float], species, t9: float, rho: float,
                        ye: Optional[float] = None, include_partition: bool = True,
                        electrons: bool = True, photons: bool = True,
                        coulomb: bool = False) -> EntropyPerNucleon:
    """Return the entropy per nucleon of the matter, split into its parts.

    ``abundances`` maps species names to abundances ``Y``; ``species`` is a
    :class:`~nucnetpy.Network` or a mapping of names to
    :class:`~nucnetpy.Species`, which supplies the spins and partition
    functions (species missing from it are parsed from their names and get a
    weight of one).  ``ye`` defaults to ``sum Z Y`` of the abundances.
    ``t9`` is the temperature in 10^9 K and ``rho`` the density in g cm^-3.

    The flags choose which parts are computed; parts that are switched off
    are reported as zero.  The Coulomb correction is off by default, as in
    the NSE solver.

    Example: pure helium-4 at ``T9 = 5`` and ``rho = 1e8`` has ion entropy
    3.24, electron entropy 1.07 and photon entropy 0.15 per nucleon.
    """
    species_map = _species_map(species)
    t9 = float(t9)
    rho = float(rho)
    if t9 <= 0.0 or rho <= 0.0:
        raise ValueError("t9 and rho must be positive")
    if ye is None:
        ye = _ye(abundances, species_map)
    s_ion = ion_entropy_per_nucleon(abundances, species_map, t9, rho, include_partition)
    s_e = _electron_entropy(t9, rho, ye) if electrons else 0.0
    s_gamma = photon_entropy_per_nucleon(t9, rho) if photons else 0.0
    s_c = 0.0
    if coulomb and ye > 0.0:
        from .coulomb import coulomb_entropy_per_nucleon
        names, sps = _present_species(abundances, species_map)
        s_c = coulomb_entropy_per_nucleon({n: abundances[n] for n in names},
                                          dict(zip(names, sps)), t9, rho, ye)
    return EntropyPerNucleon(ions=s_ion, electrons=s_e, photons=s_gamma, coulomb=s_c)


def _electron_entropy(t9: float, rho: float, ye: float) -> float:
    if ye <= 0.0:
        return 0.0
    from .electrons import electron_gas
    return electron_gas(t9, rho, ye).entropy_per_nucleon


# ---------------------------------------------------------------------------
# Finding the density or temperature for a given entropy
# ---------------------------------------------------------------------------

class _EntropyModel:
    """Entropy as a function of ``(t9, rho)`` at fixed composition."""

    def __init__(self, abundances: Mapping[str, float], species, ye: Optional[float],
                 include_partition: bool, electrons: bool, photons: bool, coulomb: bool):
        species_map = _species_map(species)
        names, sps = _present_species(abundances, species_map)
        self.y = np.array([float(abundances[n]) for n in names])
        self.ions = _IonTable(sps, include_partition)
        self.ye = _ye(abundances, species_map) if ye is None else float(ye)
        self.electrons = electrons and self.ye > 0.0
        self.photons = photons
        self.coulomb = coulomb and self.ye > 0.0
        self.coulomb_args = ({n: abundances[n] for n in names}, dict(zip(names, sps)))

    def __call__(self, t9: float, rho: float) -> float:
        s = self.ions.entropy(self.y, t9, rho)
        if self.electrons:
            s += _electron_entropy(t9, rho, self.ye)
        if self.photons:
            s += photon_entropy_per_nucleon(t9, rho)
        if self.coulomb:
            from .coulomb import coulomb_entropy_per_nucleon
            s += coulomb_entropy_per_nucleon(*self.coulomb_args, t9, rho, self.ye)
        return s


def _solve_log(func: Callable[[float], float], target: float, lo: float, hi: float,
               what: str, guess: Optional[float] = None) -> float:
    """Find ``v`` in ``[lo, hi]`` with ``func(v) = target`` by bisection-safe search in ln v."""
    from scipy.optimize import brentq

    def g(log_v: float) -> float:
        return func(math.exp(log_v)) - target

    a, b = math.log(lo), math.log(hi)
    if guess is not None and lo < guess < hi:
        # Secant steps from the guess: when the state has hardly changed since
        # the last call this needs two or three entropy evaluations, against
        # a dozen or more for a search over the whole range.
        tol = 1e-14 * max(1.0, abs(target))
        x0 = math.log(guess)
        f0 = g(x0)
        if abs(f0) <= tol:
            return guess
        x1 = x0 + (1e-6 if f0 < 0.0 else -1e-6)
        for _ in range(12):
            f1 = g(x1)
            if abs(f1) <= tol:
                return math.exp(x1)
            if f1 == f0 or not (a < x1 < b):
                break
            x0, x1, f0 = x1, x1 - f1 * (x1 - x0) / (f1 - f0), f1
            if abs(x1 - x0) <= 1e-15 * max(1.0, abs(x1)):
                return math.exp(x1)
    fa, fb = g(a), g(b)
    if fa * fb > 0.0:
        raise ValueError(f"entropy {target:g} is outside the range {min(fa, fb) + target:g} "
                         f"to {max(fa, fb) + target:g} reached for {what} between {lo:g} and {hi:g}")
    return math.exp(brentq(g, a, b, xtol=1e-14, rtol=1e-14))


def density_for_entropy(entropy: float, abundances: Mapping[str, float], species, t9: float,
                        ye: Optional[float] = None, bounds: Tuple[float, float] = (1e-12, 1e15),
                        include_partition: bool = True, electrons: bool = True,
                        photons: bool = True, coulomb: bool = False) -> float:
    """Return the density (g cm^-3) at which the matter has the given entropy.

    The composition and the temperature ``t9`` are held fixed, and the
    entropy is :func:`entropy_per_nucleon` with the same options.  The entropy
    falls steadily as the density rises, so the answer is unique; it is found
    by a bracketing search in ``ln rho`` between ``bounds``, to a relative
    precision of about 10^-14.  A ``ValueError`` reports the entropy range
    those bounds can reach if the target is outside it.
    """
    model = _EntropyModel(abundances, species, ye, include_partition, electrons, photons, coulomb)
    return _solve_log(lambda rho: model(float(t9), rho), float(entropy), *bounds, what="rho")


def t9_for_entropy(entropy: float, abundances: Mapping[str, float], species, rho: float,
                   ye: Optional[float] = None, bounds: Tuple[float, float] = (1e-3, 100.0),
                   include_partition: bool = True, electrons: bool = True,
                   photons: bool = True, coulomb: bool = False) -> float:
    """Return the temperature (in 10^9 K) at which the matter has the given entropy.

    The composition and the density ``rho`` are held fixed.  The entropy rises
    with temperature, and the search runs in ``ln T9`` between ``bounds``; see
    :func:`density_for_entropy`.
    """
    model = _EntropyModel(abundances, species, ye, include_partition, electrons, photons, coulomb)
    return _solve_log(lambda t9: model(t9, float(rho)), float(entropy), *bounds, what="t9")


def constant_entropy_thermo(entropy: float, density: Callable[[float], float], species,
                            ye: Optional[float] = None, bounds: Tuple[float, float] = (1e-3, 100.0),
                            include_partition: bool = True, electrons: bool = True,
                            photons: bool = True, coulomb: bool = False):
    """Return a thermodynamic callable for a constant-entropy expansion.

    ``density(t)`` gives the density at time ``t``.  Each time the network
    asks for the conditions, the temperature is set so that the matter, with
    its composition at that moment, has the given ``entropy`` per nucleon.
    Pass the result as the ``thermo`` argument of :func:`nucnetpy.evolve_zone`.

    This is the constant-entropy evolution of the NucNet Tools blog.  The
    entropy is held fixed, so the heat released by the reactions is not fed
    back; that is a good approximation once the burning has slowed down, as in
    the freeze-out from a hot, expanding plasma.  ``ye`` defaults to the
    composition's own ``sum Z Y`` at each call.

    The search starts from the previous temperature, so following a smooth
    expansion costs a few entropy evaluations per call.
    """
    species_map = _species_map(species)
    tables: Dict[Tuple[str, ...], _IonTable] = {}
    last = {"t9": None}

    def thermo(t: float, abundances: Mapping[str, float]) -> Tuple[float, float]:
        rho = float(density(t))
        names, sps = _present_species(abundances, species_map)
        key = tuple(names)
        table = tables.get(key)
        if table is None:
            table = tables[key] = _IonTable(sps, include_partition)
        y = np.array([float(abundances[n]) for n in names])
        ye_now = _ye(abundances, species_map) if ye is None else float(ye)

        def s_of(t9: float) -> float:
            s = table.entropy(y, t9, rho)
            if electrons and ye_now > 0.0:
                s += _electron_entropy(t9, rho, ye_now)
            if photons:
                s += photon_entropy_per_nucleon(t9, rho)
            if coulomb and ye_now > 0.0:
                from .coulomb import coulomb_entropy_per_nucleon
                s += coulomb_entropy_per_nucleon(dict(zip(names, y)), dict(zip(names, sps)),
                                                 t9, rho, ye_now)
            return s

        t9 = _solve_log(s_of, float(entropy), *bounds, what="t9", guess=last["t9"])
        last["t9"] = t9
        return t9, rho

    return thermo


# ---------------------------------------------------------------------------
# Deprecated proxies
# ---------------------------------------------------------------------------

_PROXY_WARNING = ("{name} uses a scale-free proxy, not the physical entropy, and will be "
                  "removed in a future release; use {replacement} instead")


def entropy_ideal_ions(rho: float, temperature: float, ytot: float) -> float:
    """Deprecated: return ``Ytot * (5/2 + ln(T^1.5 / rho))``, a scale-free proxy.

    This is not the physical entropy: it has no nuclear masses, partition
    functions or physical constants, and no electrons or photons.  For pure
    helium-4 at ``T = 5e9 K`` and ``rho = 1e8 g/cm^3`` it gives 4.39 against a
    true ion entropy of 3.24.  Use :func:`entropy_per_nucleon` instead.
    ``temperature`` is in kelvin.
    """
    warnings.warn(_PROXY_WARNING.format(name="entropy_ideal_ions",
                                        replacement="entropy_per_nucleon"),
                  DeprecationWarning, stacklevel=2)
    rho = max(float(rho), 1e-99); temperature = max(float(temperature), 1e-99)
    return float(ytot) * (2.5 + math.log((temperature ** 1.5) / rho))


def density_from_entropy(entropy: float, temperature: float, ytot: float) -> float:
    """Deprecated: invert :func:`entropy_ideal_ions` for the density.

    A physical entropy does not give a physical density here: inverting the
    true entropy of pure helium-4 at ``rho = 1e8 g/cm^3`` returns about
    ``1e10``.  Use :func:`density_for_entropy` instead.
    """
    warnings.warn(_PROXY_WARNING.format(name="density_from_entropy",
                                        replacement="density_for_entropy"),
                  DeprecationWarning, stacklevel=2)
    return float((temperature ** 1.5) / math.exp(float(entropy) / max(ytot, 1e-99) - 2.5))


def temperature_from_entropy(entropy: float, rho: float, ytot: float) -> float:
    """Deprecated: invert :func:`entropy_ideal_ions` for the temperature, in kelvin.

    Use :func:`t9_for_entropy` instead.
    """
    warnings.warn(_PROXY_WARNING.format(name="temperature_from_entropy",
                                        replacement="t9_for_entropy"),
                  DeprecationWarning, stacklevel=2)
    return float((rho * math.exp(float(entropy) / max(ytot, 1e-99) - 2.5)) ** (2.0/3.0))


def sound_speed_ideal_gamma(rho: float, pressure: float, gamma: float = 5.0/3.0) -> float:
    return math.sqrt(max(gamma * pressure / max(rho, 1e-99), 0.0))


def q_energy_erg_per_g(q_mev: float) -> float:
    return float(q_mev) * MEV_TO_ERG * AVOGADRO
