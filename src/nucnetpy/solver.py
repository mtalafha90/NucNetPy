"""One-zone and multi-zone evolution tools.

This layer implements the standard stoichiometric ODE interface and provides
fixed-step, SciPy BDF/Radau/LSODA, optional sparse Jacobians, positivity
projection, screening hooks, and weak-rate source terms.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np

from .core import Network, Zone
from .species import normalize_species_name

ThermoFunction = Callable[[float, Mapping[str, float]], Tuple[float, float]]
ScreeningFunction = Callable[[object, float, float, Optional[float]], float]


def constant_thermo(t9: float = 1.0, rho: float = 1.0) -> ThermoFunction:
    return lambda t, y: (float(t9), float(rho))


def zone_thermo(zone: Zone) -> ThermoFunction:
    return constant_thermo(zone.temperature9(), zone.density())


def time_grid(t0: float, t1: float, n: int, log: bool = False) -> np.ndarray:
    if log:
        return np.geomspace(max(t0, 1e-99), t1, n)
    return np.linspace(t0, t1, n)


@dataclass
class EvolutionResult:
    time: np.ndarray
    species: List[str]
    y: np.ndarray
    success: bool = True
    message: str = ""
    nfev: int = 0
    njev: int = 0

    @property
    def final_abundances(self) -> Dict[str, float]:
        return {s: float(v) for s, v in zip(self.species, self.y[-1])}

    def zone(self) -> Zone:
        return Zone(abundances=self.final_abundances)

    def mass_fraction_history(self, species_map=None) -> np.ndarray:
        factors = []
        from .species import Species
        for name in self.species:
            sp = species_map.get(name) if species_map else None
            if sp is None:
                sp = Species.parse(name)
            factors.append(sp.a)
        return self.y * np.asarray(factors, dtype=float)[None, :]


def _ye_from_vec(species, yvec, species_map=None) -> float:
    from .species import Species
    total = 0.0
    for s, y in zip(species, yvec):
        try:
            sp = species_map.get(s) if species_map else Species.parse(s)
            if sp is None: sp = Species.parse(s)
            total += sp.z * max(float(y), 0.0)
        except Exception:
            continue
    return float(total)


def rhs(network: Network, species: Sequence[str], thermo: ThermoFunction, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None):
    species = [normalize_species_name(s) for s in species]
    idx = {s: i for i, s in enumerate(species)}

    def f(t: float, yvec: np.ndarray) -> np.ndarray:
        yclip = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = {s: float(yclip[i]) for s, i in idx.items()}
        t9, rho = thermo(float(t), abund)
        ye = _ye_from_vec(species, yclip, network.species)
        dy = network.reactions.ydot(abund, t9=t9, rho=rho, screening=screening, ye=ye)
        out = np.zeros_like(yclip, dtype=float)
        for name, val in dy.items():
            if name in idx:
                out[idx[name]] += val
        if weak_rates:
            for wr in weak_rates:
                parent = getattr(wr, 'parent', None); daughter = getattr(wr, 'daughter', None)
                if parent in idx and daughter in idx:
                    rate = wr.rate(t9, rho * max(ye, 1e-30))
                    flow = rate * abund.get(parent, 0.0)
                    out[idx[parent]] -= flow
                    out[idx[daughter]] += flow
        return out
    return f


def jacobian(network: Network, species: Sequence[str], thermo: ThermoFunction, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None, eps: float = 1e-8, sparse: bool = False):
    species = [normalize_species_name(s) for s in species]
    f = rhs(network, species, thermo, screening=screening, weak_rates=weak_rates)

    def j(t: float, yvec: np.ndarray):
        y = np.asarray(yvec, dtype=float)
        base = f(t, y)
        n = len(y)
        J = np.zeros((n, n), dtype=float)
        for k in range(n):
            h = eps * max(abs(y[k]), 1.0)
            yp = y.copy(); yp[k] += h
            J[:, k] = (f(t, yp) - base) / h
        if sparse:
            try:
                from scipy.sparse import csc_matrix
                return csc_matrix(J)
            except Exception:
                return J
        return J
    return j


def jacobian_sparsity(network: Network, species: Sequence[str]):
    species = [normalize_species_name(s) for s in species]
    idx = {s: i for i, s in enumerate(species)}
    mat = np.zeros((len(species), len(species)), dtype=bool)
    for r in network.reactions.reactions:
        affected = [idx[n] for n in r.stoichiometry() if n in idx]
        deps = [idx[p.species] for p in r.reactants if p.species in idx]
        for a in affected:
            for d in deps:
                mat[a, d] = True
    try:
        from scipy.sparse import csc_matrix
        return csc_matrix(mat)
    except Exception:
        return mat


def evolve_zone(network: Network, zone: Zone, times: Sequence[float], thermo: Optional[ThermoFunction] = None, method: str = "bdf", species: Optional[Sequence[str]] = None, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None, rtol: float = 1e-6, atol: float = 1e-30, use_jacobian: bool = True, project_positive: bool = True) -> EvolutionResult:
    ts = np.asarray(times, dtype=float)
    if ts.ndim != 1 or len(ts) < 2:
        raise ValueError("times must be a one-dimensional array with at least two points")
    species = [normalize_species_name(s) for s in (species or network.species_names())]
    y0 = np.array([zone.get_abundance(s) for s in species], dtype=float)
    thermo = thermo or zone_thermo(zone)
    f = rhs(network, species, thermo, screening=screening, weak_rates=weak_rates)
    method_l = method.lower()
    if method_l in {"bdf", "radau", "lsoda", "rk45", "dop853"}:
        try:
            from scipy.integrate import solve_ivp
            kwargs = dict(method=method.upper() if method_l != "lsoda" else "LSODA", rtol=rtol, atol=atol)
            if use_jacobian and method_l in {"bdf", "radau"}:
                kwargs["jac"] = jacobian(network, species, thermo, screening=screening, weak_rates=weak_rates, sparse=True)
                kwargs["jac_sparsity"] = jacobian_sparsity(network, species)
            sol = solve_ivp(f, (ts[0], ts[-1]), y0, t_eval=ts, **kwargs)
            y = sol.y.T
            if project_positive:
                y = np.clip(y, 0.0, np.inf)
            return EvolutionResult(sol.t, species, y, bool(sol.success), str(sol.message), getattr(sol, 'nfev', 0), getattr(sol, 'njev', 0))
        except Exception as exc:
            msg = f"SciPy {method} unavailable/failed ({exc}); used rk4 fallback"
            res = _fixed_step(f, y0, ts, "rk4", project_positive=project_positive)
            return EvolutionResult(ts, species, res, True, msg)
    return EvolutionResult(ts, species, _fixed_step(f, y0, ts, method_l, project_positive=project_positive), True, f"fixed-step {method_l}")


def _fixed_step(f, y0, ts, method, project_positive=True):
    ys = np.zeros((len(ts), len(y0)), dtype=float)
    ys[0] = y0
    for i in range(1, len(ts)):
        h = ts[i] - ts[i-1]
        t = ts[i-1]
        y = ys[i-1]
        if method == "euler":
            yn = y + h * f(t, y)
        elif method in {"implicit_euler", "backward_euler"}:
            yn = _implicit_euler_step(f, t, y, h)
        else:
            k1 = f(t, y)
            k2 = f(t + 0.5*h, y + 0.5*h*k1)
            k3 = f(t + 0.5*h, y + 0.5*h*k2)
            k4 = f(t + h, y + h*k3)
            yn = y + h*(k1 + 2*k2 + 2*k3 + k4)/6.0
        ys[i] = np.clip(yn, 0.0, np.inf) if project_positive else yn
    return ys


def _implicit_euler_step(f, t, y, h, max_iter=12):
    yn = y.copy()
    for _ in range(max_iter):
        g = yn - y - h * f(t + h, yn)
        if np.linalg.norm(g) < 1e-12:
            break
        n = len(y)
        J = np.eye(n)
        eps = 1e-8
        base_f = f(t + h, yn)
        for k in range(n):
            yp = yn.copy(); yp[k] += eps * max(abs(yn[k]), 1.0)
            J[:, k] -= h * (f(t + h, yp) - base_f) / (yp[k] - yn[k])
        try:
            step = np.linalg.solve(J, -g)
        except np.linalg.LinAlgError:
            step = -g
        yn = yn + step
    return yn


def evolve_network_zones(network: Network, times: Sequence[float], thermo: Optional[ThermoFunction] = None, method: str = "bdf", **kwargs) -> List[EvolutionResult]:
    return [evolve_zone(network, z, times, thermo=thermo, method=method, **kwargs) for z in network.zones]
