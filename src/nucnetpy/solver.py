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
from .species import is_massless, normalize_species_name

ThermoFunction = Callable[[float, Mapping[str, float]], Tuple[float, float]]
ScreeningFunction = Callable[[object, float, float, Optional[float]], float]

#: Method names accepted by :func:`evolve_zone`, mapped to the exact spelling
#: that ``scipy.integrate.solve_ivp`` requires.  SciPy matches these names
#: case-sensitively, so the mapping cannot be replaced by ``str.upper``.
_SCIPY_METHODS = {
    "bdf": "BDF",
    "radau": "Radau",
    "lsoda": "LSODA",
    "rk45": "RK45",
    "dop853": "DOP853",
}


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


def analytic_jacobian(network: Network, species: Sequence[str], thermo: ThermoFunction, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None, sparse: bool = True):
    """Return an analytic Jacobian of the stoichiometric right-hand side.

    Differentiating the flow of Eq. (flow),

        F_r = (lambda_r rho^(n_r-1) / prod_i m_ir!) prod_i Y_i^{m_ir},

    with respect to a reactant abundance gives

        dF_r/dY_j = m_jr Y_j^{m_jr-1} (lambda_r rho^(n_r-1) / prod_i m_ir!)
                    prod_{i != j} Y_i^{m_ir},

    so that ``J_ij = sum_r nu_ir dF_r/dY_j``.  The partial product is formed
    explicitly rather than as ``F_r m_jr / Y_j`` so that the derivative stays
    correct when ``Y_j`` vanishes, which is the usual state of most species at
    the start of a calculation.

    The whole matrix costs one sweep over the reaction list, whereas the
    finite-difference Jacobian of :func:`jacobian` costs ``N+1`` evaluations of
    the full right-hand side.  The saving therefore grows linearly with the
    number of evolved species.

    ``lambda_r`` is treated as independent of composition.  With a screening
    function or weak-rate tables that depend on ``Ye`` the result is an
    approximate Jacobian: it changes the Newton convergence rate of an implicit
    solver but not the solution it converges to.
    """
    species = [normalize_species_name(s) for s in species]
    idx = {s: i for i, s in enumerate(species)}
    n = len(species)

    # Precompute the per-reaction index/multiplicity structure once.
    compiled = []
    for r in network.reactions.reactions:
        nuclear = r.nuclear_reactants
        reactants = [(idx[p.species], p.count) for p in nuclear if p.species in idx]
        if len(reactants) != len(nuclear):
            continue  # reaction touches species outside the evolved set
        stoich = [(idx[name], nu) for name, nu in r.stoichiometry().items() if name in idx]
        if not stoich:
            continue
        compiled.append((r, reactants, stoich, r.reactant_order, max(r.statistical_factor(), 1)))

    def j(t: float, yvec: np.ndarray):
        y = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = {s: float(y[i]) for s, i in idx.items()}
        t9, rho = thermo(float(t), abund)
        ye = _ye_from_vec(species, y, network.species)
        J = np.zeros((n, n), dtype=float)
        for r, reactants, stoich, order, stat in compiled:
            lam = r.rate(t9, rho=rho, ye=ye, screening=screening)
            if lam == 0.0:
                continue
            pref = lam * (float(rho) ** max(order - 1, 0)) / stat
            for pos, (jx, mj) in enumerate(reactants):
                # d/dY_j of the reactant product, with the j-th factor replaced
                # by m_j Y_j^{m_j - 1}.
                term = pref * mj * (y[jx] ** (mj - 1))
                for other, (kx, mk) in enumerate(reactants):
                    if other != pos:
                        term *= y[kx] ** mk
                if term == 0.0:
                    continue
                for ix, nu in stoich:
                    J[ix, jx] += nu * term
        if weak_rates:
            for wr in weak_rates:
                parent = getattr(wr, 'parent', None)
                daughter = getattr(wr, 'daughter', None)
                if parent in idx and daughter in idx:
                    rate = wr.rate(t9, rho * max(ye, 1e-30))
                    J[idx[parent], idx[parent]] -= rate
                    J[idx[daughter], idx[parent]] += rate
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


def evolve_zone(network: Network, zone: Zone, times: Sequence[float], thermo: Optional[ThermoFunction] = None, method: str = "bdf", species: Optional[Sequence[str]] = None, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None, rtol: float = 1e-6, atol: float = 1e-30, use_jacobian: bool = True, project_positive: bool = True, jac_mode: str = "analytic") -> EvolutionResult:
    """Evolve one zone over ``times``.

    ``jac_mode`` selects how the Jacobian is supplied to the implicit SciPy
    solvers.  ``"analytic"`` (the default) uses :func:`analytic_jacobian`, which
    differentiates the stoichiometric flows in closed form and costs one sweep
    over the reaction list.  ``"numerical"`` uses the finite-difference
    :func:`jacobian`, costing ``N+1`` right-hand-side evaluations.  ``"sparsity"``
    supplies no Jacobian and instead passes the stoichiometric sparsity pattern
    so that SciPy estimates the matrix by grouped finite differences.

    Note that SciPy applies ``jac_sparsity`` only when no Jacobian callable is
    given, so the pattern and an explicit Jacobian are alternatives rather than
    complements.
    """
    ts = np.asarray(times, dtype=float)
    if ts.ndim != 1 or len(ts) < 2:
        raise ValueError("times must be a one-dimensional array with at least two points")
    # Photons and leptons are part of the reaction records but not of the
    # abundance vector; evolving them would integrate a meaningless quantity.
    species = [normalize_species_name(s) for s in (species or network.species_names())]
    species = [s for s in species if not is_massless(s)]
    y0 = np.array([zone.get_abundance(s) for s in species], dtype=float)
    thermo = thermo or zone_thermo(zone)
    f = rhs(network, species, thermo, screening=screening, weak_rates=weak_rates)
    method_l = method.lower()
    if method_l in _SCIPY_METHODS:
        try:
            from scipy.integrate import solve_ivp
            kwargs = dict(method=_SCIPY_METHODS[method_l], rtol=rtol, atol=atol)
            if use_jacobian and method_l in {"bdf", "radau"}:
                mode = str(jac_mode).lower()
                if mode == "sparsity":
                    # No callable: SciPy groups columns of the pattern and
                    # estimates the Jacobian by finite differences itself.
                    kwargs["jac_sparsity"] = jacobian_sparsity(network, species)
                elif mode == "numerical":
                    kwargs["jac"] = jacobian(network, species, thermo, screening=screening, weak_rates=weak_rates, sparse=True)
                else:
                    kwargs["jac"] = analytic_jacobian(network, species, thermo, screening=screening, weak_rates=weak_rates, sparse=True)
            sol = solve_ivp(f, (ts[0], ts[-1]), y0, t_eval=ts, **kwargs)
        except Exception as exc:
            # SciPy itself is unavailable or rejected the call.  A fixed-step
            # RK4 pass keeps such environments usable, but it is an explicit
            # method on a stiff system, so it is only reported as successful
            # when it actually produced a finite trajectory.
            res = _fixed_step(f, y0, ts, "rk4", project_positive=project_positive)
            ok = bool(np.all(np.isfinite(res)))
            msg = f"SciPy {method} unavailable ({exc}); used rk4 fallback"
            if not ok:
                msg += "; rk4 fallback diverged"
            return EvolutionResult(ts, species, res, ok, msg)
        # A stiff solve that fails may return before emitting any output, in
        # which case sol.y is an empty list rather than an array.  Report the
        # failure and the solver's own diagnosis instead of substituting an
        # explicit method that cannot integrate a stiff network either.
        y = np.asarray(sol.y, dtype=float)
        if y.size == 0:
            return EvolutionResult(np.asarray(sol.t, dtype=float).reshape(0), species,
                                   np.zeros((0, len(species))), False, str(sol.message),
                                   getattr(sol, 'nfev', 0), getattr(sol, 'njev', 0))
        y = y.T
        success, message = bool(sol.success), str(sol.message)
        if project_positive:
            y, success, message = _project_positive(network, species, y, success, message)
        return EvolutionResult(sol.t, species, y, success, message, getattr(sol, 'nfev', 0), getattr(sol, 'njev', 0))
    res = _fixed_step(f, y0, ts, method_l, project_positive=project_positive)
    ok = bool(np.all(np.isfinite(res)))
    msg = f"fixed-step {method_l}" + ("" if ok else " diverged (non-finite abundances)")
    return EvolutionResult(ts, species, res, ok, msg)


#: Fraction of the initial baryon number that positivity projection may create
#: before the trajectory is reported as untrustworthy.  Clipping round-off noise
#: is harmless; clipping a genuinely negative abundance is not, because the
#: baryon number it removes is invented rather than transferred.
_PROJECTION_TOLERANCE = 1.0e-6


def _project_positive(network: Network, species: Sequence[str], y: np.ndarray,
                      success: bool, message: str):
    """Clip negative abundances, reporting how much baryon number that creates.

    The right-hand side clips its input to non-negative values, so once a
    component goes negative the derivative no longer depends on it and the
    corresponding Jacobian column vanishes.  An implicit solver can then take a
    component far negative while still satisfying its own convergence test.
    Projecting that away silently manufactures mass, so the amount created is
    measured and a materially non-conservative trajectory is reported as a
    failure rather than returned as a result.
    """
    from .species import Species
    a = np.zeros(len(species), dtype=float)
    for i, name in enumerate(species):
        sp = network.species.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        a[i] = sp.a
    clipped = np.clip(y, 0.0, np.inf)
    created = (clipped - y) @ a
    scale = max(float(np.abs(y[0] @ a)), 1e-30)
    worst = float(np.max(np.abs(created))) / scale
    if worst > _PROJECTION_TOLERANCE:
        success = False
        message = (f"{message}; positivity projection created {worst:.3e} of the "
                   f"initial baryon number, so the trajectory is not trustworthy")
    return clipped, success, message


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
