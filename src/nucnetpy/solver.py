"""One-zone and multi-zone evolution tools.

This layer implements the standard stoichiometric ODE interface and provides
fixed-step, SciPy BDF/Radau/LSODA and adaptive implicit Euler integrators,
sparse Jacobians, positivity projection, screening hooks, and weak-rate source
terms.  The right-hand side and the analytic Jacobian evaluate all the
reactions at once, as arrays.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple
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
#: Implicit Euler with its own step control, the scheme of NucNet Tools; see
#: :func:`evolve_zone`.
_ADAPTIVE_METHODS = ("adaptive_implicit_euler",)
#: Fixed-step integrators implemented here.  They take one step per interval
#: of the time grid, so the grid itself sets the accuracy.
_FIXED_STEP_METHODS = ("rk4", "euler", "implicit_euler", "backward_euler")
#: Every method name accepted by :func:`evolve_zone` (case-insensitive).
METHODS = tuple(_SCIPY_METHODS) + _ADAPTIVE_METHODS + _FIXED_STEP_METHODS
_JAC_MODES = ("analytic", "numerical", "sparsity")


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
    #: Accepted and rejected steps; counted by ``adaptive_implicit_euler`` only.
    steps: int = 0
    rejected_steps: int = 0

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


def _exp_of_fit(exponent: float) -> float:
    """``exp`` as :meth:`RateFit.rate` takes it: infinite above 700, zero below -745.

    ``math.exp`` rather than ``numpy.exp``, which can differ in the last bit.
    """
    if exponent > 700:
        return math.inf
    if exponent < -745:
        return 0.0
    return math.exp(exponent)


class _NetworkArrays:
    """The reactions among a set of evolved species, held as arrays.

    It is built once per right-hand side or Jacobian.  Each evaluation is then
    a few NumPy operations and one sparse matrix product, instead of a Python
    loop over the reactions.  The arithmetic follows :meth:`Reaction.rate` and
    :meth:`Reaction.flux` term by term and in the same order, so the two agree
    to rounding error.

    The parts of a rate that depend on the temperature alone (the ReacLib fits
    with their partition-function factors, and any tabulated, constant or
    computed rate) are evaluated once for each temperature.  The electron
    density of an electron capture and any screening factor are applied at
    every evaluation.  The reactions and their rate parameters are read when
    the arrays are built: change the network afterwards and build them again.
    """

    def __init__(self, network: Network, species: Sequence[str]):
        from scipy.sparse import csr_matrix
        from .species import Species

        self.species = list(species)
        self.index = {s: i for i, s in enumerate(self.species)}
        n = self.n = len(self.species)
        self.reactions = list(network.reactions.reactions)
        m = self.m = len(self.reactions)
        width = max([len(r.nuclear_reactants) for r in self.reactions] + [1])
        # One slot per nuclear reactant.  Unused slots point at entry n of the
        # extended abundance vector, which holds 1, with a power of 0.
        self.slot_index = np.full((m, width), n, dtype=np.int64)
        self.slot_power = np.zeros((m, width), dtype=np.int64)
        # A reaction with a reactant outside the evolved set has a zero flow,
        # because that reactant's abundance is zero.
        self.active = np.ones(m, dtype=bool)
        self.density_power = np.zeros(m, dtype=np.int64)
        self.statistical = np.ones(m)
        self.electrons = np.zeros(m, dtype=np.int64)
        rows, cols, values, coefficients, owners = [], [], [], [], []
        self.other_rates, self.partition = [], []
        for k, r in enumerate(self.reactions):
            for j, p in enumerate(r.nuclear_reactants):
                if p.species in self.index:
                    self.slot_index[k, j] = self.index[p.species]
                    self.slot_power[k, j] = p.count
                else:
                    self.active[k] = False
            self.density_power[k] = max(r.reactant_order - 1, 0)
            self.statistical[k] = max(r.statistical_factor(), 1)
            self.electrons[k] = r.electron_count
            for name, nu in r.stoichiometry().items():
                if nu and name in self.index and not is_massless(name):
                    rows.append(self.index[name])
                    cols.append(k)
                    values.append(float(nu))
            for fit in r.rate_fits:
                a = list(fit.coefficients)
                if len(a) != 7:
                    raise ValueError("ReacLib RateFit requires exactly 7 coefficients")
                coefficients.append([float(c) for c in a])
                owners.append(k)
            if (r.constant_rate is not None or r.tabular_rate is not None
                    or r.rate_function is not None):
                self.other_rates.append(k)
            if r.partition_ratio is not None:
                self.partition.append(k)
        #: ``dY_i/dt = sum_r stoichiometry[i, r] * flow_r``.
        self.stoichiometry = csr_matrix((values, (rows, cols)), shape=(n, m))
        self.coefficients = np.asarray(coefficients, dtype=float).reshape(-1, 7)
        self.owners = np.asarray(owners, dtype=np.int64)
        self.captures = np.nonzero(self.electrons)[0]
        # Charges for Ye = sum Z Y: from the nuclide record, else from the
        # name; a name that is not a nuclide contributes nothing.
        charges = np.zeros(n)
        for i, s in enumerate(self.species):
            try:
                sp = network.species.get(s) if network.species else None
                charges[i] = (sp if sp is not None else Species.parse(s)).z
            except Exception:
                charges[i] = 0.0
        self.charges = charges
        # For each slot, the reactions whose reactant there appears once, and
        # those where it appears more than once (p + p, 3 alpha, ...).  The
        # latter take Python's power, as Reaction.flux does; numpy.power can
        # differ from it in the last bit.
        self._single = [np.nonzero(self.slot_power[:, j] == 1)[0] for j in range(width)]
        self._repeated = [np.nonzero(self.slot_power[:, j] > 1)[0] for j in range(width)]
        self._t9 = None
        self._other = self._fitted = None
        self._abundance = np.ones(n + 1)

    def electron_fraction(self, y: np.ndarray) -> float:
        """``Ye = sum_i Z_i Y_i``, summed one term at a time in species order."""
        if self.n == 0:
            return 0.0
        return float(np.cumsum(self.charges * y)[-1])

    def _temperature_parts(self, t9: float):
        """Rate parts that depend on T9 alone, cached for the last T9."""
        if self._t9 is not None and t9 == self._t9:
            return self._other, self._fitted
        t9c = max(float(t9), 1e-30)
        t13 = t9c ** (1.0 / 3.0)
        a = self.coefficients
        # The sum of RateFit.rate, in its order.
        exponent = a[:, 0] + a[:, 1] / t9c
        exponent += a[:, 2] / t13
        exponent += a[:, 3] * t13
        exponent += a[:, 4] * t9c
        exponent += a[:, 5] * (t9c ** (5.0 / 3.0))
        exponent += a[:, 6] * math.log(t9c)
        values = np.array([_exp_of_fit(e) for e in exponent.tolist()])
        fitted = np.bincount(self.owners, weights=values, minlength=self.m)
        for k in self.partition:
            if fitted[k]:
                fitted[k] *= self.reactions[k].partition_ratio(t9)
        other = np.zeros(self.m)
        for k in self.other_rates:
            r = self.reactions[k]
            total = 0.0
            if r.constant_rate is not None:
                total += float(r.constant_rate)
            if r.tabular_rate is not None:
                total += r.tabular_rate.rate(t9)
            if r.rate_function is not None:
                total += float(r.rate_function(t9))
            other[k] = total
        self._t9, self._other, self._fitted = t9, other, fitted
        return other, fitted

    def prefactors(self, t9: float, rho: float, ye: float, screening=None) -> np.ndarray:
        """``lambda_r rho^(n_r - 1) / prod_i m_ir!`` for every reaction.

        Zero for a reaction with a reactant outside the evolved set.
        """
        other, fitted = self._temperature_parts(t9)
        if len(self.captures):
            fitted = fitted.copy()
            for k in self.captures:
                if fitted[k]:
                    fitted[k] *= (float(rho) * float(ye)) ** int(self.electrons[k])
        rate = other + fitted
        if screening is not None:
            for k in np.nonzero(self.active)[0]:
                r = self.reactions[k]
                try:
                    rate[k] *= float(screening(r, float(t9), float(rho), ye))
                except TypeError:
                    rate[k] *= float(screening(r))
        powers = [float(rho) ** p for p in range(int(self.density_power.max(initial=0)) + 1)]
        with np.errstate(over="ignore", invalid="ignore"):
            prefactor = rate * np.asarray(powers)[self.density_power] / self.statistical
        prefactor[~self.active] = 0.0
        return prefactor

    def extended(self, y: np.ndarray) -> np.ndarray:
        """The abundances with a final entry of 1 for the unused reactant slots."""
        self._abundance[:self.n] = y
        return self._abundance

    def _powers(self, yext: np.ndarray, k: np.ndarray, slot: int, powers) -> np.ndarray:
        """``Y ** power`` with Python's power, for reactions ``k`` at ``slot``."""
        out = []
        for base, power in zip(yext[self.slot_index[k, slot]].tolist(), list(powers)):
            try:
                out.append(base ** int(power))
            except OverflowError:
                out.append(math.inf)
        return np.array(out, dtype=float)

    def flows(self, y: np.ndarray, prefactor: np.ndarray) -> np.ndarray:
        """``F_r = prefactor_r prod_i Y_i^m_ir``, in the order of Reaction.flux."""
        yext = self.extended(y)
        flow = prefactor.copy()
        with np.errstate(over="ignore", invalid="ignore"):
            for j in range(self.slot_index.shape[1]):
                k = self._single[j]
                flow[k] *= yext[self.slot_index[k, j]]
                k = self._repeated[j]
                if len(k):
                    flow[k] *= self._powers(yext, k, j, self.slot_power[k, j])
        return flow

    def flow_derivatives(self, y: np.ndarray, prefactor: np.ndarray):
        """``dF_r/dY_j`` as a sparse (reactions x species) matrix."""
        from scipy.sparse import csr_matrix
        yext = self.extended(y)
        width = self.slot_index.shape[1]
        rows, cols, values = [], [], []
        with np.errstate(over="ignore", invalid="ignore"):
            for pos in range(width):
                power = self.slot_power[:, pos]
                k = np.nonzero((power > 0) & (prefactor != 0.0))[0]
                if not len(k):
                    continue
                # The factor of this slot replaced by m Y^(m - 1), the others
                # kept, multiplied in the same order as the flow itself.
                term = prefactor[k] * power[k]
                # Y^0 = 1 needs no multiplication and Y^1 = Y exactly; only a
                # reactant that appears three or more times needs a power.
                twice, more = power[k] == 2, power[k] > 2
                term[twice] *= yext[self.slot_index[k[twice], pos]]
                if more.any():
                    term[more] *= self._powers(yext, k[more], pos, power[k][more] - 1)
                for other in range(width):
                    if other != pos:
                        m = self.slot_power[k, other]
                        single, repeated = m == 1, m > 1
                        term[single] *= yext[self.slot_index[k[single], other]]
                        if repeated.any():
                            term[repeated] *= self._powers(yext, k[repeated], other, m[repeated])
                keep = term != 0.0
                rows.append(k[keep])
                cols.append(self.slot_index[k[keep], pos])
                values.append(term[keep])
        if rows:
            rows, cols, values = np.concatenate(rows), np.concatenate(cols), np.concatenate(values)
        return csr_matrix((values, (rows, cols)), shape=(self.m, self.n))


def rhs(network: Network, species: Sequence[str], thermo: ThermoFunction, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None):
    """Return ``f(t, y)``, the time derivative ``dY/dt`` of the abundances ``y``.

    ``y`` holds the abundances of ``species``, in that order.  Negative
    abundances are treated as zero.  ``thermo(t, abundances)`` gives T9 and
    the density.  Each call evaluates every reaction among the species at
    once, as arrays (see :class:`_NetworkArrays`), and agrees with
    :meth:`ReactionNetwork.ydot` to rounding error.
    """
    from .reactions import _refresh_screening
    species = [normalize_species_name(s) for s in species]
    arrays = _NetworkArrays(network, species)
    idx = arrays.index

    def f(t: float, yvec: np.ndarray) -> np.ndarray:
        yclip = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = dict(zip(species, yclip.tolist()))
        t9, rho = thermo(float(t), abund)
        ye = arrays.electron_fraction(yclip)
        if screening is not None:
            _refresh_screening(screening, abund, t9, rho, ye)
        prefactor = arrays.prefactors(t9, rho, ye, screening)
        out = arrays.stoichiometry @ arrays.flows(yclip, prefactor)
        if weak_rates:
            for wr in weak_rates:
                parent = getattr(wr, 'parent', None)
                daughter = getattr(wr, 'daughter', None)
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
            yp = y.copy()
            yp[k] += h
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

    The derivatives of all the flows are evaluated at once, as arrays, and
    the matrix is their product with the sparse stoichiometric matrix, so it
    is built sparse from the start.  The finite-difference Jacobian of
    :func:`jacobian` instead costs ``N+1`` evaluations of the full right-hand
    side.

    ``lambda_r`` is treated as independent of composition.  With a screening
    function or weak-rate tables that depend on ``Ye`` the result is an
    approximate Jacobian: it changes the Newton convergence rate of an implicit
    solver but not the solution it converges to.

    Returns a SciPy CSC matrix, or a dense array when ``sparse`` is false.
    """
    from scipy.sparse import coo_matrix
    species = [normalize_species_name(s) for s in species]
    arrays = _NetworkArrays(network, species)
    idx = arrays.index
    n = len(species)

    def j(t: float, yvec: np.ndarray):
        y = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = dict(zip(species, y.tolist()))
        t9, rho = thermo(float(t), abund)
        ye = arrays.electron_fraction(y)
        prefactor = arrays.prefactors(t9, rho, ye, screening)
        J = arrays.stoichiometry @ arrays.flow_derivatives(y, prefactor)
        if weak_rates:
            rows, cols, values = [], [], []
            for wr in weak_rates:
                parent = getattr(wr, 'parent', None)
                daughter = getattr(wr, 'daughter', None)
                if parent in idx and daughter in idx:
                    rate = wr.rate(t9, rho * max(ye, 1e-30))
                    rows += [idx[parent], idx[daughter]]
                    cols += [idx[parent], idx[parent]]
                    values += [-rate, rate]
            if rows:
                J = J + coo_matrix((values, (rows, cols)), shape=(n, n))
        J = J.tocsc()
        J.eliminate_zeros()
        J.sort_indices()
        return J if sparse else J.toarray()
    return j


def jacobian_sparsity(network: Network, species: Sequence[str]):
    """Return the stoichiometric sparsity pattern, as a boolean CSC matrix.

    Entry ``(i, j)`` is set when species ``j`` is a reactant of a reaction
    that changes species ``i``.
    """
    from scipy.sparse import csc_matrix
    species = [normalize_species_name(s) for s in species]
    idx = {s: i for i, s in enumerate(species)}
    rows, cols = [], []
    for r in network.reactions.reactions:
        affected = [idx[n] for n in r.stoichiometry() if n in idx]
        deps = [idx[p.species] for p in r.reactants if p.species in idx]
        for a in affected:
            for d in deps:
                rows.append(a)
                cols.append(d)
    pattern = csc_matrix((np.ones(len(rows), dtype=bool), (rows, cols)),
                         shape=(len(species), len(species)), dtype=bool)
    pattern.sum_duplicates()
    return pattern


def evolve_zone(network: Network, zone: Zone, times: Sequence[float], thermo: Optional[ThermoFunction] = None, method: str = "bdf", species: Optional[Sequence[str]] = None, screening: Optional[ScreeningFunction] = None, weak_rates: Optional[Sequence[object]] = None, rtol: float = 1e-6, atol: float = 1e-30, use_jacobian: bool = True, project_positive: bool = True, jac_mode: str = "analytic", max_change: float = 0.005, min_mass_fraction: float = 1e-15, progress: Optional[Callable[[float], None]] = None) -> EvolutionResult:
    """Evolve one zone over ``times``.

    ``progress``, if given, is called with the time the integrator has reached:
    on every right-hand-side evaluation for the SciPy methods, and after every
    accepted step for ``adaptive_implicit_euler``.  It lets a caller report on
    a long integration; it must not change anything.

    ``method="adaptive_implicit_euler"`` takes implicit (backward) Euler steps
    solved by Newton-Raphson iterations, the scheme NucNet Tools uses, and
    chooses each step so that no abundance whose mass fraction exceeds
    ``min_mass_fraction`` changes by more than ``max_change`` (as a fraction)
    in one step.  It is first-order accurate, so ``max_change`` sets the
    accuracy: halving it halves the error.  It is the robust choice for large
    networks in which a fuel runs out.  There the higher-order BDF formula
    can predict a slightly negative abundance of the fuel; the right-hand
    side treats it as zero, which stops every reaction that consumes the
    fuel, and SciPy's step size collapses.  Implicit Euler cannot turn a
    decaying abundance negative.  In explosive burning, ``max_change=0.05``
    with ``min_mass_fraction=1e-10`` is usually accurate to a fraction of a per
    cent and many times faster than the defaults; check it by repeating the
    run with a smaller ``max_change``.  ``rtol``, ``atol``, ``use_jacobian`` and
    ``jac_mode`` do not apply to it.

    ``jac_mode`` selects how the Jacobian is supplied to the implicit SciPy
    solvers.  ``"analytic"`` (the default) uses :func:`analytic_jacobian`, which
    differentiates the stoichiometric flows in closed form for all the
    reactions at once.  ``"numerical"`` uses the finite-difference
    :func:`jacobian`, costing ``N+1`` right-hand-side evaluations.  ``"sparsity"``
    supplies no Jacobian and instead passes the stoichiometric sparsity pattern
    so that SciPy estimates the matrix by grouped finite differences.

    Note that SciPy applies ``jac_sparsity`` only when no Jacobian callable is
    given, so the pattern and an explicit Jacobian are alternatives rather than
    complements.

    ``method`` is one of :data:`METHODS`, in any case.  An unrecognised name
    raises ``ValueError`` rather than falling back to another integrator.
    """
    ts = np.asarray(times, dtype=float)
    if ts.ndim != 1 or len(ts) < 2:
        raise ValueError("times must be a one-dimensional array with at least two points")
    method_l = str(method).lower()
    if method_l not in METHODS:
        raise ValueError(f"unknown method {method!r}; choose one of: {', '.join(METHODS)}")
    if str(jac_mode).lower() not in _JAC_MODES:
        raise ValueError(f"unknown jac_mode {jac_mode!r}; choose one of: {', '.join(_JAC_MODES)}")
    # Photons and leptons are part of the reaction records but not of the
    # abundance vector; evolving them would integrate a meaningless quantity.
    species = [normalize_species_name(s) for s in (species or network.species_names())]
    species = [s for s in species if not is_massless(s)]
    y0 = np.array([zone.get_abundance(s) for s in species], dtype=float)
    thermo = thermo or zone_thermo(zone)
    if method_l in _ADAPTIVE_METHODS:
        return _adaptive_implicit_euler(network, species, y0, ts, thermo, screening,
                                        weak_rates, max_change, min_mass_fraction,
                                        progress=progress)
    f = rhs(network, species, thermo, screening=screening, weak_rates=weak_rates)
    if progress is not None:
        _rhs = f

        def f(t, y):
            progress(t)
            return _rhs(t, y)
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


def _adaptive_implicit_euler(network: Network, species: Sequence[str], y0: np.ndarray,
                             ts: np.ndarray, thermo: ThermoFunction, screening,
                             weak_rates, max_change: float, min_mass_fraction: float,
                             newton_tol: float = 1e-7, max_newton: int = 12,
                             jac_every: int = 4, max_growth: float = 1.5,
                             max_rejections: int = 200,
                             progress: Optional[Callable[[float], None]] = None) -> EvolutionResult:
    """Implicit Euler steps with Newton-Raphson iterations and step control.

    Each step solves ``y - y_old - h f(t + h, y) = 0`` by Newton iterations
    with the matrix ``I - h J``.  The Jacobian is kept for up to ``jac_every``
    steps and the factorised matrix while the step changes by less than 30
    per cent; both are recomputed at the current iterate every third
    iteration while the iterations converge slowly.  The right-hand side
    treats a negative abundance as zero, so an iterate that dips below zero
    would make the iterations chatter about zero; since the implicit Euler
    solution is itself never negative, the iterates are kept at or above zero.

    A step is accepted when the iterations converge and no abundance whose
    mass fraction was above ``min_mass_fraction`` changed by more than twice
    ``max_change``; the next step then aims at ``max_change``, growing by at
    most ``max_growth``.  Only abundances already above the threshold limit
    the step, as in NucNet Tools: one just being made from nothing would
    otherwise count as a 100 per cent change.  Steps end exactly on ``ts``.
    """
    from scipy.sparse import identity
    from scipy.sparse.linalg import splu
    from .species import Species

    if not 0.0 < float(max_change) < 1.0:
        raise ValueError(f"max_change must lie between 0 and 1, not {max_change!r}")
    if float(min_mass_fraction) < 0.0:
        raise ValueError("min_mass_fraction must not be negative")
    if np.any(np.diff(ts) <= 0.0):
        raise ValueError("times must increase strictly for adaptive_implicit_euler")
    f = rhs(network, species, thermo, screening=screening, weak_rates=weak_rates)
    jac = analytic_jacobian(network, species, thermo, screening=screening,
                            weak_rates=weak_rates, sparse=True)
    mass = np.array([(network.species.get(s) or Species.parse(s)).a for s in species],
                    dtype=float)
    y_min = float(min_mass_fraction) / np.maximum(mass, 1.0)
    eye = identity(len(species), format="csc")

    def factor(matrix, step):
        newton = (eye - step * matrix).tocsc()
        try:
            # Diagonal pivots keep the fill-reducing ordering of the
            # symmetrised pattern.  The Newton matrix of a reaction network
            # needs no other pivots, and on a complete JINA database this
            # leaves a fifth of the fill-in and makes each factorisation
            # several times faster.  A vanishing pivot falls back to the
            # default threshold pivoting.
            return splu(newton, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                        options=dict(SymmetricMode=True))
        except RuntimeError:
            return splu(newton, permc_spec="MMD_AT_PLUS_A")

    y = np.clip(np.asarray(y0, dtype=float), 0.0, None)
    out = np.zeros((len(ts), len(species)))
    out[0] = y
    t = float(ts[0])
    dt = 1.0e-3 * float(ts[1] - ts[0])
    steps = rejected = nfev = njev = 0
    in_a_row = 0
    jm, jm_age = None, 0
    lu, lu_h = None, None
    k = 1
    while k < len(ts):
        target = float(ts[k])
        remaining = target - t
        # A step that would stop just short of the next output time is
        # stretched to reach it, so that no tiny leftover step remains.
        last = dt >= remaining * (1.0 - 1e-10)
        h = remaining if last else dt
        if jm is None or jm_age >= jac_every:
            jm, jm_age, lu = jac(t + h, y), 0, None
            njev += 1
        if lu is None or abs(h / lu_h - 1.0) > 0.3:
            lu, lu_h = factor(jm, h), h
        ynew = y.copy()
        converged = False
        for it in range(max_newton):
            g = ynew - y - h * f(t + h, ynew)
            nfev += 1
            delta = lu.solve(-g)
            ynew = ynew + delta
            np.maximum(ynew, 0.0, out=ynew)
            if np.max(np.abs(delta) / (np.abs(ynew) + y_min)) < newton_tol:
                converged = True
                break
            if it % 3 == 2:
                jm, jm_age = jac(t + h, ynew), 0
                njev += 1
                lu, lu_h = factor(jm, h), h
        if converged and np.all(np.isfinite(ynew)):
            big = y >= y_min
            change = (float(np.max(np.abs(ynew[big] - y[big])
                                   / np.maximum(y[big], ynew[big])))
                      if np.any(big) else 0.0)
            if change <= 2.0 * max_change:
                t = target if last else t + h
                y = ynew
                steps += 1
                if progress is not None:
                    progress(t)
                jm_age += 1
                in_a_row = 0
                grow = max_growth if change == 0.0 else min(max_growth,
                                                            0.9 * max_change / change)
                proposal = h * max(grow, 0.5)
                dt = max(proposal, dt) if h < dt else proposal
                if last:
                    out[k] = y
                    k += 1
                continue
            dt = h * max(0.1, 0.9 * max_change / change)
        else:
            dt = 0.5 * h
        rejected += 1
        in_a_row += 1
        if in_a_row > max_rejections or not np.isfinite(dt) or dt <= 0.0:
            return EvolutionResult(ts[:k], list(species), out[:k], False,
                                   f"adaptive_implicit_euler: {in_a_row} steps rejected "
                                   f"in a row at t = {t:.6e}",
                                   nfev, njev, steps, rejected)
    return EvolutionResult(ts, list(species), out, True,
                           "adaptive_implicit_euler reached the end of the interval",
                           nfev, njev, steps, rejected)


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
            yp = yn.copy()
            yp[k] += eps * max(abs(yn[k]), 1.0)
            J[:, k] -= h * (f(t + h, yp) - base_f) / (yp[k] - yn[k])
        try:
            step = np.linalg.solve(J, -g)
        except np.linalg.LinAlgError:
            step = -g
        yn = yn + step
    return yn


def evolve_network_zones(network: Network, times: Sequence[float], thermo: Optional[ThermoFunction] = None, method: str = "bdf", **kwargs) -> List[EvolutionResult]:
    return [evolve_zone(network, z, times, thermo=thermo, method=method, **kwargs) for z in network.zones]
