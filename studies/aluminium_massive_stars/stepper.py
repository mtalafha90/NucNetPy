"""A NucNet Tools style time stepper built from NucNetPy's own pieces.

``nucnetpy.evolve_zone`` hands the network to SciPy's BDF integrator.  On the
891-nuclide network of this study that fails in the first stage: as hydrogen
runs out, the higher-order BDF formula predicts a small negative hydrogen
abundance, NucNetPy's right-hand side clips it to zero and so stops all
proton captures, and the step size collapses (see README.md).

NucNet Tools, which the paper used, avoids this by taking implicit (backward)
Euler steps, solved with Newton-Raphson iterations, and by sizing each step
so that no abundance changes by more than a set fraction.  Implicit Euler
cannot drive a decaying abundance negative, which makes it robust at fuel
exhaustion; it is first-order accurate, so the step control sets the
accuracy.  This module implements that scheme around NucNetPy's right-hand
side and analytic Jacobian (``nucnetpy.solver.rhs`` and
``nucnetpy.solver.analytic_jacobian``), so the rates, partition functions,
screening and flow formula are all NucNetPy's.

Each step solves ``y - y_old - dt * f(y) = 0`` for the new abundances ``y``
by Newton iterations with the matrix ``I - dt * J``.  The Jacobian ``J`` is
kept for up to ``jac_every`` steps and recomputed sooner if the iterations
slow down.  A step is accepted when the iterations converge and no abundance
whose mass fraction was above ``x_min`` changes by more than
``2 * max_change``; the next step is then sized to aim at ``max_change``.
Steps also end exactly on the requested output times.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence

import numpy as np
from scipy.sparse import identity
from scipy.sparse.linalg import splu

from nucnetpy.solver import analytic_jacobian, rhs
from nucnetpy.species import is_massless


@dataclass
class StepperResult:
    time: np.ndarray
    species: List[str]
    y: np.ndarray
    success: bool
    message: str
    steps: int = 0
    rejected: int = 0
    nfev: int = 0
    njev: int = 0
    #: Largest mass fraction that came out negative from a converged step and
    #: was set to zero.  Implicit Euler keeps this at round-off level.
    max_negative_clipped: float = 0.0

    @property
    def final_abundances(self):
        return {s: float(v) for s, v in zip(self.species, self.y[-1])}


def evolve(network, abundances, times: Sequence[float], t9: float, rho: float,
           screening=None, max_change: float = 0.02, x_min: float = 1.0e-15,
           newton_tol: float = 1.0e-7, max_newton: int = 12, jac_every: int = 4,
           first_step: Optional[float] = None, max_growth: float = 1.5,
           progress: Optional[Callable[[float, int], None]] = None) -> StepperResult:
    """Integrate one zone at constant ``t9`` and ``rho`` over ``times``.

    ``abundances`` maps species names to abundances Y (mass fraction divided
    by mass number).  ``times`` must start at the initial time and increase.
    The first step defaults to a thousandth of the first output interval.
    """
    species = [s for s in network.species_names() if not is_massless(s)]
    a = np.array([network.species[s].a for s in species], dtype=float)
    y = np.array([float(abundances.get(s, 0.0)) for s in species])
    thermo = lambda t, ab: (float(t9), float(rho))  # noqa: E731
    f = rhs(network, species, thermo, screening=screening)
    jac = analytic_jacobian(network, species, thermo, screening=screening, sparse=True)
    eye = identity(len(species), format="csc")
    y_min = x_min / a            # abundance threshold for the step control

    times = np.asarray(times, dtype=float)
    out = np.zeros((len(times), len(species)))
    out[0] = y
    t = float(times[0])
    dt = float(first_step) if first_step else 1.0e-3 * (times[1] - times[0])
    steps = rejected = nfev = njev = 0
    worst_negative = 0.0
    jm, jm_age = None, 0
    k = 1
    while k < len(times):
        target = times[k]
        h = min(dt, target - t)
        if jm is None or jm_age >= jac_every:
            jm, jm_age = jac(t + h, y), 0
            njev += 1
        fresh = jm_age == 0
        lu = splu((eye - h * jm).tocsc())
        ynew = y.copy()
        converged = False
        for it in range(max_newton):
            g = ynew - y - h * f(t + h, ynew)
            nfev += 1
            delta = lu.solve(-g)
            ynew = ynew + delta
            if np.max(np.abs(delta) / (np.abs(ynew) + y_min)) < newton_tol:
                converged = True
                break
            if it == 3 and not fresh:
                # Slow convergence with an old Jacobian: refresh it.
                jm, jm_age, fresh = jac(t + h, np.clip(ynew, 0.0, None)), 0, True
                njev += 1
                lu = splu((eye - h * jm).tocsc())
        if not converged or not np.all(np.isfinite(ynew)):
            rejected += 1
            if not fresh:
                jm = None                 # retry the same step with a new Jacobian
                continue
            dt = 0.25 * h
            if dt < 1e-30 * max(1.0, t):
                return StepperResult(times[:k], species, out[:k], False,
                                     f"step size underflow at t = {t:.6e} s",
                                     steps, rejected, nfev, njev, worst_negative)
            continue
        neg = ynew < 0.0
        if np.any(neg):
            worst_negative = max(worst_negative, float(np.max(-ynew[neg] * a[neg])))
            ynew[neg] = 0.0
        # Only abundances already above the threshold limit the step, as in
        # NucNet Tools: one that is just being made from nothing would
        # otherwise count as a 100 per cent change and stall the integration.
        big = y >= y_min
        change = (float(np.max(np.abs(ynew[big] - y[big]) / np.maximum(y[big], ynew[big])))
                  if np.any(big) else 0.0)
        if change > 2.0 * max_change:
            dt = h * max(0.1, 0.9 * max_change / change)
            rejected += 1
            continue
        t += h
        y = ynew
        steps += 1
        jm_age += 1
        grow = max_growth if change == 0.0 else min(max_growth, 0.9 * max_change / change)
        proposal = h * max(grow, 0.5)
        dt = max(proposal, dt) if h < dt else proposal
        if abs(t - target) <= 1e-12 * max(1.0, abs(target)):
            t = float(target)
            out[k] = y
            k += 1
        if progress is not None:
            progress(t, steps)
    return StepperResult(times, species, out, True, "reached the end of the stage",
                         steps, rejected, nfev, njev, worst_negative)
