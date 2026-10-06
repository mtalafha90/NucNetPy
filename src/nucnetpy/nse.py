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
from typing import Dict, List, Optional, Sequence, Tuple
import math
import numpy as np

from .constants import AVOGADRO, KB_MEV, HBAR_C_MEV_FM, AMU_MEV
from .core import Network, Zone
from .species import Species, normalize_species_name

# hbar*c in MeV cm
_HBAR_C_MEV_CM = HBAR_C_MEV_FM * 1.0e-13


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
    kt = KB_MEV * max(float(t9), 1e-30) * 1.0e9  # kT in MeV
    # More tightly bound nuclei have lower mass excess and are enhanced.  The
    # Z*ME(p)+N*ME(n) piece of the binding energy is linear in Z and N and is
    # absorbed by the proton/neutron chemical potentials, so only -ME_i remains.
    return _log_quantum_abundance(sp, t9, rho, include_partition) - float(sp.mass_excess) / kt


def _log_quantum_abundance(sp: Species, t9: float, rho: float, include_partition: bool = True) -> float:
    """Natural log of ``Y_Q = (2J+1) G(T) n_Q / (rho N_A)``.

    ``n_Q = (A m_u kT / 2 pi hbar^2)^(3/2)`` is the quantum concentration.  This
    is the NSE prefactor without the binding-energy term, and it is also the
    quantity in the Sackur-Tetrode entropy (:func:`nucnetpy.thermo.entropy_per_nucleon`).
    """
    t9 = max(float(t9), 1e-30)
    rho = max(float(rho), 1e-300)
    kt = KB_MEV * t9 * 1.0e9  # kT in MeV
    a = max(int(sp.a), 1)
    g = _partition(sp, t9) if include_partition else 1.0
    # libnucnet/JINA files store the *normalised* partition function
    # G(T) = Z(T)/(2J_0+1), which tends to one as T -> 0 for every nuclide.
    # The ground-state spin degeneracy is therefore not contained in the table
    # and has to be supplied here to build the full statistical weight
    # Z(T) = (2J_0+1) G(T).  Omitting it biases the NSE abundance of every
    # nuclide with non-zero ground-state spin by that factor.  A species whose
    # spin is unknown (``None``) keeps a weight of one rather than a guess.
    if sp.spin is not None:
        g *= 2.0 * float(sp.spin) + 1.0
    # Quantum concentration (cm^-3) for a nucleus of mass A*m_u.  All energies
    # are in MeV and hbar*c is in MeV*cm, so the units are consistent; the
    # nucleon mass enters as its rest energy m_u c^2 = AMU_MEV, not its mass in
    # grams.  The species mass is folded in through the a**1.5 term.
    log_theta = 1.5 * math.log((AMU_MEV * kt) / (2.0 * math.pi * (_HBAR_C_MEV_CM ** 2)))
    return math.log(max(g, 1e-300)) + 1.5 * math.log(a) + log_theta - math.log(rho * AVOGADRO)


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


def _equilibrium_species(network: Network, species: Optional[Sequence[str]],
                         require_nuclear_data: bool) -> List[Species]:
    """Return the nuclides an equilibrium is solved over.

    Photons and leptons carry no baryon number and are left out.  Species
    invented to satisfy a reaction record carry a placeholder mass excess of
    zero; including them would let an unbound nuclide compete with the iron
    peak, so they are dropped unless ``require_nuclear_data`` is false.
    """
    names = [normalize_species_name(s) for s in (species or network.species_names())]
    if require_nuclear_data:
        placeholders = set(network.species_without_nuclear_data())
        names = [n for n in names if n not in placeholders]
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
    return sps


def _log_sum_exp(values: np.ndarray) -> Tuple[float, np.ndarray]:
    """Return ``log(sum(exp(values)))`` and the normalised weights, without overflow."""
    top = float(np.max(values))
    scaled = np.exp(values - top)
    total = float(np.sum(scaled))
    return top + math.log(total), scaled / total


def _solve_equilibrium(log_pref: np.ndarray, z: np.ndarray, a: np.ndarray, ye: float,
                       cluster: Optional[np.ndarray] = None,
                       log_constraints: Sequence[float] = (),
                       tol: float = 1e-8, max_iter: int = 200):
    """Find the chemical potentials of an equilibrium, optionally constrained.

    This is the solver behind both :func:`solve_nse` and
    :func:`nucnetpy.qse.solve_qse`.  Species ``i`` has

        ln Y_i = log_pref_i + B_i . x,     x = (mu_p, mu_n, lambda_1, ..., lambda_k),

    where ``B_i = (Z_i, N_i, e_c)`` and ``e_c`` marks the cluster, if any, that
    species ``i`` belongs to.  The residuals are all of order one:

        ln(sum_i A_i Y_i)                      so that sum A Y = 1
        ln(P) - ln(M)                          so that sum Z Y = Ye
        ln(sum_{i in c} Y_i) - ln(Y_c)         one per cluster

    Charge balance is written as ``P = M``, where ``P`` sums the charge excess
    ``(Z_i - Ye A_i) Y_i`` of the species with ``Z_i > Ye A_i`` and ``M`` the
    deficit of those with ``Z_i < Ye A_i``.  Species with ``Z_i = Ye A_i``
    cannot affect the balance and drop out exactly.  This matters at
    ``Ye = 0.5``: the bulk is then in N = Z nuclei, and the balance is decided
    by trace free nucleons whose effect on ``sum Z Y - Ye`` is far below
    floating-point resolution, so that form leaves their abundances
    undetermined.  The logarithmic form stays of order one however small they
    are.  When only one side exists, ``Ye`` cannot be met exactly and the
    linear form ``sum Z Y / sum A Y - Ye`` is used instead.

    The Jacobian is computed exactly.  A Levenberg-Marquardt solve from several
    starting points is followed by Gauss-Newton polishing, so the result is
    converged to machine precision rather than to wherever the optimiser's
    step-size test happened to stop.

    The potentials need not be unique.  If every species has the same ratio
    N/Z -- an alpha chain without free nucleons, for example -- only the
    combination ``Z mu_p + N mu_n`` is fixed by the constraints.  Any split
    gives the same abundances, and which split an optimiser lands on depends
    on its version.  The minimum-norm solution is therefore returned whenever
    ``B`` is rank-deficient, which makes the reported potentials reproducible.

    Returns ``(x, ln_y, residual_norm, message)``.
    """
    from scipy.optimize import least_squares

    k = len(log_constraints)
    cluster = np.full(len(z), -1, dtype=int) if cluster is None else np.asarray(cluster)
    basis = np.zeros((len(z), 2 + k), dtype=float)
    basis[:, 0] = z
    basis[:, 1] = a - z
    members = [cluster == c for c in range(k)]
    for c, mask in enumerate(members):
        basis[mask, 2 + c] = 1.0
    log_a = np.log(a)
    charge_per_nucleon = z / a
    log_constraints = np.asarray(log_constraints, dtype=float)
    excess = z - ye * a
    excess[np.abs(excess) < 1e-12 * a] = 0.0     # round-off from a Ye such as 0.45
    surplus, deficit = excess > 0.0, excess < 0.0
    two_sided = bool(surplus.any() and deficit.any())
    if two_sided:
        log_surplus, log_deficit = np.log(excess[surplus]), np.log(-excess[deficit])

    def residual_and_jacobian(x):
        log_y = log_pref + basis @ x
        log_mass, w = _log_sum_exp(log_y + log_a)        # w_i = A_i Y_i / sum_j A_j Y_j
        mean_b = w @ basis
        if two_sided:
            log_p, p_weights = _log_sum_exp(log_y[surplus] + log_surplus)
            log_m, m_weights = _log_sum_exp(log_y[deficit] + log_deficit)
            charge_residual = log_p - log_m
            charge_row = p_weights @ basis[surplus] - m_weights @ basis[deficit]
        else:
            charge = float(charge_per_nucleon @ w)
            charge_residual = charge - ye
            charge_row = (charge_per_nucleon * w) @ basis - charge * mean_b
        residuals = [log_mass, charge_residual]
        rows = [mean_b, charge_row]
        for c, mask in enumerate(members):
            log_cluster, u = _log_sum_exp(log_y[mask])
            residuals.append(log_cluster - log_constraints[c])
            rows.append(u @ basis[mask])
        return np.array(residuals), np.array(rows)

    best = None
    for guess in (0.0, -1.0, -5.0, -10.0, -20.0):
        x0 = np.concatenate([[guess, guess], np.zeros(k)])
        sol = least_squares(lambda x: residual_and_jacobian(x)[0], x0,
                            jac=lambda x: residual_and_jacobian(x)[1],
                            method="lm", xtol=1e-14, ftol=1e-14, max_nfev=max_iter * 10)
        x, norm = _gauss_newton_polish(residual_and_jacobian, np.asarray(sol.x, dtype=float))
        if best is None or norm < best[1]:
            best = (x, norm, str(sol.message))
        if norm < tol:
            break
    x, norm, message = best

    if np.linalg.matrix_rank(basis) < basis.shape[1]:
        x = np.linalg.pinv(basis) @ (basis @ x)
        message += ("; the potentials are not unique for this species set, so the "
                    "minimum-norm values are reported")
    return x, log_pref + basis @ x, norm, message


def _gauss_newton_polish(residual_and_jacobian, x: np.ndarray, iterations: int = 20):
    """Refine a solution with Gauss-Newton steps while the residual keeps falling.

    Each step solves ``J dx = -r`` in the least-squares sense, which also copes
    with a rank-deficient Jacobian.  Returns ``(x, residual_norm)``.
    """
    r, j = residual_and_jacobian(x)
    norm = float(np.linalg.norm(r))
    for _ in range(iterations):
        step = np.linalg.lstsq(j, -r, rcond=None)[0]
        r_new, j_new = residual_and_jacobian(x + step)
        norm_new = float(np.linalg.norm(r_new))
        if not norm_new < norm:
            break
        x, r, j, norm = x + step, r_new, j_new, norm_new
    return x, norm


def solve_nse(network: Network, t9: float, rho: float, ye: float, species: Optional[Sequence[str]] = None, include_partition: bool = True, tol: float = 1e-8, max_iter: int = 200, nse_correction=None, require_nuclear_data: bool = True) -> NSEResult:
    """Solve an NSE composition for a network.

    The constraints are ``sum(A_i Y_i) = 1`` and ``sum(Z_i Y_i) = Ye``.  The
    solve works on well-scaled residuals with a numerically stable log-sum-exp,
    which keeps the iron-peak prefactors representable; see
    :func:`_solve_equilibrium` for the method.  ``success`` is true only when
    the residual norm is below ``tol``.

    When every species in the solve has the same ratio N/Z, only the
    combination ``Z mu_p + N mu_n`` is determined.  The abundances are still
    unique, and the reported ``mu_p`` and ``mu_n`` are the minimum-norm pair, so
    they do not depend on the SciPy version.

    ``nse_correction`` is an optional callable ``(species, t9, rho, ye) ->
    f_corr`` added to each species' NSE exponent, the libnucnet NSE correction
    factor hook.  Pass :func:`nucnetpy.coulomb.nse_correction` for the Bravo &
    Garcia-Senz Coulomb correction.
    """
    sps = _equilibrium_species(network, species, require_nuclear_data)
    if not sps:
        raise ValueError("No valid species available for NSE solve")
    log_pref = np.array([_log_prefactor(sp, t9, rho, include_partition=include_partition) for sp in sps], dtype=float)
    if nse_correction is not None:
        log_pref = log_pref + np.array([float(nse_correction(sp, t9, rho, float(ye))) for sp in sps], dtype=float)
    z = np.array([sp.z for sp in sps], dtype=float)
    a = np.array([sp.a for sp in sps], dtype=float)
    ye = float(ye)

    x, log_y, norm, message = _solve_equilibrium(log_pref, z, a, ye, tol=tol, max_iter=max_iter)
    abund = {sp.name: float(math.exp(min(v, 700.0))) for sp, v in zip(sps, log_y)}
    return NSEResult(t9, rho, ye, float(x[0]), float(x[1]), abund, norm < tol, message)


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
