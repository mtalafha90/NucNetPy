"""Quasi-statistical equilibrium (QSE) with constrained clusters.

A pure-Python port of the libnuceq cluster workflow used by the NucNet Tools
examples (``examples/analysis/compare_equil.cpp``) and the blog posts on
cluster flows, (n,gamma)-(gamma,n) equilibrium, and comparing network
calculations to equilibrium.

In libnuceq an equilibrium *cluster* is a subset of species (selected by an
XPath expression such as ``[z >= 6]``) whose total abundance is constrained to
a prescribed value.  Each cluster gets one extra chemical-potential multiplier
``lambda_c`` on top of the global proton/neutron potentials:

    Y_i = pref_i * exp(Z_i mu_p + N_i mu_n + lambda_{c(i)})

subject to ``sum A_i Y_i = 1``, ``sum Z_i Y_i = Ye``, and
``sum_{i in c} Y_i = Y_c`` for every cluster ``c``.  Full NSE is the special
case with no clusters; an (n,gamma)-(gamma,n) r-process equilibrium is one
cluster holding all heavy nuclei with ``Y_c`` equal to the current heavy-nuclei
abundance.

Cluster membership is by species name; use
:func:`nucnetpy.network_limiter.select_species` for XPath-like Z/A selection.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .core import Network, Zone
from .nse import _log_prefactor
from .species import Species, normalize_species_name


@dataclass
class QSECluster:
    """A constrained equilibrium cluster: species names and their total Y."""
    species: Sequence[str]
    constraint: float
    label: str = ""


@dataclass
class QSEResult:
    t9: float
    rho: float
    ye: float
    mu_p: float
    mu_n: float
    lambdas: List[float]
    abundances: Dict[str, float]
    success: bool
    message: str = ""

    def zone(self) -> Zone:
        return Zone(abundances=dict(self.abundances),
                    properties={"t9": str(self.t9), "rho": str(self.rho), "ye": str(self.ye)})

    @property
    def xsum(self) -> float:
        return float(sum(Species.parse(k).a * v for k, v in self.abundances.items()))

    @property
    def computed_ye(self) -> float:
        return float(sum(Species.parse(k).z * v for k, v in self.abundances.items()))


def solve_qse(network: Network, t9: float, rho: float, ye: float, clusters: Sequence[QSECluster], species: Optional[Sequence[str]] = None, include_partition: bool = True, nse_correction=None, tol: float = 1e-8, max_iter: int = 300, require_nuclear_data: bool = True) -> QSEResult:
    """Solve a constrained (QSE) equilibrium for a network.

    ``clusters`` is a sequence of :class:`QSECluster`; species not in any
    cluster follow unconstrained NSE.  A species may belong to at most one
    cluster.  With ``clusters=[]`` the result equals :func:`~nucnetpy.nse.solve_nse`,
    and the same solver is used, so ``success`` means the residual norm is
    below ``tol`` and non-unique potentials are reported as the minimum-norm
    solution.  Species without nuclear data are left out, as in
    :func:`~nucnetpy.nse.solve_nse`, unless ``require_nuclear_data`` is false.
    """
    from .nse import _equilibrium_species, _solve_equilibrium

    sps = _equilibrium_species(network, species, require_nuclear_data)
    if not sps:
        raise ValueError("No valid species available for QSE solve")

    idx = {sp.name: k for k, sp in enumerate(sps)}
    n_c = len(clusters)
    member = np.full(len(sps), -1, dtype=int)  # cluster index per species, -1 = free
    for c, cl in enumerate(clusters):
        if float(cl.constraint) <= 0.0:
            raise ValueError(f"cluster {c} constraint must be positive")
        for name in cl.species:
            k = idx.get(normalize_species_name(name))
            if k is None:
                continue
            if member[k] >= 0:
                raise ValueError(f"species {sps[k].name} assigned to two clusters")
            member[k] = c
        if not np.any(member == c):
            raise ValueError(f"cluster {c} contains none of the species being solved for")

    log_pref = np.array([_log_prefactor(sp, t9, rho, include_partition=include_partition) for sp in sps], dtype=float)
    if nse_correction is not None:
        log_pref = log_pref + np.array([float(nse_correction(sp, t9, rho, float(ye))) for sp in sps], dtype=float)
    z = np.array([sp.z for sp in sps], dtype=float)
    a = np.array([sp.a for sp in sps], dtype=float)
    ye = float(ye)
    log_constraints = [math.log(float(cl.constraint)) for cl in clusters]

    x, log_y, norm, message = _solve_equilibrium(log_pref, z, a, ye, cluster=member,
                                                 log_constraints=log_constraints,
                                                 tol=tol, max_iter=max_iter)
    abund = {sp.name: float(math.exp(min(float(v), 700.0))) for sp, v in zip(sps, log_y)}
    return QSEResult(t9, rho, ye, float(x[0]), float(x[1]), [float(v) for v in x[2:]],
                     abund, norm < tol, message)


def cluster_abundance(abundances: Mapping[str, float], cluster_species: Sequence[str]) -> float:
    """Return ``Y_c``, the summed abundance of the cluster members."""
    wanted = {normalize_species_name(s) for s in cluster_species}
    return float(sum(y for name, y in abundances.items() if name in wanted))


def cluster_ydot(network: Network, abundances: Mapping[str, float], cluster_species: Sequence[str], t9: float, rho: float = 1.0, use_reverse: bool = False) -> Dict[str, float]:
    """Return per-reaction contributions to ``dY_c/dt`` for a cluster.

    Each reaction contributes ``flux * (net change in cluster member count)``;
    reactions that only shuffle abundance inside the cluster contribute nothing
    (blog: "Calculating cluster flows", C++ ``compute_Ycdot``).  With
    ``use_reverse=True`` the net (forward minus detailed-balance reverse) flux
    is used.  ``sum(result.values())`` is dY_c/dt.
    """
    wanted = {normalize_species_name(s) for s in cluster_species}
    flows = None
    if use_reverse:
        from .detailed_balance import net_flows
        flows = net_flows(network, abundances, t9=t9, rho=rho)
    out: Dict[str, float] = {}
    for r in network.reactions.reactions:
        delta = 0
        for p in r.products:
            if p.species in wanted:
                delta += p.count
        for p in r.reactants:
            if p.species in wanted:
                delta -= p.count
        if delta == 0:
            continue
        flux = flows[r.string][2] if flows is not None else r.flux(abundances, t9=t9, rho=rho)
        out[r.string] = float(delta * flux)
    return out
