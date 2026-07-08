"""Reverse reaction rates from detailed balance.

This replaces the libnucnet behaviour of computing reverse rates from the
forward rate, nuclear masses, and partition functions (blog workflows
"Comparing forward and reverse reaction rates", "Computing reaction flows",
and the (n,gamma)-(gamma,n) equilibrium studies).

For a reaction ``sum_r nu_r R -> sum_p nu_p P`` that conserves Z and A, the
equilibrium abundance ratio follows from the same Saha prefactors used by the
NSE solver: ``ln K = sum_p nu_p ln(pref_p) - sum_r nu_r ln(pref_r)`` where
``Y_i^eq = pref_i * exp(Z_i mu_p + N_i mu_n)`` and the chemical-potential terms
cancel for a balanced reaction.  Matching forward and reverse fluxes at
equilibrium gives

    lambda_rev = lambda_fwd * rho**(n_r - n_p) * (s_rev / s_fwd) / K

with ``n`` the reactant/product orders and ``s`` the duplicate-particle
statistical factors.  Because ``ln(pref)`` contains ``-ln(rho * N_A)``, the
explicit density powers cancel and the reverse rate is a pure function of
temperature, as it must be.

Photons and leptons (``gamma``, ``electron``, ...) are ignored on either side,
so the reverse of a radiative capture is the photodisintegration rate.
"""
from __future__ import annotations

import math
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from .nse import _log_prefactor
from .reactions import Reaction, TabularRate
from .species import Species, normalize_species_name

_LOG_MAX = 700.0


def _nuclear(species_map: Mapping[str, Species], name: str) -> Optional[Species]:
    """Return the Species for a participant, or None for photons/leptons."""
    sp = species_map.get(normalize_species_name(name))
    if sp is None:
        try:
            sp = Species.parse(name)
        except Exception:
            return None
    return sp if sp.a > 0 and sp.z >= 0 else None


def _side_stat_and_order(parts, species_map) -> Tuple[int, int]:
    stat = 1
    order = 0
    for p in parts:
        if _nuclear(species_map, p.species) is None:
            continue
        stat *= math.factorial(p.count)
        order += p.count
    return stat, order


def log_equilibrium_constant(reaction: Reaction, species_map: Mapping[str, Species], t9: float, rho: float = 1.0, include_partition: bool = True) -> float:
    """Return ``ln K`` with ``K = prod Y_p^eq / prod Y_r^eq`` for the reaction.

    Requires mass excesses (and optionally partition functions) on the species;
    the Q-value enters through the mass excesses as ``exp(Q/kT)``.
    """
    total = 0.0
    for p in reaction.products:
        sp = _nuclear(species_map, p.species)
        if sp is not None:
            total += p.count * _log_prefactor(sp, t9, rho, include_partition=include_partition)
    for r in reaction.reactants:
        sp = _nuclear(species_map, r.species)
        if sp is not None:
            total -= r.count * _log_prefactor(sp, t9, rho, include_partition=include_partition)
    return float(total)


def reverse_rate(reaction: Reaction, species_map: Mapping[str, Species], t9: float, rho: float = 1.0, forward: Optional[float] = None, include_partition: bool = True) -> float:
    """Return the detailed-balance reverse rate for ``reaction`` at ``t9``.

    ``forward`` overrides the forward rate (default ``reaction.rate(t9)``).
    The result is in the same convention as forward rates: to obtain a flux it
    must be combined with ``rho**(n_p - 1)``, product abundances, and the
    product-side statistical factor (see :func:`net_flows`).
    """
    lam_f = float(forward if forward is not None else reaction.rate(t9, rho=rho))
    if lam_f <= 0.0:
        return 0.0
    s_fwd, n_r = _side_stat_and_order(reaction.reactants, species_map)
    s_rev, n_p = _side_stat_and_order(reaction.products, species_map)
    log_k = log_equilibrium_constant(reaction, species_map, t9, rho, include_partition=include_partition)
    log_lam = math.log(lam_f) + (n_r - n_p) * math.log(max(float(rho), 1e-300)) + math.log(s_rev / s_fwd) - log_k
    if log_lam > _LOG_MAX:
        return float("inf")
    if log_lam < -_LOG_MAX:
        return 0.0
    return float(math.exp(log_lam))


def reverse_reaction(reaction: Reaction, species_map: Mapping[str, Species], t9_grid: Optional[Sequence[float]] = None, include_partition: bool = True) -> Reaction:
    """Return a Reaction for the reverse process with a tabulated rate.

    The rate is sampled from :func:`reverse_rate` on ``t9_grid`` (default 30
    points, T9 = 0.1 .. 10) so the reverse process — e.g. the (gamma,n)
    partner of an (n,gamma) capture — can be added to a network like any other
    reaction.
    """
    grid = np.geomspace(0.1, 10.0, 30) if t9_grid is None else np.asarray(t9_grid, dtype=float)
    rates = [reverse_rate(reaction, species_map, float(t9), include_partition=include_partition) for t9 in grid]
    return Reaction(
        reactants=list(reaction.products),
        products=list(reaction.reactants),
        tabular_rate=TabularRate(list(grid), rates),
        q_value=-reaction.q_value,
        source="detailed_balance",
        label=(reaction.label + "_reverse") if reaction.label else "reverse",
    )


def net_flows(network, abundances: Mapping[str, float], t9: float, rho: float = 1.0, include_partition: bool = True) -> Dict[str, Tuple[float, float, float]]:
    """Return ``{reaction: (forward, reverse, net)}`` fluxes for a network.

    The forward flux follows :meth:`Reaction.flux`; the reverse flux uses the
    detailed-balance rate with the product abundances.  At NSE abundances the
    net flux of every balanced reaction vanishes.
    """
    out: Dict[str, Tuple[float, float, float]] = {}
    species_map = network.species
    for r in network.reactions.reactions:
        fwd = r.flux(abundances, t9=t9, rho=rho)
        lam_r = reverse_rate(r, species_map, t9, rho=rho, include_partition=include_partition)
        s_rev, n_p = _side_stat_and_order(r.products, species_map)
        rev = lam_r * (float(rho) ** max(n_p - 1, 0)) / max(s_rev, 1)
        for p in r.products:
            if _nuclear(species_map, p.species) is None:
                continue
            y = max(float(abundances.get(p.species, 0.0)), 0.0)
            rev *= y ** p.count
        out[r.string] = (float(fwd), float(rev), float(fwd - rev))
    return out
