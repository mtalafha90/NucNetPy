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
from .reactions import PartitionRatio, Reaction, TabularRate
from .species import MASSLESS_SPECIES, Species, normalize_species_name

_LOG_MAX = 700.0

#: ReacLib writes the forward-to-reverse conversion with 1/(k in MeV/GK).
_REACLIB_INV_K = 11.6045

#: Participants that rule out a detailed-balance pair.  A photon may take part
#: (a capture and its photodisintegration); a lepton may not (a beta decay and
#: an electron capture are not each other's detailed-balance reverse).
_LEPTONS = frozenset(name for name in MASSLESS_SPECIES if name != "gamma")


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


def _q_value(reaction: Reaction, species_map: Mapping[str, Species],
             placeholders: Sequence[str] = ()) -> float:
    """Q of ``reaction`` in MeV from the nuclide masses, else its own record.

    libnucnet reaction files carry no Q value, so ``reaction.q_value`` is zero
    for every reaction read from one; the masses are the only reliable source.
    The record's value is used only when a nuclear participant has no mass in
    ``species_map`` or carries a placeholder mass.
    """
    total = 0.0
    for sign, side in ((1.0, reaction.reactants), (-1.0, reaction.products)):
        for p in side:
            if normalize_species_name(p.species) in MASSLESS_SPECIES:
                continue
            sp = species_map.get(normalize_species_name(p.species))
            if sp is None or sp.name in placeholders:
                return float(reaction.q_value)
            total += sign * p.count * float(sp.mass_excess)
    return total


def _is_library_reverse(candidate: Reaction, partner: Reaction,
                        rtol: float) -> bool:
    """True if every fit of ``candidate`` is ReacLib's reverse of a fit of ``partner``.

    ReacLib derives a reverse rate from its forward rate by detailed balance:
    the conversion keeps a2..a5, lowers a1 by 11.6045 Q (Q > 0 for the forward
    reaction), raises a6 by 3/2 times the change in the number of nuclear
    particles, and shifts a0 by the logarithm of the spin and mass factors.
    That signature identifies a derived reverse rate without the ``v`` flag,
    which libnucnet files do not carry.  An independently fitted reverse rate
    does not share a2..a5 with its forward rate and is left alone.
    """
    if not candidate.rate_fits or not partner.rate_fits:
        return False
    if len(candidate.rate_fits) != len(partner.rate_fits):
        return False
    if (candidate.tabular_rate is not None or candidate.constant_rate is not None
            or candidate.rate_function is not None):
        return False
    for r in (candidate, partner):
        if any(normalize_species_name(p.species) in _LEPTONS
               for p in r.reactants + r.products):
            return False

    def nuclear_count(side) -> int:
        return sum(p.count for p in side
                   if normalize_species_name(p.species) not in MASSLESS_SPECIES)

    a6_shift = 1.5 * (nuclear_count(partner.reactants)
                      - nuclear_count(partner.products))
    cand = [list(fit.coefficients) for fit in candidate.rate_fits]
    part = [list(fit.coefficients) for fit in partner.rate_fits]

    def compatible(c, f) -> bool:
        return (np.allclose(c[2:6], f[2:6], rtol=rtol, atol=1e-8)
                and abs((c[6] - f[6]) - a6_shift) <= 1e-4)

    # Every set shares one a1 shift, -11.6045 Q.  The sets need not be stored
    # in the same order in the two records, and sets with a2..a5 = 0 match each
    # other on those coefficients alone, so try each a1 shift the first set
    # allows and accept the one under which all the sets pair up.
    for first in part:
        if not compatible(cand[0], first):
            continue
        a1_shift = cand[0][1] - first[1]
        tol = 1e-4 * max(1.0, abs(a1_shift))
        unused = list(part)
        for c in cand:
            match = next((f for f in unused if compatible(c, f)
                          and abs((c[1] - f[1]) - a1_shift) <= tol), None)
            if match is None:
                break
            unused.remove(match)
        else:
            # The reverse of an exothermic forward rate is the member with the
            # lower a1.
            return a1_shift < -1e-6
    return False


def _partition_table(species_map: Mapping[str, Species], name: str):
    sp = species_map.get(normalize_species_name(name))
    if sp is None or not sp.partition:
        return (), ()
    keys = tuple(float(k) for k in sorted(sp.partition))
    return keys, tuple(float(sp.partition[k]) for k in keys)


def apply_reverse_partition_functions(network, rtol: float = 1e-6) -> list:
    r"""Multiply library reverse rates by the partition functions they omit.

    ReacLib supplies a reverse rate for most reactions, flagged ``v`` in its
    own format.  The library derives it from the forward rate by detailed
    balance using the ground-state spins only, so the temperature-dependent
    partition functions are left out and must be applied by the user: the
    reverse rate of ``A -> B`` is the fitted value times
    :math:`\prod_B G_i(T)/\prod_A G_i(T)`.  Evaluated without that factor, a
    network relaxes to an equilibrium in which no nuclide has excited states,
    which at :math:`T_9=5` differs from nuclear statistical equilibrium by up
    to 48 per cent in the alpha-chain mass fractions.

    libnucnet files do not carry the ``v`` flag, so each derived reverse rate
    is recognised from its coefficients (see :func:`_is_library_reverse`) and
    given a :class:`~nucnetpy.reactions.PartitionRatio` built from the
    partition functions in ``network.species``.  Species without a table count
    as :math:`G=1`.  The XML readers call this function by default; call it
    yourself for a network assembled in memory, or pass
    ``reverse_partition_functions=False`` to a reader to keep the uncorrected
    rates.

    Returns the list of reactions that were corrected.
    """
    reactions = network.reactions.reactions
    by_key = {r.key: r for r in reactions}
    corrected = []
    for reaction in reactions:
        partner = by_key.get((reaction.key[1], reaction.key[0]))
        if partner is None or partner is reaction:
            continue
        if not _is_library_reverse(reaction, partner, rtol):
            continue
        numerator = tuple(_partition_table(network.species, p.species) + (p.count,)
                          for p in reaction.products
                          if normalize_species_name(p.species) not in MASSLESS_SPECIES)
        denominator = tuple(_partition_table(network.species, p.species) + (p.count,)
                            for p in reaction.reactants
                            if normalize_species_name(p.species) not in MASSLESS_SPECIES)
        reaction.partition_ratio = PartitionRatio(numerator, denominator)
        corrected.append(reaction)
    return corrected


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


def _exact_reverse_reaction(forward: Reaction, species_map: Mapping[str, Species],
                            include_partition: bool) -> Reaction:
    """Reverse of ``forward`` whose rate is recomputed at each temperature."""
    def rate_at(t9: float, _f=forward, _m=species_map, _p=include_partition) -> float:
        return reverse_rate(_f, _m, float(t9), include_partition=_p)

    return Reaction(
        reactants=list(forward.products),
        products=list(forward.reactants),
        rate_function=rate_at,
        q_value=-forward.q_value,
        source="detailed_balance",
        label=(forward.label + "_reverse") if forward.label else "reverse",
    )


def consistent_reverse_network(network, t9_grid: Optional[Sequence[float]] = None,
                               include_partition: bool = True,
                               tabulate: bool = False):
    r"""Return a copy of ``network`` whose reverse rates obey detailed balance.

    Rate libraries such as JINA ReacLib derive each reverse rate from its
    forward rate by detailed balance, but with the library's own Q value, with
    ground-state spins only, and as a fitted expression.  Once the partition
    functions are restored (:func:`apply_reverse_partition_functions`, applied
    by the XML readers), the library reverse rates of the alpha chain at
    :math:`T_9=5` agree with the detailed-balance value implied by the nuclide
    file to about one per cent, which still leaves the stationary state of the
    network slightly away from the nuclear statistical equilibrium of those same
    data.

    This function pairs every reaction with its inverse, keeps the exothermic
    member of each pair, judged from the nuclide masses because libnucnet files
    carry no Q value, and regenerates the endothermic member from
    :func:`reverse_rate`.  Reactions with no inverse in the network are carried
    over unchanged.  The resulting network relaxes to the composition returned
    by :func:`nucnetpy.solve_nse` to about one part in :math:`10^{7}`, so the
    integrator and the equilibrium solver become mutually consistent.

    By default the reverse rate is attached as a callable evaluated exactly at
    whatever temperature the integrator asks for.  This matters: the reverse
    rate carries a factor :math:`\exp(-Q/kT)`, which is far too steep to
    interpolate accurately, and tabulating it on 200 logarithmically spaced
    points degrades the agreement with equilibrium from about one part in
    :math:`10^{7}` to about one part in :math:`10^{3}`.  Pass ``tabulate=True``
    to sample the rate on ``t9_grid`` instead, which is less accurate but can be
    written to XML.

    Note that this replaces the library's reverse rates, which carry the
    library's Q values and the error of its fits, with values derived from the
    masses in the nuclide file, so it enforces thermodynamic consistency with
    that file rather than improving individual rates.
    """
    import copy as _copy

    out = _copy.deepcopy(network)
    grid = (np.geomspace(0.01, 10.0, 200) if t9_grid is None
            else np.asarray(t9_grid, dtype=float))
    by_key = {r.key: r for r in out.reactions.reactions}
    placeholders = out.species_without_nuclear_data()

    rebuilt = []
    handled = set()
    for reaction in out.reactions.reactions:
        if reaction.key in handled:
            continue
        inverse = by_key.get((reaction.key[1], reaction.key[0]))
        if inverse is None:
            rebuilt.append(reaction)
            handled.add(reaction.key)
            continue
        handled.add(reaction.key)
        handled.add(inverse.key)
        forward = (reaction if _q_value(reaction, out.species, placeholders)
                   >= _q_value(inverse, out.species, placeholders) else inverse)
        rebuilt.append(forward)
        if tabulate:
            rebuilt.append(reverse_reaction(forward, out.species, t9_grid=grid,
                                            include_partition=include_partition))
        else:
            rebuilt.append(_exact_reverse_reaction(forward, out.species,
                                                   include_partition))

    out.reactions.reactions = rebuilt
    return out


def reverse_flux(reaction: Reaction, species_map: Mapping[str, Species],
                 abundances: Mapping[str, float], t9: float, rho: float = 1.0,
                 include_partition: bool = True) -> float:
    """Return the detailed-balance reverse flux of ``reaction`` at this composition.

    This is the flux the inverse reaction would carry with the rate from
    :func:`reverse_rate`, in the same units as :meth:`Reaction.flux`.
    """
    product = 1.0
    for p in reaction.products:
        if _nuclear(species_map, p.species) is None:
            continue
        product *= max(float(abundances.get(p.species, 0.0)), 0.0) ** p.count
    if product == 0.0:
        return 0.0          # avoids inf * 0 when the reverse rate overflows
    lam_r = reverse_rate(reaction, species_map, t9, rho=rho, include_partition=include_partition)
    s_rev, n_p = _side_stat_and_order(reaction.products, species_map)
    return float(lam_r * (float(rho) ** max(n_p - 1, 0)) / max(s_rev, 1) * product)


def net_flows(network, abundances: Mapping[str, float], t9: float, rho: float = 1.0, include_partition: bool = True) -> Dict[str, Tuple[float, float, float]]:
    """Return ``{reaction: (forward, reverse, net)}`` fluxes for a network.

    The forward flux follows :meth:`Reaction.flux`; the reverse flux uses the
    detailed-balance rate with the product abundances.  At NSE abundances the
    net flux of every balanced reaction vanishes.

    Every listed reaction is given its own detailed-balance reverse, including
    one whose inverse the network also lists explicitly, as JINA networks do.
    The entries for such a pair therefore describe the same process twice; do
    not add them up.  :func:`nucnetpy.analysis.entropy_generation_rate` and
    :func:`nucnetpy.analysis.integrated_currents` handle pairs correctly.
    """
    out: Dict[str, Tuple[float, float, float]] = {}
    species_map = network.species
    for r in network.reactions.reactions:
        fwd = r.flux(abundances, t9=t9, rho=rho)
        rev = reverse_flux(r, species_map, abundances, t9, rho=rho,
                           include_partition=include_partition)
        out[r.string] = (float(fwd), float(rev), float(fwd - rev))
    return out
