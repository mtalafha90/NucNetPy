"""Analysis helpers equivalent to common NucNet Tools example programs."""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np

from .core import Network, Zone
from .species import Species, normalize_species_name


def largest_mass_fractions(zone: Zone, species_map: Mapping[str, Species], n: int = 10, min_x: float = 0.0):
    items = [(name, x) for name, x in zone.mass_fractions(species_map).items() if x >= min_x]
    return sorted(items, key=lambda kv: kv[1], reverse=True)[:n]


def abundances_for_element(zone: Zone, element: str, species_map: Mapping[str, Species]) -> Dict[str, float]:
    element = element.lower()
    return {name: y for name, y in zone.abundances.items() if species_map.get(name, Species.parse(name)).element == element}


def abundances_vs_nucleon_number(zone: Zone, species_map: Mapping[str, Species]) -> Dict[int, float]:
    out = defaultdict(float)
    for name, y in zone.abundances.items():
        out[species_map.get(name, Species.parse(name)).a] += y
    return dict(sorted(out.items()))


def element_abundances(zone: Zone, species_map: Mapping[str, Species]) -> Dict[str, float]:
    out = defaultdict(float)
    for name, y in zone.abundances.items():
        out[species_map.get(name, Species.parse(name)).element] += y
    return dict(sorted(out.items()))


def abundance_moment(zone: Zone, species_map: Mapping[str, Species], moment: int = 1, kind: str = "a") -> float:
    total = 0.0
    for name, y in zone.abundances.items():
        sp = species_map.get(name, Species.parse(name))
        q = sp.a if kind.lower() == "a" else sp.z if kind.lower() == "z" else sp.n
        total += (q ** moment) * y
    return float(total)


def species_history(network: Network, species: str) -> List[Tuple[int, Tuple[str, str, str], float]]:
    key = normalize_species_name(species)
    return [(i, z.label, z.get_abundance(key)) for i, z in enumerate(network.zones)]


def flows(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None):
    z = network.zone(zone_index)
    return network.reactions.flows(z.abundances, t9=t9 or z.temperature9(), rho=rho or z.density())


def ydot(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None):
    z = network.zone(zone_index)
    return network.reactions.ydot(z.abundances, t9=t9 or z.temperature9(), rho=rho or z.density())


def energy_generation_rate(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> float:
    # MeV per nucleon-ish proxy: sum Q * reaction flux. Users should calibrate
    # to their original libnucnet build for publication-grade energetics.
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    total = 0.0
    for r in network.reactions.reactions:
        total += r.q_value * r.flux(z.abundances, t9=t9, rho=rho)
    return float(total)


def compare_rates(net_a: Network, net_b: Network, t9: float, rho: float = 1.0) -> List[Tuple[str, float, float, float]]:
    ra = net_a.reactions.rates(t9, rho)
    rb = net_b.reactions.rates(t9, rho)
    keys = sorted(set(ra) | set(rb))
    return [(k, ra.get(k, 0.0), rb.get(k, 0.0), rb.get(k, 0.0) - ra.get(k, 0.0)) for k in keys]


def separation_energy(species_name: str, species_map: Mapping[str, Species], particle: str = "n") -> Optional[float]:
    """Return S_n or S_p in MeV, or None if the needed nuclides are absent.

    ``S_n(Z, A) = ME(Z, A-1) + ME(n) - ME(Z, A)`` and
    ``S_p(Z, A) = ME(Z-1, A-1) + ME(p) - ME(Z, A)``.
    """
    from .species import species_from_za
    sp = species_map.get(normalize_species_name(species_name), Species.parse(species_name))
    if sp.a <= 1:
        return None
    try:
        if particle.lower() == "n":
            daughter_name = species_from_za(sp.z, sp.a - 1).name
            particle_name = "n"
        else:
            if sp.z < 1:
                return None
            daughter_name = species_from_za(sp.z - 1, sp.a - 1).name
            particle_name = "h1"
    except ValueError:
        return None
    if daughter_name not in species_map or particle_name not in species_map:
        return None
    return species_map[daughter_name].mass_excess + species_map[particle_name].mass_excess - sp.mass_excess


def charge_changing_flows(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> Dict[str, float]:
    """Return each reaction's contribution to dYe/dt (flux times net ΔZ).

    Only reactions with a nonzero nuclear charge change (weak reactions such as
    beta decays and electron captures) appear in the result; strong reactions
    conserve Z and contribute nothing.  ``sum(result.values())`` is dYe/dt.
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    out: Dict[str, float] = {}
    for r in network.reactions.reactions:
        dz = 0
        for name, coeff in r.stoichiometry().items():
            sp = network.species.get(name)
            if sp is None:
                try:
                    sp = Species.parse(name)
                except Exception:
                    continue
            if sp.a > 0 and sp.z >= 0:
                dz += coeff * sp.z
        if dz != 0:
            out[r.string] = dz * r.flux(z.abundances, t9=t9, rho=rho)
    return out


def system_timescales(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None) -> Dict[str, float]:
    """Return per-species timescales ``Y / |dY/dt|`` in seconds.

    Species with zero derivative get ``inf``.  The shortest timescales identify
    the stiffest components of the system (blog: "Computing system timescales").
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    dy = network.reactions.ydot(z.abundances, t9=t9, rho=rho)
    out: Dict[str, float] = {}
    for name in set(z.abundances) | set(dy):
        y = z.get_abundance(name)
        rate = abs(dy.get(name, 0.0))
        out[name] = float(y / rate) if rate > 0.0 else float("inf")
    return out


def heavy_nuclei_abundance(zone: Zone, species_map: Mapping[str, Species], zmin: int = 3) -> float:
    """Return Y_h, the total abundance of nuclei with Z >= ``zmin``.

    The photon-to-heavy-nucleus ratio and the r-process neutron-to-seed ratio
    are built from this quantity (blog: "Computing the number of heavy nuclei").
    """
    total = 0.0
    for name, y in zone.abundances.items():
        sp = species_map.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        if sp.z >= zmin:
            total += y
    return float(total)


def neutron_exposure(result, thermo) -> float:
    """Return the s-process neutron exposure tau = ∫ n_n v_T dt in mb^-1.

    ``result`` is an :class:`~nucnetpy.solver.EvolutionResult` whose species
    include the neutron; ``thermo`` is the same ``(t, abundances) -> (t9, rho)``
    function used for the evolution.  ``n_n = rho N_A Y_n`` and the thermal
    velocity is ``v_T = sqrt(2 k T / m_n)`` (blog: "Computing the s-process
    neutron exposure").
    """
    from .constants import AVOGADRO, KB_CGS, MN_G
    if "n" not in result.species:
        return 0.0
    j = result.species.index("n")
    integrand = np.zeros(len(result.time))
    for i, t in enumerate(result.time):
        abund = {s: float(v) for s, v in zip(result.species, result.y[i])}
        t9, rho = thermo(float(t), abund)
        n_n = max(float(rho), 0.0) * AVOGADRO * max(float(result.y[i, j]), 0.0)
        v_t = np.sqrt(2.0 * KB_CGS * max(float(t9), 1e-30) * 1.0e9 / MN_G)
        integrand[i] = n_n * v_t
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 compat
    tau_cm2 = float(trapezoid(integrand, np.asarray(result.time, dtype=float)))
    return tau_cm2 * 1.0e-27  # cm^-2 -> mb^-1


#: Abundance at or below which a species counts as absent, as in NucNet Tools.
_ABSENT_ABUNDANCE = 1.0e-100


def _all_present(reaction, abundances: Mapping[str, float]) -> bool:
    """True if every nuclide the reaction touches has a non-negligible abundance."""
    from .species import is_massless
    return all(float(abundances.get(p.species, 0.0)) > _ABSENT_ABUNDANCE
               for p in reaction.reactants + reaction.products
               if not is_massless(p.species))


def _forward_fluxes(network: Network, abundances: Mapping[str, float],
                    t9: float, rho: float) -> Dict[tuple, float]:
    """Forward flux of each distinct reaction, keyed by ``Reaction.key``.

    Reactions listed more than once with the same reactants and products are
    added together.
    """
    from collections import defaultdict
    forward: Dict[tuple, float] = defaultdict(float)
    for r in network.reactions.reactions:
        forward[r.key] += r.flux(abundances, t9=t9, rho=rho)
    return dict(forward)


def _reverse_fluxes(network: Network, forward: Mapping[tuple, float],
                    abundances: Mapping[str, float], t9: float, rho: float) -> Dict[tuple, float]:
    """Reverse flux of each distinct reaction, keyed by ``Reaction.key``.

    This is the flux of the inverse reaction when the network lists it, which
    is what the calculation actually evolves, and the detailed-balance value
    otherwise.
    """
    from collections import defaultdict
    from .detailed_balance import reverse_flux
    reverse: Dict[tuple, float] = defaultdict(float)
    for r in network.reactions.reactions:
        inverse = (r.key[1], r.key[0])
        if inverse in forward and inverse != r.key:
            reverse[r.key] = forward[inverse]
        else:
            reverse[r.key] += reverse_flux(r, network.species, abundances, t9, rho=rho)
    return dict(reverse)


def reaction_entropy_changes(network: Network, abundances: Mapping[str, float], t9: float, rho: float) -> Dict[str, float]:
    """Return per-reaction entropy change ΔS in units of k_B per reaction.

    ``ΔS = -sum_i nu_i mu_i/kT`` with Maxwell–Boltzmann chemical potentials
    ``mu_i/kT = ln Y_i - ln pref_i``; this equals the C++ per-reaction
    ``Q/kT + sum_r ln(Y/Y_Q) - sum_p ln(Y/Y_Q)`` of NucNet Tools
    ``flow_utilities.cpp`` and, under detailed balance, ``ln(f/r)``.

    ΔS is formally infinite for a reaction that involves an absent species,
    because ``ln Y`` diverges as ``Y`` tends to zero.  Such an abundance is
    replaced by 1e-300 here, so the value returned for that reaction reflects
    the floor rather than the physics; :func:`entropy_generation_rate` leaves
    those reactions out.
    """
    from .nse import _log_prefactor
    import math
    out: Dict[str, float] = {}
    for r in network.reactions.reactions:
        ds = 0.0
        for name, nu in r.stoichiometry().items():
            sp = network.species.get(name)
            if sp is None:
                try:
                    sp = Species.parse(name)
                except Exception:
                    continue
            if sp.a <= 0:
                continue
            y = max(float(abundances.get(name, 0.0)), 1e-300)
            ds -= nu * (math.log(y) - _log_prefactor(sp, t9, rho))
        out[r.string] = float(ds)
    return out


def entropy_generation_rate(network: Network, zone_index: int = 0, t9: Optional[float] = None, rho: Optional[float] = None, use_reverse: bool = True) -> float:
    """Return dS/dt per nucleon in units of k_B per second.

    Ports NucNet Tools ``compute_entropy_generation_rate``: the sum over
    reactions of ``(f - r) * ΔS``, where ``ΔS`` is the per-reaction entropy
    change of :func:`reaction_entropy_changes` (blog series "Computing the
    entropy generation rate").  Each forward/reverse pair is counted once.
    When the network lists the inverse of a reaction, as JINA networks do,
    ``r`` is that inverse's own flux, so the result is exactly
    ``-sum_i (mu_i/kT) dY_i/dt`` for the network being evolved.  Otherwise
    ``r`` is the detailed-balance reverse flux, as in NucNet Tools.  With
    reverse rates that obey detailed balance every term equals
    ``(f - r) ln(f/r)``, so the total is non-negative and vanishes at NSE.

    With ``use_reverse=False`` only the listed reactions are used, each with
    its forward flux: ``sum_r f_r ΔS_r = -sum_i (mu_i/kT) dY_i/dt`` for exactly
    the network as given, with no detailed-balance reverses implied.

    Reactions that involve an absent species (``Y <= 1e-100``, the NucNet
    Tools threshold) are left out.  Creating a species from zero abundance
    produces formally infinite entropy, so any finite value for such a
    reaction would be set by an arbitrary floor rather than by the physics;
    treat the total for a composition with absent species as a lower bound.
    (NucNet Tools skips a reaction when its forward flux vanishes, and drops
    the logarithm of an absent product, which can give a term of the wrong
    sign.)

    Electron and neutrino chemical-potential terms of the C++ version are not
    included; nucnetpy networks carry weak reactions with their own tabulated
    rates instead.
    """
    z = network.zone(zone_index)
    t9 = t9 or z.temperature9(); rho = rho or z.density()
    abundances = z.abundances
    ds = reaction_entropy_changes(network, abundances, t9, rho)
    reactions = [r for r in network.reactions.reactions if _all_present(r, abundances)]
    if not use_reverse:
        return float(sum(r.flux(abundances, t9=t9, rho=rho) * ds[r.string] for r in reactions))
    forward = _forward_fluxes(network, abundances, t9, rho)
    reverse = _reverse_fluxes(network, forward, abundances, t9, rho)
    total = 0.0
    counted = set()
    for r in reactions:
        if r.key in counted:
            continue
        counted.update({r.key, (r.key[1], r.key[0])})
        total += (forward[r.key] - reverse[r.key]) * ds[r.string]
    return float(total)


def integrated_currents(network: Network, result, thermo, use_reverse: bool = True) -> Dict[str, float]:
    """Return per-reaction time-integrated net currents over an evolution.

    Ports NucNet Tools ``update_flow_currents``: for each reaction the current
    accumulates ``(forward - reverse) * dt`` across the time grid of
    ``result`` (an :class:`~nucnetpy.solver.EvolutionResult`), integrated here
    with the trapezoidal rule.  The integrated current of a reaction is the
    net number of transitions per nucleon it produced over the calculation
    (blog: "Creating integrated currents diagrams", "Analyzing integrated
    currents quantitatively").

    The reverse flux is the flux of the inverse reaction when the network
    lists it, and the detailed-balance value otherwise.  A pair listed in both
    directions therefore appears twice, with equal and opposite currents.
    With ``use_reverse=False`` each entry is the reaction's own forward
    current.  Reactions listed more than once are added together.
    """
    times = np.asarray(result.time, dtype=float)
    keys = list(dict.fromkeys(r.key for r in network.reactions.reactions))
    strings = {r.key: r.string for r in network.reactions.reactions}
    rates = np.zeros((len(times), len(keys)))
    for i, t in enumerate(times):
        abund = {s: float(v) for s, v in zip(result.species, result.y[i])}
        t9, rho = thermo(float(t), abund)
        forward = _forward_fluxes(network, abund, t9, rho)
        if use_reverse:
            reverse = _reverse_fluxes(network, forward, abund, t9, rho)
            rates[i] = [forward[k] - reverse[k] for k in keys]
        else:
            rates[i] = [forward[k] for k in keys]
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2.0 compat
    return {strings[k]: float(trapezoid(rates[:, j], times)) for j, k in enumerate(keys)}


def nuclear_energy_generation_rate(network: Network, zone_index: int = 0,
                                   t9: Optional[float] = None,
                                   rho: Optional[float] = None,
                                   screening=None,
                                   abundances: Optional[Mapping[str, float]] = None,
                                   ydot_values: Optional[Mapping[str, float]] = None) -> float:
    """Return the nuclear energy generation rate in erg per gram per second.

    The rate follows from the change in total mass excess of the composition,

        eps = -N_A * C * sum_i (dY_i/dt) * Delta M_i,

    with ``Delta M_i`` the mass excess in MeV and
    ``C = 1.602176634e-6 erg/MeV`` (``constants.MEV_TO_ERG``) converting the
    sum, which carries units of MeV/g/s, into erg/g/s.  This is exact, and it is
    preferable to summing ``Q_r`` times the flow of each reaction: it does not
    use the reaction Q-values stored in the rate library, so it is internally
    consistent with the adopted mass-excess table -- the same table the
    equilibrium solver uses -- and energy release and equilibrium are therefore
    built from one set of nuclear data.  It may differ from library Q-values
    when the two data sources use different nuclear masses.

    A positive value denotes energy released to the plasma.  Neutrino losses
    are not subtracted: energy carried away by neutrinos from weak reactions
    must be accounted for separately.

    ``abundances`` and ``ydot_values`` override the zone composition and its
    derivative, which is useful when the rate is wanted along a trajectory that
    has already been integrated.
    """
    from .constants import AVOGADRO, MEV_TO_ERG

    zone = network.zone(zone_index)
    t9 = t9 or zone.temperature9()
    rho = rho or zone.density()
    if abundances is None:
        abundances = zone.abundances
    if ydot_values is None:
        ydot_values = network.reactions.ydot(abundances, t9=t9, rho=rho,
                                             screening=screening)

    total = 0.0
    for name, dy in ydot_values.items():
        sp = network.species.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        if sp.a <= 0:
            continue
        total -= float(dy) * float(sp.mass_excess)
    return float(total * AVOGADRO * MEV_TO_ERG)


def nuclear_energy_release(network: Network,
                           initial: Mapping[str, float],
                           final: Mapping[str, float]) -> float:
    """Return the energy released between two compositions, in erg per gram.

    This is the time integral of :func:`nuclear_energy_generation_rate` and
    depends only on the endpoints, being the difference in total mass excess.
    """
    from .constants import AVOGADRO, MEV_TO_ERG

    total = 0.0
    for name in set(initial) | set(final):
        sp = network.species.get(name)
        if sp is None:
            try:
                sp = Species.parse(name)
            except Exception:
                continue
        if sp.a <= 0:
            continue
        dy = float(final.get(name, 0.0)) - float(initial.get(name, 0.0))
        total -= dy * float(sp.mass_excess)
    return float(total * AVOGADRO * MEV_TO_ERG)
