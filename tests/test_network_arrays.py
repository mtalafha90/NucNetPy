"""The array right-hand side and Jacobian against the per-reaction calculation.

``rhs`` and ``analytic_jacobian`` evaluate every reaction at once, as arrays.
Each test here compares them with the plain per-reaction calculation, kept
below as a reference, on a network that holds every kind of rate and reactant
the solver has to handle.  The two agree to rounding error; on most platforms
they agree bit for bit.
"""
import numpy as np
import pytest

from nucnetpy import SkyNetScreening, consistent_reverse_network
from nucnetpy.core import Network, Zone
from nucnetpy.io.xml import read_xml
from nucnetpy.reactions import PartitionRatio, RateFit, Reaction, TabularRate
from nucnetpy.solver import (analytic_jacobian, constant_thermo, evolve_zone,
                             jacobian_sparsity, rhs, time_grid)
from nucnetpy.species import Species, is_massless
from nucnetpy.weak import read_weak_table

GOLDEN = "tests/golden/"


# --- the per-reaction reference --------------------------------------------

def _ye(network, species, y):
    total = 0.0
    for s, value in zip(species, y):
        sp = network.species.get(s) or Species.parse(s)
        total += sp.z * max(float(value), 0.0)
    return total


def reference_rhs(network, species, thermo, screening=None, weak_rates=None):
    idx = {s: i for i, s in enumerate(species)}

    def f(t, yvec):
        y = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = {s: float(y[i]) for s, i in idx.items()}
        t9, rho = thermo(t, abund)
        ye = _ye(network, species, y)
        dy = network.reactions.ydot(abund, t9=t9, rho=rho, screening=screening, ye=ye)
        out = np.zeros_like(y)
        for name, value in dy.items():
            if name in idx:
                out[idx[name]] += value
        for wr in weak_rates or ():
            if wr.parent in idx and wr.daughter in idx:
                flow = wr.rate(t9, rho * max(ye, 1e-30)) * abund[wr.parent]
                out[idx[wr.parent]] -= flow
                out[idx[wr.daughter]] += flow
        return out
    return f


def reference_jacobian(network, species, thermo, screening=None, weak_rates=None):
    idx = {s: i for i, s in enumerate(species)}

    def j(t, yvec):
        y = np.clip(np.asarray(yvec, dtype=float), 0.0, np.inf)
        abund = {s: float(y[i]) for s, i in idx.items()}
        t9, rho = thermo(t, abund)
        ye = _ye(network, species, y)
        J = np.zeros((len(y), len(y)))
        for r in network.reactions.reactions:
            nuclear = r.nuclear_reactants
            reactants = [(idx[p.species], p.count) for p in nuclear if p.species in idx]
            if len(reactants) != len(nuclear):
                continue
            lam = r.rate(t9, rho=rho, ye=ye, screening=screening)
            pref = lam * rho ** max(r.reactant_order - 1, 0) / max(r.statistical_factor(), 1)
            for pos, (jx, mj) in enumerate(reactants):
                term = pref * mj * y[jx] ** (mj - 1)
                for other, (kx, mk) in enumerate(reactants):
                    if other != pos:
                        term *= y[kx] ** mk
                for name, nu in r.stoichiometry().items():
                    if name in idx:
                        J[idx[name], jx] += nu * term
        for wr in weak_rates or ():
            if wr.parent in idx and wr.daughter in idx:
                rate = wr.rate(t9, rho * max(ye, 1e-30))
                J[idx[wr.parent], idx[wr.parent]] -= rate
                J[idx[wr.daughter], idx[wr.parent]] += rate
        return J
    return j


# --- a network with every kind of rate --------------------------------------

def _fit(a0, a1=0.0, a6=0.0):
    return RateFit([a0, a1, 0.0, 0.0, 0.0, 0.0, a6])


def mixed_network():
    net = Network()
    for name, z, a, dm, spin in [("n", 0, 1, 8.071, 0.5), ("h1", 1, 1, 7.289, 0.5),
                                 ("h2", 1, 2, 13.136, 1.0), ("he3", 2, 3, 14.931, 0.5),
                                 ("he4", 2, 4, 2.425, 0.0), ("li7", 3, 7, 14.907, 1.5),
                                 ("be7", 4, 7, 15.769, 1.5), ("c12", 6, 12, 0.0, 0.0),
                                 ("n13", 7, 13, 5.346, 0.5), ("o16", 8, 16, -4.737, 0.0)]:
        net.add_species(Species(name, z, a, mass_excess=dm, spin=spin))
    add = net.reactions.add
    # Two fit sets, two identical reactants.
    add(Reaction.from_names(["h1", "h1"], ["h2", "positron", "neutrino_e"],
                            rate_fits=[_fit(-34.8, -3.5, -0.6), _fit(-30.0, -4.0)]))
    # Electron captures: one with two identical nuclear reactants.
    add(Reaction.from_names(["h1", "h1", "electron"], ["h2", "neutrino_e"],
                            rate_fits=[_fit(-43.6, -0.002, -0.2)]))
    add(Reaction.from_names(["be7", "electron"], ["li7", "neutrino_e"],
                            rate_fits=[_fit(-23.8, a6=-0.5)]))
    # Three identical reactants.
    add(Reaction.from_names(["he4", "he4", "he4"], ["c12"], rate_fits=[_fit(-11.0, -4.4)]))
    # A capture and its reverse, with a partition-function factor.
    add(Reaction.from_names(["c12", "h1"], ["n13"], rate_fits=[_fit(17.0, -3.0, -0.7)]))
    reverse = Reaction.from_names(["n13"], ["c12", "h1"], rate_fits=[_fit(40.0, -22.0, 0.8)])
    reverse.partition_ratio = PartitionRatio(
        numerator=(((0.1, 1.0, 10.0), (1.0, 1.02, 1.3), 1),),
        denominator=(((0.1, 1.0, 10.0), (1.0, 1.05, 1.6), 1),))
    add(reverse)
    # A photodisintegration: the photon is not an abundance.
    add(Reaction.from_names(["o16", "gamma"], ["c12", "he4"], rate_fits=[_fit(50.0, -83.0)]))
    # Tabulated, constant and computed rates, one of them with a fit as well.
    add(Reaction.from_names(["he3", "he4"], ["be7"],
                            tabular_rate=TabularRate([0.1, 1.0, 10.0], [1e-8, 1e-2, 1e2])))
    add(Reaction.from_names(["n"], ["h1", "electron", "anti-neutrino_e"], constant_rate=1.0 / 880.0))
    add(Reaction.from_names(["h2", "h1"], ["he3"], rate_fits=[_fit(5.0, -3.7)],
                            rate_function=lambda t9: 1e-3 * t9 ** 2))
    # A reaction with a reactant that a partial species list leaves out.
    add(Reaction.from_names(["li7", "h1"], ["he4", "he4"], rate_fits=[_fit(20.0, -8.5)]))
    zone = Zone(abundances={"n": 1e-6, "h1": 0.7, "h2": 2e-5, "he3": 1e-5, "he4": 0.07,
                            "li7": 1e-9, "be7": 1e-10, "c12": 2e-4, "n13": 1e-8, "o16": 5e-4})
    net.add_zone(zone)
    return net


def species_of(net):
    return [s for s in net.species_names() if not is_massless(s)]


def states(net, species, count=4):
    """The zone's abundances, then perturbed ones with some negative entries."""
    rng = np.random.default_rng(3)
    y0 = np.array([net.zone(0).get_abundance(s) for s in species])
    out = [y0]
    for _ in range(count - 1):
        y = y0 + rng.random(len(y0)) * 10.0 ** rng.uniform(-10, -2, len(y0))
        y[rng.random(len(y)) < 0.2] *= -1.0
        out.append(y)
    return out


class UpdatingScreening:
    """A screening model that, like SkyNetScreening, is refreshed with the state."""

    def __init__(self):
        self.ye = None

    def update(self, abundances, t9, rho, ye):
        self.ye = ye

    def __call__(self, reaction, t9, rho, ye):
        return 1.0 + 0.01 * len(reaction.nuclear_reactants) * (self.ye or 0.0)


CASES = {
    "plain": dict(),
    "screening": dict(screening=lambda r, t9, rho, ye: 1.0 + 0.1 * len(r.nuclear_reactants)),
    "updating screening": dict(screening=UpdatingScreening()),
}


# --- the tests ---------------------------------------------------------------

@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("partial", [False, True])
def test_rhs_matches_the_per_reaction_calculation(case, partial):
    net = mixed_network()
    species = species_of(net)
    if partial:
        species = [s for s in species if s != "li7"]
    thermo = constant_thermo(t9=1.5, rho=1.0e5)
    kwargs = CASES[case]
    new, old = rhs(net, species, thermo, **kwargs), reference_rhs(net, species, thermo, **kwargs)
    for y in states(net, species):
        np.testing.assert_allclose(new(0.0, y), old(0.0, y), rtol=1e-14, atol=0.0)


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("partial", [False, True])
def test_jacobian_matches_the_per_reaction_calculation(case, partial):
    net = mixed_network()
    species = species_of(net)
    if partial:
        species = [s for s in species if s != "li7"]
    thermo = constant_thermo(t9=1.5, rho=1.0e5)
    kwargs = CASES[case]
    sparse = analytic_jacobian(net, species, thermo, **kwargs)
    dense = analytic_jacobian(net, species, thermo, sparse=False, **kwargs)
    old = reference_jacobian(net, species, thermo, **kwargs)
    for y in states(net, species):
        expected = old(0.0, y)
        assert sparse(0.0, y).format == "csc"
        np.testing.assert_allclose(sparse(0.0, y).toarray(), expected, rtol=1e-14, atol=0.0)
        np.testing.assert_allclose(dense(0.0, y), expected, rtol=1e-14, atol=0.0)
        # No stored zeros: the structure is that of the nonzero entries.
        assert np.array_equal(sparse(0.0, y).toarray() != 0.0, expected != 0.0)


def test_weak_tables_and_screening_on_the_golden_network():
    net = read_xml(GOLDEN + "golden_network.xml")
    species = species_of(net)
    weak = [read_weak_table(GOLDEN + "golden_weak_rates.txt")]
    thermo = constant_thermo(t9=3.0, rho=1.0e7)
    for kwargs in (dict(weak_rates=weak), dict(screening=SkyNetScreening(net.species), weak_rates=weak)):
        new, old = rhs(net, species, thermo, **kwargs), reference_rhs(net, species, thermo, **kwargs)
        jn, jo = analytic_jacobian(net, species, thermo, **kwargs), reference_jacobian(net, species, thermo, **kwargs)
        for y in states(net, species):
            np.testing.assert_allclose(new(0.0, y), old(0.0, y), rtol=1e-14, atol=0.0)
            np.testing.assert_allclose(jn(0.0, y).toarray(), jo(0.0, y), rtol=1e-14, atol=0.0)


def test_reverse_rates_from_detailed_balance():
    net = consistent_reverse_network(read_xml(GOLDEN + "golden_network.xml"))
    species = species_of(net)
    thermo = constant_thermo(t9=3.0, rho=1.0e7)
    new, old = rhs(net, species, thermo), reference_rhs(net, species, thermo)
    for y in states(net, species):
        np.testing.assert_allclose(new(0.0, y), old(0.0, y), rtol=1e-14, atol=0.0)


def test_rates_follow_a_changing_temperature():
    # The temperature-dependent part of each rate is cached; it must be
    # recomputed when T9 changes and reused when it comes back.
    net = mixed_network()
    species = species_of(net)
    def thermo(t, abundances):
        return 1.0 + t, 1.0e5 * (1.0 + t)
    new, old = rhs(net, species, thermo), reference_rhs(net, species, thermo)
    y = states(net, species)[1]
    for t in (0.0, 0.5, 2.0, 0.5, 0.0):
        np.testing.assert_allclose(new(t, y), old(t, y), rtol=1e-14, atol=0.0)


def test_rhs_conserves_baryon_number():
    net = mixed_network()
    species = species_of(net)
    mass = np.array([net.species[s].a for s in species], dtype=float)
    f = rhs(net, species, constant_thermo(t9=1.5, rho=1.0e5))
    for y in states(net, species):
        dy = f(0.0, y)
        assert abs(mass @ dy) <= 1e-12 * np.sum(np.abs(mass * dy))


def test_sparsity_pattern_matches_the_reactions():
    net = mixed_network()
    species = species_of(net)
    idx = {s: i for i, s in enumerate(species)}
    expected = np.zeros((len(species), len(species)), dtype=bool)
    for r in net.reactions.reactions:
        for name in r.stoichiometry():
            for p in r.reactants:
                if name in idx and p.species in idx:
                    expected[idx[name], idx[p.species]] = True
    pattern = jacobian_sparsity(net, species)
    assert pattern.dtype == bool and pattern.format == "csc"
    assert np.array_equal(pattern.toarray(), expected)


def test_evolution_is_unchanged_by_the_array_form():
    # The same integration with the reference right-hand side and Jacobian,
    # patched in, must give the same trajectory.
    import nucnetpy.solver as solver
    net = mixed_network()
    zone = net.zone(0)
    times = time_grid(0.0, 1.0e3, 12)
    thermo = constant_thermo(t9=0.3, rho=1.0e4)
    fast = evolve_zone(net, zone, times, thermo=thermo, atol=1e-20)
    saved = solver.rhs, solver.analytic_jacobian
    try:
        solver.rhs = lambda n, s, th, screening=None, weak_rates=None: reference_rhs(n, s, th, screening, weak_rates)
        from scipy.sparse import csc_matrix
        solver.analytic_jacobian = (lambda n, s, th, screening=None, weak_rates=None, sparse=True:
                                    (lambda f: (lambda t, y: csc_matrix(f(t, y))))(reference_jacobian(n, s, th, screening, weak_rates)))
        slow = evolve_zone(net, zone, times, thermo=thermo, atol=1e-20)
    finally:
        solver.rhs, solver.analytic_jacobian = saved
    assert fast.success and slow.success
    assert (fast.nfev, fast.njev) == (slow.nfev, slow.njev)
    np.testing.assert_allclose(fast.y, slow.y, rtol=1e-9, atol=1e-25)
