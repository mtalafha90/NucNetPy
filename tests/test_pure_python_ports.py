import numpy as np
from nucnetpy.mathutils import linear_interpolation, bilinear_interpolation, two_d_interpolation
from nucnetpy.decay import decay_constant_from_half_life, DecayRecord
from nucnetpy.hydro import exponential_expansion
from nucnetpy.matrix_solver import solve_linear
from nucnetpy.network_limiter import limit_network
from nucnetpy import Network, Species, Reaction, ReactionNetwork, Zone


def test_mathutils_interpolation():
    assert linear_interpolation([0,1], [10,20], -1) == 10
    assert linear_interpolation([0,1], [10,20], 2) == 20
    assert linear_interpolation([0,1], [10,20], 0.5) == 15
    val, diff = bilinear_interpolation([0,1], [0,1], [[0,1],[2,3]], 0.5, 0.5)
    assert val == 1.5 and diff == 3
    val2, diff2 = two_d_interpolation([0,1], [0,1], [[0,1],[2,3]], -1, 0.5)
    assert diff2 == 9999.0


def test_decay_hydro_matrix_and_limiter():
    assert decay_constant_from_half_life(2.0) > 0
    rxn = DecayRecord('ni56', 'co56', 6.0).reaction()
    assert rxn.constant_rate > 0
    tr = exponential_expansion(0, 1.0, 10.0, tau=2.0, n=5)
    assert tr.thermo(0)[0] == 1.0
    x = solve_linear(np.eye(2), np.array([1.0, 2.0]))
    assert np.allclose(x, [1,2])
    net = Network()
    for n in ['h1','he2','he4']:
        net.add_species(Species.parse(n))
    net.reactions = ReactionNetwork([Reaction.from_names(['h1','h1'], ['he2']), Reaction.from_names(['he2','he2'], ['he4'])])
    net.add_zone(Zone(abundances={'h1':1,'he2':0,'he4':0}))
    limit_network(net, ['h1','he2'])
    assert 'he4' not in net.species
    assert len(net.reactions.reactions) == 1
