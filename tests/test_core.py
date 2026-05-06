from nucnetpy import Species, Reaction, RateFit, Network, Zone, evolve_zone, constant_thermo, time_grid


def test_species_parse():
    assert Species.parse('ni56').z == 28
    assert Species.parse('alpha').name == 'he4'


def test_reaclib_rate():
    r = RateFit([0,0,0,0,0,0,0])
    assert abs(r.rate(1.0) - 1.0) < 1e-12


def test_reaction_flux_and_ydot():
    rxn = Reaction.from_names(['h1','h1'], ['he2'], rate_fits=[RateFit([0,0,0,0,0,0,0])])
    f = rxn.flux({'h1': 1.0}, t9=1.0, rho=2.0)
    assert f == 1.0  # rate*rho/factorial(2)
    assert rxn.stoichiometry()['h1'] == -2


def test_evolve_zone():
    net = Network()
    for n in ['h1','he2']:
        net.add_species(Species.parse(n))
    net.reactions.add(Reaction.from_names(['h1','h1'], ['he2'], rate_fits=[RateFit([0,0,0,0,0,0,0])]))
    z = Zone(); z.set_abundance('h1', 1.0); z.set_abundance('he2', 0.0); net.add_zone(z)
    res = evolve_zone(net, z, time_grid(0, 1e-4, 5), constant_thermo(1,1), method='rk4')
    assert res.final_abundances['h1'] < 1.0
