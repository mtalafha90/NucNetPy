from nucnetpy import Network, Zone, Species, Reaction, RateFit, evolve_zone, constant_thermo, time_grid

net = Network()
for name in ["h1", "he4"]:
    net.add_species(Species.parse(name))
net.reactions.add(Reaction.from_names(["h1", "h1", "h1", "h1"], ["he4"], rate_fits=[RateFit([0,0,0,0,0,0,0])], q_value=26.7))
z = Zone(properties={"t9": "1.0", "rho": "1e5"})
z.set_abundance("h1", 0.7)
z.set_abundance("he4", 0.0)
net.add_zone(z)
res = evolve_zone(net, z, time_grid(0, 1e-8, 20), constant_thermo(1.0, 1e5), method="rk4")
print(res.final_abundances)
