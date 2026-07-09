"""A tour of the main nucnetpy features on a small helium-burning network.

Run with:  python examples/use_nucnetpy.py
"""
from nucnetpy import (
    Network, Zone, Species, Reaction, RateFit,
    evolve_zone, constant_thermo, time_grid,
    solve_nse, net_flows, system_timescales, entropy_generation_rate,
)

# --- Build a network: species need mass excesses (MeV) for NSE and
# --- detailed-balance reverse rates.
net = Network()
for name, mass_excess in [("he4", 2.4249), ("be8", 4.9416),
                          ("c12", 0.0), ("o16", -4.7370)]:
    net.add_species(Species.parse(name, mass_excess=mass_excess))

# Toy ReacLib fits: rate = exp(a0) at T9 = 1 (constant here for clarity).
net.reactions.add(Reaction.from_names(["he4", "he4"], ["be8"],
                                      rate_fits=[RateFit([2, 0, 0, 0, 0, 0, 0])],
                                      q_value=-0.0918))
net.reactions.add(Reaction.from_names(["be8", "he4"], ["c12"],
                                      rate_fits=[RateFit([3, 0, 0, 0, 0, 0, 0])],
                                      q_value=7.3666))
net.reactions.add(Reaction.from_names(["c12", "he4"], ["o16"],
                                      rate_fits=[RateFit([1, 0, 0, 0, 0, 0, 0])],
                                      q_value=7.1616))

# --- A zone of pure helium: Y(he4) = 0.25 means X(he4) = 1.
zone = Zone(abundances={"he4": 0.25})
net.add_zone(zone)
t9, rho = 2.0, 1.0e5

# --- Evolve one zone with the stiff BDF solver.
result = evolve_zone(net, zone, time_grid(0.0, 1.0e-2, 100),
                     thermo=constant_thermo(t9, rho), method="bdf")
print("Final abundances after 10 ms of helium burning:")
for name, y in sorted(result.final_abundances.items()):
    print(f"  {name:5s} Y = {y:.6e}   X = {net.species[name].a * y:.6e}")

# --- Forward / detailed-balance-reverse / net fluxes at the final state.
print("\nNet flows at the final state:")
for reaction, (fwd, rev, net_flux) in net_flows(net, result.final_abundances,
                                                t9=t9, rho=rho).items():
    print(f"  {reaction:25s} fwd={fwd:.3e}  rev={rev:.3e}  net={net_flux:.3e}")

# --- Shortest species timescales identify the stiff components.
print("\nSpecies timescales Y/|dY/dt| (s):")
for name, tau in sorted(system_timescales(net, 0, t9=t9, rho=rho).items(),
                        key=lambda kv: kv[1]):
    print(f"  {name:5s} {tau:.3e}")

# --- Entropy generation rate (k_B per nucleon per second, always >= 0).
print(f"\nEntropy generation rate: "
      f"{entropy_generation_rate(net, 0, t9=t9, rho=rho):.6e} k_B/nucleon/s")

# --- The NSE composition this network is driving toward.
nse = solve_nse(net, t9=4.0, rho=rho, ye=0.5)
print(f"\nNSE at T9=4, rho=1e5, Ye=0.5 (success={nse.success}):")
for name, y in sorted(nse.abundances.items(), key=lambda kv: -kv[1]):
    print(f"  {name:5s} Y = {y:.6e}")
