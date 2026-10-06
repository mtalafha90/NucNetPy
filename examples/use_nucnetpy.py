"""A tour of the main nucnetpy features on a small helium-burning network.

The network burns helium to carbon and oxygen.  Its reverse rates are rebuilt
from detailed balance, so the burn settles on exactly the composition that the
equilibrium solver predicts, which this script shows.

Run with:  python examples/use_nucnetpy.py
"""
import copy

from nucnetpy import (
    Network, Zone, Species, Reaction,
    evolve_zone, constant_thermo, time_grid, solve_nse,
    consistent_reverse_network, net_flows, system_timescales, entropy_generation_rate,
)

# --- Species.  Mass excesses (MeV) are needed for equilibrium and for
# --- reverse rates derived from detailed balance.
net = Network()
for name, mass_excess in [("he4", 2.4249), ("be8", 4.9416),
                          ("c12", 0.0), ("o16", -4.7370)]:
    net.add_species(Species.parse(name, mass_excess=mass_excess, spin=0.0))

# --- Reactions.  A rate library lists each reaction in both directions with
# --- separately fitted rates.  The constant rates here are toy values.
for reactants, products, rate, q in [(["be8"], ["he4", "he4"], 1e2, 0.0918),
                                     (["be8", "he4"], ["c12"], 1e9, 7.3666),
                                     (["c12", "he4"], ["o16"], 1e4, 7.1616)]:
    net.reactions.add(Reaction.from_names(reactants, products, constant_rate=rate, q_value=q))
    net.reactions.add(Reaction.from_names(products, reactants, constant_rate=rate * 1e-3, q_value=-q))

# Library reverse rates do not match the nuclear masses exactly, so the
# network would settle away from equilibrium.  Rebuild them from detailed
# balance so that the burn and the equilibrium solver agree.
net = consistent_reverse_network(net)

# --- Burn pure helium (Y = 0.25, so X = 1) at constant temperature and density.
t9, rho = 3.5, 1.0e6
zone = Zone(abundances={"he4": 0.25})
net.add_zone(zone)
result = evolve_zone(net, zone, time_grid(1e-8, 1e6, 40, log=True),
                     thermo=constant_thermo(t9, rho), method="bdf", rtol=1e-8, atol=1e-23)
assert result.success, result.message

# --- The burn ends at the NSE composition for the same T9, rho and Ye, to a
# --- few parts per million (the limit of SciPy's BDF Newton iteration).
nse = solve_nse(net, t9=t9, rho=rho, ye=0.5)
print(f"After {result.time[-1]:.0e} s at T9 = {t9}, rho = {rho:.0e} g/cm^3:")
print("  species   X (burn)       X (NSE)")
for name in ["he4", "be8", "c12", "o16"]:
    a = net.species[name].a
    print(f"  {name:5s}   {a * result.final_abundances[name]:.6e}   {a * nse.abundances[name]:.6e}")

# --- Every net flow vanishes at equilibrium.  The network lists each reaction
# --- in both directions, so show each pair once, by its exothermic member.
print("\nNet flows at the end (forward minus reverse):")
flows_at_end = net_flows(net, result.final_abundances, t9=t9, rho=rho)
for reaction in net.reactions.reactions:
    if reaction.q_value > 0:
        fwd, rev, net_flux = flows_at_end[reaction.string]
        print(f"  {reaction.string:22s} {net_flux / fwd:+.1e} of the forward flow")

# --- Diagnostics part-way through the burn, where every species is present.
middle = copy.deepcopy(net)
middle.zones = [Zone(abundances=dict(zip(result.species, result.y[12])))]
print(f"\nPart-way through the burn (t = {result.time[12]:.1e} s):")
print("  species timescales Y/|dY/dt| (s):")
for name, tau in sorted(system_timescales(middle, 0, t9=t9, rho=rho).items(), key=lambda kv: kv[1]):
    print(f"    {name:5s} {tau:.3e}")
print(f"  entropy generation rate: {entropy_generation_rate(middle, 0, t9=t9, rho=rho):.4e} "
      f"k_B per nucleon per second")
end = copy.deepcopy(net)
end.zones = [Zone(abundances=result.final_abundances)]
print(f"At the end the entropy generation rate is "
      f"{entropy_generation_rate(end, 0, t9=t9, rho=rho):.1e}: equilibrium.")
