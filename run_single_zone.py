from nucnetpy import read_jina_xml
from nucnetpy.solver import evolve_zone, time_grid
from nucnetpy.thermo import constant_thermo


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------
nuclides_xml = "nuclides.xml"
reactions_xml = "reaction_data.xml"
zone_xml = "zone(1).xml"


# ------------------------------------------------------------
# 2. Load JINA nuclides + reactions + initial zone
# ------------------------------------------------------------
net = read_jina_xml(
    nuclides_xml=nuclides_xml,
    reactions_xml=reactions_xml,
    zones_xml=zone_xml,
)

zone = net.zone(0)

print("Number of species:", len(net.species))
print("Number of reactions:", len(net.reactions.reactions))
print("Initial zone properties:")
for k, v in zone.properties.items():
    print(f"  {k} = {v}")


# ------------------------------------------------------------
# 3. Get thermodynamic parameters from the zone file
# ------------------------------------------------------------
t9_0 = float(zone.properties.get("t9_0", 0.20))
rho_0 = float(zone.properties.get("rho_0", 1.5e4))
tend = float(zone.properties.get("tend", 100.0))

# For a first test, use a simple fixed-temperature/fixed-density run.
thermo = constant_thermo(
    t9=t9_0,
    rho=rho_0,
)


# ------------------------------------------------------------
# 4. Build time grid
# ------------------------------------------------------------
times = time_grid(
    t0=0.0,
    t1=tend,
    n=200,
)


# ------------------------------------------------------------
# 5. Evolve the single zone
# ------------------------------------------------------------
result = evolve_zone(
    network=net,
    zone=zone,
    times=times,
    thermo=thermo,
    method="bdf",      # stiff solver
    rtol=1e-8,
    atol=1e-20,
)


# ------------------------------------------------------------
# 6. Print final abundances
# ------------------------------------------------------------
print("\nFinal abundances:")
for name, y in sorted(result.final_abundances.items()):
    if y > 1e-20:
        print(f"{name:8s}  Y = {y:.8e}")


# ------------------------------------------------------------
# 7. Print final mass fractions
# ------------------------------------------------------------
print("\nFinal mass fractions:")
for name, x in sorted(result.final_mass_fractions.items()):
    if x > 1e-20:
        print(f"{name:8s}  X = {x:.8e}")