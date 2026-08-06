# Common workflows

## Reading a JINA/libnucnet database

```python
import nucnetpy as nn

net = nn.read_jina_xml("nuclides.xml", "reaction_data.xml", zones_xml="zone.xml")
report = net.validate()
print(len(net.species), "species in memory,",
      len(net.reactions.reactions), "reactions")
print("without nuclear data:", len(report["species_without_nuclear_data"]))
```

Count carefully here. `len(net.species)` is the size of the in-memory network
after parsing, not the number of entries in the nuclide file. Anything named by
a reaction record with no entry in that file gets a synthesised placeholder, so
it is present in memory whether or not the file described it.

On the JINA database used throughout this manual the nuclide file holds 7852
nuclides, while `len(net.species)` is 8131: those 7852 plus 279 placeholders.
`validate()` reports no missing species precisely because the placeholders
already exist, so `species_without_nuclear_data` is the field to read, not
`missing_species`.

The two containers also differ. `net.species` holds nuclides only, while
`net.species_names()` returns 8136 — the same 8131 plus the massless reaction
participants `gamma`, `electron`, `positron` and the two neutrinos, which
balance a reaction record but never enter the abundance vector.

Reading the database takes 20 to 30 s and holds it in about 0.9 GB, for a
process peak near 1 GB. That is the practical ceiling for a pure-Python
implementation, which is why the next step matters.

## Cutting a network down to size

A one-zone calculation rarely needs the whole chart of nuclides.

```python
from nucnetpy.network_limiter import select_species, limit_network

limit_network(net, select_species(net, zmax=20, amax=44) + ["gamma"])
```

**Keep `gamma`.** `limit_network` retains only reactions whose participants all
survive the cut, and `select_species` returns nuclides, so dropping the photon
silently deletes every photodisintegration. See {doc}`pitfalls`.

## Evolving one zone

```python
result = nn.evolve_zone(
    net, net.zone(0), nn.time_grid(0.0, 10.0, 120),
    thermo=nn.constant_thermo(t9=5.0, rho=1.0e8),
    method="bdf", rtol=1e-8, atol=1e-14,
)
assert result.success, result.message
```

`thermo` is any callable `(t, abundances) -> (T9, rho)`, so a prescribed
temperature and density history drops straight in;
`nucnetpy.hydro.Trajectory` wraps a tabulated one with interpolation.

Supported methods are `bdf`, `radau`, `lsoda`, `rk45` and `dop853`, plus
fixed-step `rk4`, `euler` and `implicit_euler` fallbacks. Reaction networks are
stiff, so the first three are the practical choices.

### Choosing the Jacobian

The flows are monomials in the abundances, so the Jacobian is available in
closed form and costs one sweep over the reaction list, against `N+1`
right-hand-side evaluations for finite differences.

```python
nn.evolve_zone(..., jac_mode="analytic")   # default; ~275x faster at 100 species
nn.evolve_zone(..., jac_mode="numerical")  # finite differences
nn.evolve_zone(..., jac_mode="sparsity")   # SciPy estimates it from the pattern
```

SciPy applies a sparsity pattern only when no Jacobian callable is supplied, so
the pattern and an explicit Jacobian are alternatives rather than complements.

## Equilibrium

```python
nse = nn.solve_nse(net, t9=5.0, rho=1.0e8, ye=0.5)
print(nse.success, nse.xsum, nse.computed_ye)
```

The solve shares no machinery with the time integration, which is what makes it
a genuine check on a network calculation. Species carrying no nuclear data are
excluded automatically; pass `require_nuclear_data=False` to override that, and
read {doc}`pitfalls` before doing so.

Constrained cluster equilibria (the libnuceq workflow) use `solve_qse` with
`QSECluster` constraints. Plasma corrections are available through
`nucnetpy.coulomb.nse_correction`.

## Thermodynamic consistency

Rate libraries fit the forward and reverse directions of a reaction separately,
so a network built from them relaxes to a state that is not the NSE of the same
nuclear data.

```python
consistent = nn.consistent_reverse_network(net)
```

This pairs each reaction with its inverse, keeps the exothermic member, and
regenerates the other from detailed balance. On the JINA alpha chain it improves
agreement between the integrated composition and an independent NSE solve from a
median of 4.5 per cent to 4.5 parts per million.

It is an option rather than a default, because it replaces measured reverse
rates with values derived from mass differences. Pass `tabulate=True` if the
network must remain writable to XML, at some cost in accuracy.

## Screening

```python
from nucnetpy import SkyNetScreening

nn.evolve_zone(..., screening=SkyNetScreening(net.species))
```

Each nuclear charge is assigned a Coulomb chemical potential `mu(Z)/kT`, and a
reaction is enhanced by the difference between the separated reactants and the
fused compound charge. Being a chemical potential, it applies consistently to
forward and reverse rates, so screening does not undo detailed balance. A
pairwise Salpeter factor is also available through `reaction_screening_factor`.

## Energy release

```python
from nucnetpy import nuclear_energy_generation_rate, nuclear_energy_release
```

The rate follows from the change in total mass excess,
`eps = -N_A * sum_i (dY_i/dt) * dM_i`, so it needs no `Q`-values and cannot
disagree with them where a rate library and a nuclide file are inconsistent.
Neutrino losses are not subtracted.

## Analysis

`nucnetpy.analysis` carries the diagnostics of the NucNet Tools workflows:
largest mass fractions, element abundances, abundance moments, per-species
timescales `Y/|dY/dt|`, per-reaction contributions to `dYe/dt`, entropy
generation, `s`-process neutron exposure, integrated reaction currents, and
separation energies.

## Command line

Every operation above has a command-line equivalent; see {doc}`../reference/cli`.
