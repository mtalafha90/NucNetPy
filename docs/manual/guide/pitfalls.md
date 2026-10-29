# Pitfalls

Everything in this chapter produces a plausible-looking answer rather than an
error message. Each item was found by running the package against a production
database rather than a test network, and each is now covered by a regression
test — but the underlying trap belongs to the problem, not to this
implementation, so it is worth knowing whichever code you use.

## Dropping the photon when you cut a network

A code that keeps a reaction only if *all* its participants survive a cut, and
cuts by a list of nuclides, discards every photodisintegration, because the
photon is not a nuclide. Nothing warns you. The network simply cannot reach
equilibrium, and a silicon-burning calculation will not start at all, because
it begins with $^{28}\mathrm{Si}(\gamma,\alpha)^{24}\mathrm{Mg}$.

`limit_network` avoids this: only nuclides decide whether a reaction survives,
so photons and leptons are kept automatically and listing `"gamma"` is
harmless. Development versions of NucNetPy before the 1.0.0 release fell into
the trap too, and with `include_linked=True` they also treated the photon as a
link between otherwise unrelated reactions:

```python
limit_network(net, select_species(net, zmax=20))   # keeps the photodisintegrations
```

## Using library reverse rates without their partition functions

A network relaxes to the equilibrium implied by the **ratios of the rates it is
given**, not to the equilibrium of the nuclear data. ReacLib derives each
reverse rate from its forward rate by detailed balance with the ground-state
spins only, and expects the user to multiply it by the ratio of partition
functions, $\prod G_{\rm products}/\prod G_{\rm reactants}$ of the reverse
reaction.

Left out, that factor is not small. On the alpha chain at $T_9=5$ the reverse
fits then differ from detailed balance by a median factor of 1.10 and up to 1.72
(for $^{56}$Ni$(\gamma,\alpha)$), and the network relaxes to a composition a
median 4.5 per cent and up to 48 per cent away from NSE. Versions of NucNetPy
before the fix did exactly that; the comparison with `solve_nse` is what
exposed it.

The XML readers now recognise the derived reverse rates from their coefficients,
because libnucnet files do not carry ReacLib's `v` flag, and apply the ratio:
the alpha chain then matches NSE to a median 0.34 per cent. For a network built
in memory, call `apply_reverse_partition_functions(net)` yourself.

What is left comes from the library's Q values and the accuracy of its fits.
Use `consistent_reverse_network` if you need the network and the equilibrium
solver to agree exactly, and understand that you are then trusting the nuclide
file's masses over the library's Q values.

## Placeholder nuclear data in an equilibrium solve

A production reaction file names nuclides that its nuclide file does not
describe — 279 of them in the database used for the manuscript. NucNetPy
synthesises placeholders so the reactions remain usable, with a mass excess of
zero.

A mass excess of zero is not a neutral default: it makes an unbound nuclide look
as tightly bound as the most stable one. Left in an NSE solve, $^{5}$Li reached
a mass fraction of 0.28 and displaced the iron peak. They are therefore recorded
in `network.validate()["species_without_nuclear_data"]` and excluded from
equilibrium solves by default.

## Trusting `success` without checking it

```python
result = nn.evolve_zone(...)
assert result.success, result.message
```

`success` is false when the solver failed, and also when positivity projection
had to invent enough baryon number to break conservation by more than the
requested tolerance.

That second case is worth understanding. The right-hand side clips its input to
non-negative values, so once a component goes negative the derivative stops
depending on it and its Jacobian column vanishes; an implicit solver is then
free to carry that component far negative while still passing its own
convergence test. Clipping it back to zero at the end converts a large negative
number into invented mass. In one long integration this returned
$\sum_i X_i = 8.56$ with $X(^4\mathrm{He}) = 7.18$, while SciPy reported that it
had reached the end of the interval successfully.

Checking conservation costs one inner product per output step. Do it.

## Forgetting the spin degeneracy

libnucnet partition-function tables are *normalised*: they tend to one as
$T \rightarrow 0$ for every nuclide, high-spin ones included. The statistical
weight is $(2J+1)\,G(T)$, and the degeneracy has to be supplied separately.
Omitting it makes every nuclide with $J \neq 0$ too rare by that factor — up to
90 per cent for $^{50}$V.

NucNetPy applies it automatically when a species carries a `spin`. If your
species were built with `Species.parse`, they do not.

## Tolerances

The default `atol` is an absolute floor on abundances. Set it far below the
smallest abundance you care about, but not so far that the solver chases
numerical noise in species that are effectively absent — that is a common cause
of a stiff solve failing to start at all. Confirm the answer has stopped moving:

```python
for rtol in [1e-4, 1e-6, 1e-8, 1e-10]:
    r = nn.evolve_zone(..., rtol=rtol, atol=rtol * 1e-8)
```

## Interpolating a detailed-balance reverse rate

The reverse rate carries $\exp(-Q/kT)$, which is far too steep to interpolate.
Tabulating it on 200 logarithmically spaced points costs three orders of
magnitude of accuracy against evaluating it exactly. `consistent_reverse_network`
attaches a callable by default for this reason; `tabulate=True` is for when the
network must be serialisable.
