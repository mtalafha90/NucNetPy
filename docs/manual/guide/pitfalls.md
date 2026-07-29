# Pitfalls

Everything in this chapter produces a plausible-looking answer rather than an
error message. Each item was found by running the package against a production
database rather than a test network, and each is now covered by a regression
test — but the underlying trap belongs to the problem, not to this
implementation, so it is worth knowing whichever code you use.

## Dropping the photon when you cut a network

`limit_network` retains only reactions whose participants *all* survive the
cut, and `select_species` returns nuclides. The photon is not a nuclide, so:

```python
limit_network(net, select_species(net, zmax=20))            # wrong
limit_network(net, select_species(net, zmax=20) + ["gamma"]) # right
```

Without it, every photodisintegration is discarded. Nothing warns you. The
network simply cannot reach equilibrium, and a silicon-burning calculation will
not start at all, because it begins with
$^{28}\mathrm{Si}(\gamma,\alpha)^{24}\mathrm{Mg}$.

## Assuming reverse rates are consistent with the masses

A network relaxes to the equilibrium implied by the **ratios of the rates it is
given**, not to the equilibrium of the nuclear data. Rate libraries fit the two
directions separately, so the two need not agree.

Evaluated at the NSE composition, JINA's reverse fits leave the net flux of a
reaction far from zero; rebuilt from detailed balance, the net flux vanishes to
machine precision. On the alpha chain at $T_9=5$ the fits differ from detailed
balance by a median factor of 1.10 and up to 1.72, which displaces equilibrium
mass fractions by a median of 4.5 per cent.

Use `consistent_reverse_network` if you need the network and the equilibrium
solver to agree, and understand that you are then trusting the mass table over
the library's reverse fit.

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
