# Pitfalls

Everything in this chapter produces a plausible-looking answer rather than an
error message. Each item shows up with a production database rather than a
test network. NucNetPy guards against each one, and the test suite checks
those guards, but the traps belong to the problem, not to this implementation,
so they are worth knowing whichever code you use.

## Dropping the photon when you cut a network

A code that keeps a reaction only if *all* its participants survive a cut, and
cuts by a list of nuclides, discards every photodisintegration, because the
photon is not a nuclide. Nothing warns you. The network simply cannot reach
equilibrium, and a silicon-burning calculation will not start at all, because
it begins with $^{28}\mathrm{Si}(\gamma,\alpha)^{24}\mathrm{Mg}$.

`limit_network` avoids this: only nuclides decide whether a reaction survives,
so photons and leptons are kept automatically and listing `"gamma"` is
harmless. With `include_linked=True`, which adds every nuclide that shares a
reaction with the kept ones, the photon does not count as a link either;
otherwise every photon reaction would join the network:

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
median 4.5 per cent and up to 48 per cent away from NSE. No error message
shows it; only the comparison with `solve_nse` does.

The XML readers recognise the derived reverse rates from their coefficients,
because libnucnet files do not usually carry ReacLib's `v` flag, and apply the
ratio: the alpha chain then matches NSE to a median 0.34 per cent. For a
network built in memory, call `apply_reverse_partition_functions(net)` yourself.

The coefficients show which two rates form a pair, but not which of the two was
derived. The reader takes ReacLib's `v` flag from the fit labels when they carry
it. Otherwise it takes the photodisintegration of a capture/photodisintegration
pair, and the endothermic member of any other pair. Choosing the endothermic
member everywhere would put the factor on the wrong rate in 1,862 of the
32,676 pairs of the ReacLib snapshot of 29 March 2022. Those are
captures into nuclei beyond the drip lines, which are endothermic, yet ReacLib
still fits the capture and derives the photodisintegration. Eighteen two-body
pairs in that snapshot can be told apart only by the flag. Write the labels
with the flag (for example `ls09rv`) if you build a libnucnet file from ReacLib
yourself.

What is left comes from the library's Q values and the accuracy of its fits.
Use `consistent_reverse_network` if you need the network and the equilibrium
solver to agree exactly, and understand that you are then trusting the nuclide
file's masses over the library's Q values.

## Electron captures without the electron density

ReacLib fits its three electron captures (on $^3$He and $^7$Be, and the pep
reaction p + p + e⁻) per unit electron density, $\rho Y_e$. Treated like a
photodisintegration, whose rate already contains the photon density, they run
too slowly by that factor: by about 3 in the core of a massive star on the
main sequence, and by about 100 in the Sun. NucNetPy multiplies the ReacLib
fits of a reaction by $\rho Y_e$ once for each electron among its reactants.
`Reaction.rate` and `ReactionNetwork.rates` need `ye` for that. `Reaction.flux`,
`ReactionNetwork.flows`, `ReactionNetwork.ydot` and the network's right-hand
side work it out from the abundances.

## Running out of fuel with BDF

As a fuel runs out, the higher-order BDF formula can predict a slightly
negative abundance for it. The right-hand side treats a negative abundance as
zero, which switches off every reaction that consumes the fuel, and SciPy's
step size then collapses. On an 891-nuclide network from hydrogen to nickel this
stopped hydrogen burning soon after 1.8 of its 6.7 Myr. Neither removing the clipping
nor making the Jacobian agree with it cured the failure: each let the
integration go further and then fail in a later burning stage.

`method="adaptive_implicit_euler"` avoids the problem. It takes implicit Euler
steps, which cannot turn a decaying abundance negative, as NucNet Tools does,
and sizes each step so that no abundance changes by more than `max_change`
(0.5 per cent by default). It is first-order accurate: halving `max_change`
halves the error. Use it when a large network runs a fuel down.

## Integrating a whole database

`nucnetpy evolve-zone` integrates every species in the file it is given. For a
complete JINA database that is 8131 species and 81758 reactions, and one
evaluation of the right-hand side takes most of a second, because it is a
Python loop over the reactions. A stiff integration needs thousands of
evaluations, and far more when the conditions are explosive. A hydrogen-rich
zone at $T_9=2$ and $\rho=10^{6}\,\mathrm{g\,cm^{-3}}$ is such a case: proton
captures there take about $10^{-10}$ s, and on a 898-species network BDF had
covered less than $10^{-10}$ s of the requested second after more than two
minutes. The command therefore warns when a network has more than 400 species,
and reports how far it has got every 30 s on standard error (`--progress`).

Cut the network to the nuclides the problem needs with `--zmax` and `--amax`,
which keep the reactions among the retained nuclides, the weak decays and
photodisintegrations included. Do not cut below the nuclides the burning
reaches, or the material piles up at the edge of the network. `--atol 1e-12`
is often enough, and much faster than the default of $10^{-30}$, which asks
for every trace abundance to the relative tolerance. In Python, use
`limit_network` and the `atol` argument of `evolve_zone`.

## Placeholder nuclear data in an equilibrium solve

A production reaction file names nuclides that its nuclide file does not
describe — 279 of them in the JINA database behind the published results. NucNetPy
synthesises placeholders so the reactions remain usable, with a mass excess of
zero.

A mass excess of zero is not a neutral default: it makes an unbound nuclide look
as tightly bound as the most stable one. Left in an NSE solve, $^{5}$Li reached
a mass fraction of 0.28 and displaced the iron peak. They are therefore recorded
in `network.validate()["species_without_nuclear_data"]` and excluded from
equilibrium solves by default. `nucnetpy validate` reports them in a single
`no_nuclear_data` warning.

The record has to survive a file. Written out as an ordinary nuclide record
with a mass excess of zero, a placeholder could not be told from real data when
the file is read back, and an equilibrium solve on the combined file would
include it again. `nucnetpy jina-combine` and `write_xml` therefore leave
placeholders out of the nuclide records, and both readers recreate and record
them from the reactions that name them.

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
