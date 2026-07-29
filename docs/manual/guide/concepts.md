# Concepts

Five objects carry everything.

## Species

A nuclide: name, charge `z`, mass number `a`, mass excess, ground-state spin,
and an optional partition-function table.

```python
import nucnetpy as nn

ni56 = nn.Species("ni56", z=28, a=56, mass_excess=-53.907, spin=0.0)
```

Mass excess and spin are not decoration. Equilibrium, detailed balance and
energy release are all computed from mass excesses, and the statistical weight
in an equilibrium solve is `(2J+1) * G(T)`, so a placeholder value produces a
placeholder answer. A species whose spin is unknown is given a weight of one
rather than a guess.

`Species.parse` builds one from a name alone, normalising aliases
(`p` → `h1`, `alpha` → `he4`), but it cannot invent nuclear data.

### Participants that are not nuclides

`gamma`, `electron`, `positron` and the neutrinos appear in reaction records so
that charge and lepton number balance, but they carry no baryon number and are
never part of the abundance vector. `nucnetpy.species.is_massless` identifies
them. See {doc}`pitfalls`.

## Reaction

Reactants, products, a rate, a `Q`-value, and metadata. The rate may be a
ReacLib seven-coefficient fit, a tabulated rate, a constant, or a callable; a
reaction may carry more than one, and they add.

```python
r = nn.Reaction.from_names(["c12", "he4"], ["o16", "gamma"],
                           rate_fits=[nn.RateFit([10, 0, 0, 0, 0, 0, 0])],
                           q_value=7.162)
```

The flow of a reaction is

$$F_r = \frac{\lambda_r\,\rho^{\max(n_r-1,0)}}{\prod_i m_{ir}!}\prod_i Y_i^{m_{ir}},$$

where the products run over the reactants that carry baryon number, and $n_r$
is that reactant order. `Reaction.nuclear_reactants` is the list actually used.

## ReactionNetwork

A collection of reactions, with duplicate removal, filtering, and the
construction of rates, flows, and $dY_i/dt$.

## Zone

One composition together with its thermodynamic properties. Abundances are
molar: `Y` is moles per gram, and the mass fraction is `X = A*Y`. A normalised
composition has `sum(X) = 1`.

```python
zone = nn.Zone(abundances={"si28": 1.0 / 28.0},
               properties={"t9": "5.0", "rho": "1e8"})
```

## Network

Species, reactions and zones in one container, plus `validate()`, which reports
species named by a reaction but absent from the species map, reactions that do
not balance, and species carrying no nuclear data.

## Units

| Quantity | Unit |
|---|---|
| Temperature `t9` | $10^9$ K |
| Density `rho` | g cm⁻³ |
| Abundance `Y` | mol g⁻¹ |
| Mass excess | MeV |
| `Q`-value | MeV |
| Rate | s⁻¹ scaled by the density powers of the flow |
| Energy generation | erg g⁻¹ s⁻¹ |
| Time | s |
