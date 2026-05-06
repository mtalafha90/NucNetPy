"""Reaction and rate machinery for a pure-Python NucNet Tools replacement.

The implementation supports common ReacLib 7-coefficient rates, tabulated rates,
manual constants, duplicate handling, conservation checks, and stoichiometric
flux evaluation suitable for one-zone evolution and analysis commands.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from collections import Counter, defaultdict
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Any
import numpy as np

from .species import Species, normalize_species_name

@dataclass(frozen=True)
class ReactionParticipant:
    species: str
    count: int = 1

    def __post_init__(self):
        object.__setattr__(self, "species", normalize_species_name(self.species))
        object.__setattr__(self, "count", int(self.count))

@dataclass
class RateFit:
    """A standard seven-parameter ReacLib rate fit.

    Rate = exp(a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3) + a4*T9
               + a5*T9^(5/3) + a6*ln(T9))
    """
    coefficients: Sequence[float]
    label: str = ""
    chapter: Optional[int] = None

    def rate(self, t9: float) -> float:
        t9 = max(float(t9), 1e-30)
        a = list(self.coefficients)
        if len(a) != 7:
            raise ValueError("ReacLib RateFit requires exactly 7 coefficients")
        t13 = t9 ** (1.0 / 3.0)
        expo = a[0] + a[1] / t9 + a[2] / t13 + a[3] * t13 + a[4] * t9 + a[5] * (t9 ** (5.0 / 3.0)) + a[6] * math.log(t9)
        if expo > 700:
            return float("inf")
        if expo < -745:
            return 0.0
        return float(math.exp(expo))

@dataclass
class TabularRate:
    t9: Sequence[float]
    rates: Sequence[float]
    log_interp: bool = True

    def rate(self, t9: float) -> float:
        x = np.asarray(self.t9, dtype=float)
        y = np.asarray(self.rates, dtype=float)
        if len(x) == 0:
            return 0.0
        if len(x) == 1:
            return float(y[0])
        if self.log_interp:
            xx = np.log(np.clip(x, 1e-99, None))
            yy = np.log(np.clip(y, 1e-300, None))
            return float(np.exp(np.interp(math.log(max(t9, 1e-99)), xx, yy)))
        return float(np.interp(float(t9), x, y))

@dataclass
class Reaction:
    reactants: List[ReactionParticipant]
    products: List[ReactionParticipant]
    rate_fits: List[RateFit] = field(default_factory=list)
    tabular_rate: Optional[TabularRate] = None
    q_value: float = 0.0
    source: str = ""
    label: str = ""
    metadata: Dict[str, str] = field(default_factory=dict)
    constant_rate: Optional[float] = None

    @classmethod
    def from_names(cls, reactants: Sequence[str], products: Sequence[str], **kwargs) -> "Reaction":
        return cls(_participants_from_names(reactants), _participants_from_names(products), **kwargs)

    @property
    def key(self) -> Tuple[Tuple[Tuple[str, int], ...], Tuple[Tuple[str, int], ...]]:
        return (_side_key(self.reactants), _side_key(self.products))

    @property
    def string(self) -> str:
        return f"{_side_string(self.reactants)} -> {_side_string(self.products)}"

    @property
    def latex(self) -> str:
        return f"{_side_latex(self.reactants)} \\rightarrow {_side_latex(self.products)}"

    @property
    def reactant_order(self) -> int:
        return sum(p.count for p in self.reactants)

    def bare_rate(self, t9: float, rho: float = 1.0, ye: Optional[float] = None) -> float:
        total = 0.0
        if self.constant_rate is not None:
            total += float(self.constant_rate)
        if self.tabular_rate is not None:
            total += self.tabular_rate.rate(t9)
        for fit in self.rate_fits:
            total += fit.rate(t9)
        return float(total)

    def rate(self, t9: float, rho: float = 1.0, ye: Optional[float] = None, screening: Optional[Callable[["Reaction", float, float, Optional[float]], float]] = None) -> float:
        total = self.bare_rate(t9, rho=rho, ye=ye)
        if screening is not None:
            try:
                total *= float(screening(self, float(t9), float(rho), ye))
            except TypeError:
                total *= float(screening(self))
        return float(total)

    def statistical_factor(self) -> int:
        c = Counter()
        for p in self.reactants:
            c[p.species] += p.count
        out = 1
        for n in c.values():
            out *= math.factorial(n)
        return out

    def flux(self, abundances: Mapping[str, float], t9: float, rho: float = 1.0, screening: Optional[Callable[["Reaction", float, float, Optional[float]], float]] = None, ye: Optional[float] = None) -> float:
        lam = self.rate(t9, rho=rho, ye=ye, screening=screening)
        order = self.reactant_order
        flux = lam * (float(rho) ** max(order - 1, 0)) / max(self.statistical_factor(), 1)
        for p in self.reactants:
            y = max(float(abundances.get(p.species, 0.0)), 0.0)
            flux *= y ** p.count
        return float(flux)

    def stoichiometry(self) -> Dict[str, int]:
        d: Dict[str, int] = defaultdict(int)
        for p in self.products:
            d[p.species] += p.count
        for p in self.reactants:
            d[p.species] -= p.count
        return dict(d)

    def conserves_a_z(self, species_map: Mapping[str, Species]) -> Tuple[bool, int, int]:
        da = dz = 0
        for name, coeff in self.stoichiometry().items():
            sp = species_map.get(normalize_species_name(name))
            if sp is None:
                continue
            da += coeff * sp.a
            dz += coeff * sp.z
        return (da == 0 and dz == 0), da, dz

@dataclass
class ReactionNetwork:
    reactions: List[Reaction] = field(default_factory=list)

    def add(self, reaction: Reaction) -> None:
        self.reactions.append(reaction)

    def extend(self, reactions: Iterable[Reaction]) -> None:
        self.reactions.extend(reactions)

    def remove_duplicates(self) -> "ReactionNetwork":
        seen = set()
        unique = []
        for r in self.reactions:
            if r.key not in seen:
                unique.append(r); seen.add(r.key)
        self.reactions = unique
        return self

    def invalid_reactions(self, species_map: Mapping[str, Species]) -> List[Reaction]:
        return [r for r in self.reactions if not r.conserves_a_z(species_map)[0]]

    def filter_valid(self, species_map: Mapping[str, Species]) -> "ReactionNetwork":
        self.reactions = [r for r in self.reactions if r.conserves_a_z(species_map)[0]]
        return self

    def rates(self, t9: float, rho: float = 1.0) -> Dict[str, float]:
        return {r.string: r.rate(t9, rho=rho) for r in self.reactions}

    def flows(self, abundances: Mapping[str, float], t9: float, rho: float = 1.0, screening: Optional[Callable[[Reaction, float, float, Optional[float]], float]] = None, ye: Optional[float] = None) -> Dict[str, float]:
        return {r.string: r.flux(abundances, t9=t9, rho=rho, screening=screening, ye=ye) for r in self.reactions}

    def ydot(self, abundances: Mapping[str, float], t9: float, rho: float = 1.0, screening: Optional[Callable[[Reaction, float, float, Optional[float]], float]] = None, ye: Optional[float] = None) -> Dict[str, float]:
        out: Dict[str, float] = defaultdict(float)
        for r in self.reactions:
            f = r.flux(abundances, t9=t9, rho=rho, screening=screening, ye=ye)
            for name, nu in r.stoichiometry().items():
                out[name] += nu * f
        return dict(out)

    def species_names(self) -> List[str]:
        s = set()
        for r in self.reactions:
            for p in r.reactants + r.products:
                s.add(p.species)
        return sorted(s)

def _participants_from_names(names: Sequence[str]) -> List[ReactionParticipant]:
    c = Counter(normalize_species_name(n) for n in names if str(n).strip())
    return [ReactionParticipant(k, v) for k, v in sorted(c.items())]

def _side_key(parts: Sequence[ReactionParticipant]) -> Tuple[Tuple[str, int], ...]:
    c = Counter()
    for p in parts:
        c[p.species] += p.count
    return tuple(sorted(c.items()))

def _side_string(parts: Sequence[ReactionParticipant]) -> str:
    atoms = []
    for name, count in _side_key(parts):
        atoms.extend([name] * count)
    return " + ".join(atoms) if atoms else "∅"

def _side_latex(parts: Sequence[ReactionParticipant]) -> str:
    atoms = []
    for name, count in _side_key(parts):
        try:
            latex = Species.parse(name).latex
        except Exception:
            latex = name
        atoms.extend([latex] * count)
    return " + ".join(atoms) if atoms else r"\varnothing"


@dataclass
class DetailedBalanceOptions:
    """Options for estimating reverse rates from forward rates.

    Exact NucNet matching requires identical masses, partition functions, and
    Coulomb corrections.  This helper provides a practical replacement hook.
    """
    t9: float
    rho: float
    ye: float
    equilibrium_ratio: float = 1.0


def reverse_rate_from_detailed_balance(forward_rate: float, equilibrium_ratio: float) -> float:
    """Return reverse rate using ``lambda_rev = lambda_fwd / K_eq``."""
    return float(forward_rate) / max(float(equilibrium_ratio), 1e-300)


def read_reaclib_text(text: str) -> ReactionNetwork:
    """Read a small ReacLib-like text block.

    This parser accepts common two-line ReacLib records and a compact one-line
    form: ``reactants -> products ; a0 ... a6``.  It is intentionally permissive
    for migration from C++ examples and hand-written test networks.
    """
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith('#')]
    net = ReactionNetwork()
    i = 0
    while i < len(lines):
        line = lines[i]
        if '->' in line:
            lhs, rest = line.split('->', 1)
            rhs, coeff_text = (rest.split(';', 1) + [''])[:2] if ';' in rest else (rest, '')
            coeffs = []
            for tok in coeff_text.replace(',', ' ').split():
                try: coeffs.append(float(tok))
                except Exception: pass
            fits = [RateFit(coeffs)] if len(coeffs) == 7 else []
            net.add(Reaction(_participants_from_names(lhs.replace('+',' ').split()), _participants_from_names(rhs.replace('+',' ').split()), rate_fits=fits))
            i += 1
            continue
        # minimal two-line record fallback: species line then coefficients line
        toks = line.split()
        coeffs = []
        if i + 1 < len(lines):
            for tok in lines[i+1].replace(',', ' ').split():
                try: coeffs.append(float(tok))
                except Exception: pass
        if len(coeffs) >= 7 and len(toks) >= 2:
            mid = max(1, len(toks)//2)
            net.add(Reaction(_participants_from_names(toks[:mid]), _participants_from_names(toks[mid:]), rate_fits=[RateFit(coeffs[:7])]))
            i += 2
        else:
            i += 1
    return net
