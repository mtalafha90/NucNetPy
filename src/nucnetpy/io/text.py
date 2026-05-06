from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, List
from ..core import Zone
from ..reactions import Reaction, RateFit


def read_abundance_text(path) -> Zone:
    z = Zone()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) >= 2:
                z.set_abundance(parts[0], float(parts[1]))
    return z


def read_simple_reactions(path) -> List[Reaction]:
    rxns = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            if ';' in line:
                spec, coeff = line.split(';', 1)
            else:
                spec, coeff = line, ''
            left, right = spec.split('->')
            coeffs = [float(x) for x in coeff.split()] if coeff.strip() else []
            fits = [RateFit(coeffs)] if len(coeffs) == 7 else []
            rxns.append(Reaction.from_names(left.replace('+',' ').split(), right.replace('+',' ').split(), rate_fits=fits))
    return rxns
