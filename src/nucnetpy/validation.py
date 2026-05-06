"""Validation and regression utilities for replacing NucNet Tools workflows."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence
import math
import numpy as np

from .core import Network, Zone
from .species import Species, normalize_species_name

@dataclass
class ValidationIssue:
    level: str
    code: str
    message: str


def validate_network(network: Network, strict: bool = False) -> List[ValidationIssue]:
    issues: List[ValidationIssue] = []
    for name, sp in network.species.items():
        if sp.a < 1 or sp.z < 0 or sp.z > sp.a:
            issues.append(ValidationIssue('error', 'bad_species', f'{name}: invalid Z/A'))
    missing = set()
    for r in network.reactions.reactions:
        for p in r.reactants + r.products:
            if p.species not in network.species:
                missing.add(p.species)
        ok, da, dz = r.conserves_a_z(network.species)
        if not ok:
            issues.append(ValidationIssue('error' if strict else 'warning', 'nonconserving_reaction', f'{r.string}: ΔA={da}, ΔZ={dz}'))
        if not (r.rate_fits or r.tabular_rate or r.constant_rate is not None):
            issues.append(ValidationIssue('warning', 'missing_rate', f'{r.string}: no rate data'))
    for m in sorted(missing):
        issues.append(ValidationIssue('warning', 'missing_species', f'{m}: appears in reactions but not in species list'))
    return issues


def validate_zone(zone: Zone, network: Optional[Network] = None, atol: float = 1e-6) -> List[ValidationIssue]:
    issues: List[ValidationIssue] = []
    if any(v < 0 for v in zone.abundances.values()):
        issues.append(ValidationIssue('error', 'negative_abundance', 'zone contains negative abundances'))
    try:
        xs = zone.xsum(network.species if network else None)
        if abs(xs - 1.0) > atol:
            issues.append(ValidationIssue('warning', 'mass_fraction_sum', f'sum X = {xs:.8g}, not 1'))
    except Exception as exc:
        issues.append(ValidationIssue('warning', 'mass_fraction_sum_failed', str(exc)))
    return issues


def compare_abundances(a: Mapping[str, float], b: Mapping[str, float], rtol: float = 1e-6, atol: float = 1e-30) -> Dict[str, float]:
    names = sorted(set(map(normalize_species_name, a)) | set(map(normalize_species_name, b)))
    out = {}
    for n in names:
        av = float(a.get(n, 0.0)); bv = float(b.get(n, 0.0))
        out[n] = abs(av - bv) / max(abs(bv), atol) if abs(av - bv) > atol else 0.0
    return out


def max_relative_abundance_error(a: Mapping[str, float], b: Mapping[str, float], atol: float = 1e-30) -> float:
    errs = compare_abundances(a, b, atol=atol)
    return max(errs.values()) if errs else 0.0


def regression_summary(reference: Mapping[str, float], candidate: Mapping[str, float]) -> Dict[str, object]:
    errs = compare_abundances(candidate, reference)
    worst = sorted(errs.items(), key=lambda kv: kv[1], reverse=True)[:10]
    return {"max_relative_error": max(errs.values()) if errs else 0.0, "worst_species": worst}
