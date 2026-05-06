"""Weak-rate and neutrino helper layer.

Supports two-dimensional T9/rhoYe tables, simple FFN/Oda/Langanke-style text
blocks after normalization, bilinear interpolation, and ODE source-term hooks.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np
from .species import normalize_species_name

@dataclass
class WeakRateTable:
    parent: str
    daughter: str
    t9: Sequence[float]
    rho_ye: Sequence[float]
    rates: np.ndarray
    nu_loss: Optional[np.ndarray] = None
    label: str = ""

    def __post_init__(self):
        self.parent = normalize_species_name(self.parent)
        self.daughter = normalize_species_name(self.daughter)
        self.t9 = np.asarray(self.t9, dtype=float)
        self.rho_ye = np.asarray(self.rho_ye, dtype=float)
        self.rates = np.asarray(self.rates, dtype=float)
        if self.rates.shape != (len(self.t9), len(self.rho_ye)):
            self.rates = np.reshape(self.rates, (len(self.t9), len(self.rho_ye)))
        if self.nu_loss is not None:
            self.nu_loss = np.asarray(self.nu_loss, dtype=float)

    def _interp2(self, values: np.ndarray, t9: float, rho_ye: float, log_values: bool = True) -> float:
        x = np.log(np.clip(self.t9, 1e-99, None))
        y = np.log(np.clip(self.rho_ye, 1e-99, None))
        z = np.log(np.clip(values, 1e-300, None)) if log_values else np.asarray(values, dtype=float)
        lx = np.log(max(t9, 1e-99)); ly = np.log(max(rho_ye, 1e-99))
        if len(x) == 1 and len(y) == 1:
            val = z[0, 0]
        elif len(x) == 1:
            iy = np.searchsorted(y, ly) - 1; iy = max(0, min(iy, len(y)-2))
            ty = (ly-y[iy]) / max(y[iy+1]-y[iy], 1e-99)
            val = (1-ty)*z[0,iy] + ty*z[0,iy+1]
        elif len(y) == 1:
            ix = np.searchsorted(x, lx) - 1; ix = max(0, min(ix, len(x)-2))
            tx = (lx-x[ix]) / max(x[ix+1]-x[ix], 1e-99)
            val = (1-tx)*z[ix,0] + tx*z[ix+1,0]
        else:
            ix = np.searchsorted(x, lx) - 1; iy = np.searchsorted(y, ly) - 1
            ix = max(0, min(ix, len(x)-2)); iy = max(0, min(iy, len(y)-2))
            tx = (lx - x[ix]) / max(x[ix+1] - x[ix], 1e-99)
            ty = (ly - y[iy]) / max(y[iy+1] - y[iy], 1e-99)
            val = (1-tx)*(1-ty)*z[ix,iy] + tx*(1-ty)*z[ix+1,iy] + (1-tx)*ty*z[ix,iy+1] + tx*ty*z[ix+1,iy+1]
        return float(np.exp(val) if log_values else val)

    def rate(self, t9: float, rho_ye: float) -> float:
        return self._interp2(self.rates, t9, rho_ye, log_values=True)

    def neutrino_loss(self, t9: float, rho_ye: float) -> float:
        if self.nu_loss is None:
            return 0.0
        return self._interp2(self.nu_loss, t9, rho_ye, log_values=False)


def read_weak_table(path: str | Path, parent: Optional[str] = None, daughter: Optional[str] = None) -> WeakRateTable:
    """Read a normalized weak-rate table.

    Accepted format::

        # parent fe56 daughter mn56
        t9 rhoYe rate [nu_loss]
        1.0 1e7 2.3e-4 1.0e-5

    Repeated rows are reshaped into a rectangular T9 x rhoYe grid.
    """
    path = Path(path)
    rows = []
    meta: Dict[str, str] = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            toks = s[1:].split()
            for i, tok in enumerate(toks[:-1]):
                if tok.lower() in {'parent','daughter','label'}:
                    meta[tok.lower()] = toks[i+1]
            continue
        toks = s.replace(',', ' ').split()
        vals = []
        for tok in toks[:4]:
            try: vals.append(float(tok))
            except Exception: pass
        if len(vals) >= 3:
            rows.append(vals)
    if not rows:
        raise ValueError(f"No weak-rate rows found in {path}")
    parent = parent or meta.get('parent')
    daughter = daughter or meta.get('daughter')
    if not parent or not daughter:
        raise ValueError("parent and daughter must be supplied or present in a header")
    arr = np.asarray(rows, dtype=float)
    t9 = np.unique(arr[:,0]); rho = np.unique(arr[:,1])
    rates = np.zeros((len(t9), len(rho)))
    nuloss = np.zeros_like(rates) if arr.shape[1] >= 4 else None
    ti = {v:i for i,v in enumerate(t9)}; ri = {v:i for i,v in enumerate(rho)}
    for row in arr:
        i = ti[row[0]]; j = ri[row[1]]
        rates[i,j] = row[2]
        if nuloss is not None and len(row) >= 4:
            nuloss[i,j] = row[3]
    return WeakRateTable(parent, daughter, t9, rho, rates, nuloss, label=meta.get('label',''))


def read_weak_tables(paths: Iterable[str | Path]) -> List[WeakRateTable]:
    return [read_weak_table(p) for p in paths]


def compute_yedot(abundances: Mapping[str, float], weak_rates: Sequence[WeakRateTable], species_map, t9: float, rho: float, ye: float) -> float:
    total = 0.0
    for wr in weak_rates:
        parent = species_map.get(wr.parent)
        daughter = species_map.get(wr.daughter)
        if not parent or not daughter:
            continue
        dz = daughter.z - parent.z
        total += dz * abundances.get(wr.parent, 0.0) * wr.rate(t9, rho*ye)
    return float(total)


def weak_source_terms(abundances: Mapping[str, float], weak_rates: Sequence[WeakRateTable], t9: float, rho: float, ye: float) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for wr in weak_rates:
        flow = wr.rate(t9, rho * max(ye, 1e-30)) * abundances.get(wr.parent, 0.0)
        out[wr.parent] = out.get(wr.parent, 0.0) - flow
        out[wr.daughter] = out.get(wr.daughter, 0.0) + flow
    return out
