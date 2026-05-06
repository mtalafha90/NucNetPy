"""Hydrodynamic trajectory utilities replacing user/hydro* helpers."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, Tuple
import numpy as np
from .mathutils import linear_interpolation

@dataclass
class Trajectory:
    time: np.ndarray
    t9: np.ndarray
    rho: np.ndarray

    def thermo(self, t: float, y=None) -> Tuple[float, float]:
        return (linear_interpolation(self.time, self.t9, t), linear_interpolation(self.time, self.rho, t))

    @classmethod
    def from_columns(cls, time: Sequence[float], t9: Sequence[float], rho: Sequence[float]) -> "Trajectory":
        return cls(np.asarray(time, dtype=float), np.asarray(t9, dtype=float), np.asarray(rho, dtype=float))


def read_trajectory(path: str | Path) -> Trajectory:
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        vals = [float(x) for x in s.replace(',', ' ').split()[:3]]
        if len(vals) == 3:
            rows.append(vals)
    if not rows:
        raise ValueError(f"No trajectory rows found in {path}")
    arr = np.asarray(rows, dtype=float)
    return Trajectory(arr[:,0], arr[:,1], arr[:,2])


def exponential_expansion(t0: float, t9_0: float, rho_0: float, tau: float, n: int = 200, t_end: float | None = None) -> Trajectory:
    t_end = float(t_end if t_end is not None else 10.0 * tau)
    time = np.linspace(t0, t_end, int(n))
    fac = np.exp(-(time - t0) / max(float(tau), 1e-99))
    return Trajectory(time, t9_0 * fac, rho_0 * fac)
