"""Sparse/dense linear solver helpers replacing user/matrix_solver.* APIs."""
from __future__ import annotations
from typing import Any
import numpy as np


def solve_linear(a, b, method: str = "auto", **kwargs):
    method = method.lower()
    if method in {"spsolve", "auto", "ilu", "bicgstab", "gmres"}:
        try:
            from scipy.sparse import issparse
            from scipy.sparse.linalg import spsolve, gmres, bicgstab, spilu, LinearOperator
            if issparse(a):
                if method in {"auto", "spsolve"}:
                    return spsolve(a, b)
                if method == "gmres":
                    x, info = gmres(a, b, **kwargs)
                    if info != 0: raise RuntimeError(f"gmres failed with info={info}")
                    return x
                if method == "bicgstab":
                    x, info = bicgstab(a, b, **kwargs)
                    if info != 0: raise RuntimeError(f"bicgstab failed with info={info}")
                    return x
                if method == "ilu":
                    ilu = spilu(a.tocsc())
                    m = LinearOperator(a.shape, ilu.solve)
                    x, info = gmres(a, b, M=m, **kwargs)
                    if info != 0: raise RuntimeError(f"ilu-gmres failed with info={info}")
                    return x
        except ImportError:
            pass
    return np.linalg.solve(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
