"""Numerical utility routines ported from NucNet Tools helper math code."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence, Tuple
import numpy as np

OUT_OF_TABLE_DIFF = 9999.0


def linear_interpolation(x: Sequence[float], y: Sequence[float], x0: float) -> float:
    """C++ compatible clamped 1-D interpolation.

    Matches the behavior of nnt::linear_interpolation: values below the first
    grid point return y[0], values at/above the last point return y[-1].
    """
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    if xv.ndim != 1 or yv.ndim != 1 or len(xv) != len(yv) or len(xv) == 0:
        raise ValueError("x and y must be non-empty one-dimensional arrays of equal length")
    if x0 < xv[0]:
        return float(yv[0])
    if x0 >= xv[-1]:
        return float(yv[-1])
    j = int(np.searchsorted(xv, x0, side="right") - 1)
    return float(yv[j] + (yv[j+1] - yv[j]) * (x0 - xv[j]) / (xv[j+1] - xv[j]))


def bilinear_interpolation(x1: Sequence[float], x2: Sequence[float], matrix, v1: float, v2: float) -> Tuple[float, float]:
    """Strict in-table bilinear interpolation.

    Returns ``(value, max_corner_difference)`` and raises if outside the table,
    following the C++ routine that exits on out-of-range input.
    """
    x1 = np.asarray(x1, dtype=float); x2 = np.asarray(x2, dtype=float)
    m = np.asarray(matrix, dtype=float)
    if not (x1[0] <= v1 <= x1[-1]) or not (x2[0] <= v2 <= x2[-1]):
        raise ValueError("interpolation point outside table")
    j = min(max(int(np.searchsorted(x1, v1, side="right") - 1), 0), len(x1)-2)
    k = min(max(int(np.searchsorted(x2, v2, side="right") - 1), 0), len(x2)-2)
    y1 = m[j, k]; y2 = m[j+1, k]; y3 = m[j+1, k+1]; y4 = m[j, k+1]
    t = (v1 - x1[j]) / (x1[j+1] - x1[j])
    u = (v2 - x2[k]) / (x2[k+1] - x2[k])
    val = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4
    diff = max(abs(y1-y2), abs(y1-y3), abs(y1-y4), abs(y2-y3), abs(y2-y4), abs(y3-y4))
    return float(val), float(diff)


def two_d_interpolation(x1: Sequence[float], x2: Sequence[float], matrix, v1: float, v2: float) -> Tuple[float, float]:
    """2-D interpolation with boundary clamping.

    This mirrors the public behavior of the NucNet helper: interpolate inside
    the table, extrapolate/clamp on edges, and report ``9999`` as uncertainty
    when the requested point lies outside the available grid.
    """
    x1 = np.asarray(x1, dtype=float); x2 = np.asarray(x2, dtype=float)
    m = np.asarray(matrix, dtype=float)
    outside = v1 < x1[0] or v1 > x1[-1] or v2 < x2[0] or v2 > x2[-1]
    vv1 = float(np.clip(v1, x1[0], x1[-1])); vv2 = float(np.clip(v2, x2[0], x2[-1]))
    if len(x1) == 1 and len(x2) == 1:
        return float(m[0,0]), OUT_OF_TABLE_DIFF if outside else 0.0
    if len(x1) == 1:
        return linear_interpolation(x2, m[0, :], vv2), OUT_OF_TABLE_DIFF if outside else 0.0
    if len(x2) == 1:
        return linear_interpolation(x1, m[:, 0], vv1), OUT_OF_TABLE_DIFF if outside else 0.0
    val, diff = bilinear_interpolation(x1, x2, m, vv1, vv2)
    return val, OUT_OF_TABLE_DIFF if outside else diff
