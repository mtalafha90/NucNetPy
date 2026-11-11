"""NucNetPy's network equations at fixed temperature and density, as arrays.

Within one burning stage the temperature and density are fixed, so the rate
of almost every reaction is a constant.  ``nucnetpy.solver.rhs`` nevertheless
evaluates all 10,787 rates and flows one by one in Python at every call,
which costs about 0.1 s; the implicit Euler stepper needs tens of thousands
of calls per stage.

:class:`FrozenNetwork` asks NucNetPy for each reaction's rate once, through
``Reaction.bare_rate`` (so the fitted rates and the partition-function
factors of the reverse rates are NucNetPy's), and then evaluates NucNetPy's
flow formula,

    F_r = lambda_r * rho**(n_r - 1) / prod_i(m_ir!) * prod_i Y_i**m_ir,

for all reactions at once with sparse matrices.  Abundances are clipped at
zero exactly as ``nucnetpy.solver.rhs`` clips them.  Two things depend on
the composition and are re-evaluated at every call, as NucNetPy does: the
rates of the three electron captures, which NucNetPy multiplies by rho * Ye
(``Reaction.bare_rate`` is asked for them again with the current Ye), and,
when screening is used, the factor of NucNetPy's ``SkyNetScreening``.
``check()`` compares the result with ``nucnetpy.solver.rhs`` and
``nucnetpy.solver.analytic_jacobian``.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

SLOTS = 4          # no ReacLib rate has more than four reactant nuclei


class FrozenNetwork:
    def __init__(self, network, species: Sequence[str], t9: float, rho: float,
                 screening=None):
        """``screening`` is the screening model NucNetPy would use, or None."""
        self.network = network
        self.species = list(species)
        self.t9, self.rho = float(t9), float(rho)
        self.screening = screening
        n = len(self.species)
        idx = {s: i for i, s in enumerate(self.species)}
        self.z = np.array([network.species[s].z for s in self.species], dtype=float)

        reactions = network.reactions.reactions
        r_count = len(reactions)
        slot_idx = np.full((r_count, SLOTS), n, dtype=np.int64)   # n = the "1" slot
        slot_m = np.zeros((r_count, SLOTS), dtype=np.int64)
        base = np.zeros(r_count)
        captures = []
        rows, cols, vals = [], [], []
        charge = np.zeros((r_count, SLOTS), dtype=np.int64)
        for r, reaction in enumerate(reactions):
            for s, p in enumerate(reaction.nuclear_reactants):
                slot_idx[r, s] = idx[p.species]
                slot_m[r, s] = p.count
                charge[r, s] = network.species[p.species].z
            for name, nu in reaction.stoichiometry().items():
                if name in idx:
                    rows.append(idx[name])
                    cols.append(r)
                    vals.append(nu)
            # The bare rate, with the partition-function factor of a library
            # reverse rate.  Screening is applied per call below, and so is
            # the rate of an electron capture, which depends on Ye.
            lam = reaction.bare_rate(self.t9, rho=self.rho)
            order = reaction.reactant_order
            scale = self.rho ** max(order - 1, 0) / max(reaction.statistical_factor(), 1)
            base[r] = lam * scale
            if reaction.electron_count:
                captures.append((r, reaction, scale))
        self.slot_idx, self.slot_m, self.base = slot_idx, slot_m, base
        self.captures = captures
        self.stoich = csr_matrix((vals, (rows, cols)), shape=(n, r_count))
        self.charge = charge
        self.screened = screening is not None
        self.n = n

    # ------------------------------------------------------------------
    def coefficients(self, y: np.ndarray) -> np.ndarray:
        """Rate coefficients at composition ``y`` (already clipped at zero)."""
        k = self.base.copy()
        ye = float(self.z @ y)
        for r, reaction, scale in self.captures:
            k[r] = reaction.bare_rate(self.t9, rho=self.rho, ye=ye) * scale
        if self.screened:
            abund = {s: float(v) for s, v in zip(self.species, y)}
            self.screening.update(abund, self.t9, self.rho, ye)
            k *= self.screening_factors()
        return k

    def screening_factors(self) -> np.ndarray:
        """NucNetPy's SkyNet factor for every reaction, vectorised."""
        sc = self.screening
        mu = np.array([sc.chemical_potential(z) for z in range(sc.max_z + 1)])
        zs = np.clip(self.charge, 0, sc.max_z)
        total_mu = (mu[zs] * np.where(self.charge > 0, self.slot_m, 0)).sum(axis=1)
        z_tot = np.clip(np.where(self.charge > 0, self.charge * self.slot_m, 0).sum(axis=1),
                        0, sc.max_z)
        exponent = np.clip(total_mu - mu[z_tot], -200.0, 200.0)
        return np.where(z_tot > 0, np.exp(exponent), 1.0)

    def _powers(self, y: np.ndarray) -> np.ndarray:
        yext = np.append(y, 1.0)
        return yext[self.slot_idx] ** self.slot_m

    def rhs(self, y: np.ndarray) -> np.ndarray:
        yc = np.clip(np.asarray(y, dtype=float), 0.0, None)
        flux = self.coefficients(yc) * np.prod(self._powers(yc), axis=1)
        return self.stoich @ flux

    def jacobian(self, y: np.ndarray):
        yc = np.clip(np.asarray(y, dtype=float), 0.0, None)
        k = self.coefficients(yc)
        yext = np.append(yc, 1.0)
        pw = yext[self.slot_idx] ** self.slot_m
        rows, cols, vals = [], [], []
        r_all = np.arange(len(k))
        for s in range(SLOTS):
            m = self.slot_m[:, s]
            use = m > 0
            if not np.any(use):
                continue
            others = np.prod(np.delete(pw, s, axis=1), axis=1)
            deriv = k * m * yext[self.slot_idx[:, s]] ** np.maximum(m - 1, 0) * others
            rows.append(r_all[use])
            cols.append(self.slot_idx[use, s])
            vals.append(deriv[use])
        dfdy = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(len(k), self.n)).tocsr()
        return (self.stoich @ dfdy).tocsc()

    # ------------------------------------------------------------------
    def check(self, y: np.ndarray) -> dict:
        """Largest relative differences from NucNetPy's own functions at ``y``."""
        from nucnetpy.solver import analytic_jacobian, rhs
        thermo = lambda t, ab: (self.t9, self.rho)  # noqa: E731
        ref_f = rhs(self.network, self.species, thermo, screening=self.screening)(0.0, y)
        ref_j = analytic_jacobian(self.network, self.species, thermo,
                                  screening=self.screening, sparse=True)(0.0, y).toarray()
        f = self.rhs(y)
        j = self.jacobian(y).toarray()
        scale_f = np.max(np.abs(ref_f))
        scale_j = np.abs(ref_j).max(axis=1, keepdims=True) + 1e-300
        return {"rhs": float(np.max(np.abs(f - ref_f)) / scale_f),
                "jacobian": float(np.max(np.abs(j - ref_j) / scale_j))}
