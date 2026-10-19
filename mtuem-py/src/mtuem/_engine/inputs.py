"""Model inputs: the database, start values, control settings and latent class parameters."""

from collections.abc import Callable, Mapping
from types import SimpleNamespace

import numpy as np
import pandas as pd

from ..config import Control


class Inputs:
    """Database, control, parameters and latent class parameters (``lc_pars``) of a model."""

    def __init__(
        self,
        database: pd.DataFrame,
        beta: Mapping[str, float],
        fixed: list[str],
        control: Control,
        lc_pars: Callable | None = None,
        autodiff: bool = False,
    ):
        self.database = database.reset_index(drop=True)
        # True when the probabilities function is written with mtuem.likelihood_jax (exact derivatives)
        self.autodiff = autodiff
        # The engine always works in logs: likelihood functions return log-likelihoods, so tiny
        # densities never underflow to zero.
        self.work_in_logs = True
        self.beta = pd.Series(beta, dtype=float)
        self.fixed = list(fixed)
        self.control = control
        self.lc_pars = lc_pars
        self.obs_matrix = None
        self._columns = {}
        for col in self.database.columns:
            s = self.database[col]
            if pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s):
                self._columns[col] = s.to_numpy(dtype=float)
            else:
                self._columns[col] = s.to_numpy()
        ids = self.database[control.indiv_id].to_numpy() if control.indiv_id else np.arange(len(database))
        _, first, inverse = np.unique(ids, return_index=True, return_inverse=True)
        order = np.argsort(first)
        rank = np.empty_like(order)
        rank[order] = np.arange(order.size)
        self.indiv_index = rank[inverse]  # 0..nIndivs-1, in order of first appearance
        self.nIndivs = int(order.size)
        self.nObs = len(self.database)
        self.first_rows = np.sort(first)
        self.ids = ids

    def matrix(self, columns: list[str]) -> np.ndarray:
        """Database columns as a float matrix, cached because likelihoods ask for it every call."""
        key = tuple(columns)
        cache = self.__dict__.setdefault("_matrix_cache", {})
        if key not in cache:
            cache[key] = np.column_stack([np.asarray(self._columns[c], dtype=float) for c in columns])
        return cache[key]

    @property
    def panel(self) -> bool:
        return self.nIndivs < self.nObs

    def attach(self, beta: Mapping[str, float]) -> SimpleNamespace:
        """Database columns, parameter values and (for latent class) lc_pars in one namespace."""
        values = {k: (float(v) if isinstance(v, (int, float, np.number)) else v) for k, v in dict(beta).items()}
        ns = SimpleNamespace(**self._columns, **values)
        if self.lc_pars is not None:
            for k, v in self.lc_pars(ns, self).items():
                setattr(ns, k, v)
        return ns


def make_inputs(
    database: pd.DataFrame,
    beta: Mapping[str, float],
    fixed: list[str] | None,
    control: Control | None = None,
    lc_pars: Callable | None = None,
    autodiff: bool = False,
    silent: bool = True,
) -> Inputs:
    """Check and bundle the model inputs.

    With ``autodiff=True`` the model is written with ``mtuem.likelihood_jax``: the probabilities
    function returns log-likelihoods and estimation uses exact (JAX) gradients, scores and
    Hessians.
    """
    if autodiff:
        from .. import likelihood_jax as _lj  # noqa: F401  (fails early with a clear message without JAX)
    control = Control() if control is None else control
    fixed = list(fixed or [])
    missing = [f for f in fixed if f not in beta]
    if missing:
        raise ValueError(f"fixed contains names not in beta: {', '.join(missing)}")
    if control.indiv_id and control.indiv_id not in database.columns:
        raise ValueError(f"indiv_id '{control.indiv_id}' is not a column of the database")
    for name, value in dict(beta).items():
        if not np.isfinite(float(value)):
            raise ValueError(f"Starting value of '{name}' is not finite")
    inputs = Inputs(database, beta, fixed, control, lc_pars, autodiff)
    if not silent:
        print("All checks on control completed.")
        print("All checks on database completed.")
    return inputs
