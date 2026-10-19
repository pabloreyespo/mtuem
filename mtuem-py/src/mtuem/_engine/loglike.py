"""Log-likelihood wrapper over the free parameters, with numerical or exact (JAX) derivatives."""

from collections.abc import Callable

import numpy as np
import pandas as pd

from .. import numderiv
from .inputs import Inputs
from .latent import ProbFn, _loglik_per_indiv


class _LogLike:
    """Log-likelihood wrapper over the free parameters, optionally scaled."""

    def __init__(self, beta: pd.Series, fixed: list[str], prob_fn: ProbFn, inputs: Inputs, scale=None):
        self.names = list(beta.index)
        self.free = [n for n in self.names if n not in fixed]
        self.base = beta.astype(float).copy()
        self.prob_fn = prob_fn
        self.inputs = inputs
        self.scale = np.ones(len(self.free)) if scale is None else np.asarray(scale, dtype=float)
        self.autodiff = inputs.autodiff
        if self.autodiff:
            self._build_autodiff()

    def _build_autodiff(self) -> None:
        """Jitted per-individual log-likelihood, exact scores (jacfwd) and exact Hessian."""
        import jax
        import jax.numpy as jnp

        base = {n: float(v) for n, v in self.base.items()}
        free, inputs, prob_fn = self.free, self.inputs, self.prob_fn
        idx = jnp.asarray(inputs.indiv_index)

        def per_indiv(x):
            b = dict(base)
            for i, n in enumerate(free):
                b[n] = x[i]
            ll = jnp.ravel(jnp.asarray(_model_component(prob_fn(b, inputs, "estimate"))))
            if inputs.panel and ll.shape[0] == inputs.nObs:
                ll = jax.ops.segment_sum(ll, idx, num_segments=inputs.nIndivs)
            return ll

        self._ad_per_indiv = jax.jit(per_indiv)
        self._ad_scores = jax.jit(jax.jacfwd(per_indiv))
        self._ad_hessian = jax.jit(jax.hessian(lambda x: jnp.sum(per_indiv(x))))

    def scores(self, x_free: np.ndarray) -> np.ndarray:
        """Individual scores d log L_n / d beta (N x K): exact with autodiff, else numDeriv."""
        if self.autodiff:
            return np.asarray(self._ad_scores(np.asarray(x_free, dtype=float)))
        return numderiv.jacobian(lambda z: self.per_indiv(z), x_free)

    def hessian(self, x_free: np.ndarray, mode: str | None = None) -> np.ndarray:
        """Hessian of the total log-likelihood: ``"analytic"`` (exact, autodiff) or ``"numeric"``.

        The default is exact with autodiff and numDeriv otherwise.
        """
        if mode is None:
            mode = "analytic" if self.autodiff else "numeric"
        if mode == "analytic":
            if not self.autodiff:
                raise ValueError("hessian='analytic' requires autodiff inputs")
            return np.asarray(self._ad_hessian(np.asarray(x_free, dtype=float)))
        return _hessian(self.total, x_free)

    def grad_and_scores_scaled(self, z: np.ndarray):
        """Exact gradient and scores with respect to the scaled parameters (for the optimiser)."""
        S = self.scores(np.asarray(z) * self.scale) * self.scale
        if not np.all(np.isfinite(S)):
            return None, None
        return S.sum(axis=0), S

    def full(self, x_free: np.ndarray, scaled: bool = False) -> dict:
        b = self.base.copy()
        b[self.free] = np.asarray(x_free) * (self.scale if scaled else 1.0)
        return b.to_dict()

    def per_indiv(self, x_free: np.ndarray, scaled: bool = False) -> np.ndarray:
        if self.autodiff:
            x = np.asarray(x_free, dtype=float) * (self.scale if scaled else 1.0)
            return np.asarray(self._ad_per_indiv(x))
        P = self.prob_fn(self.full(x_free, scaled), self.inputs, "estimate")
        return _loglik_per_indiv(_model_component(P), self.inputs)

    def total(self, x_free: np.ndarray, scaled: bool = False) -> float:
        ll = self.per_indiv(x_free, scaled)
        return float(np.sum(ll)) if np.all(np.isfinite(ll)) else -np.inf


def _model_component(P: dict):
    """The 'model' entry of the dict a probabilities function returns."""
    if "model" not in P:
        raise ValueError("The dict of probabilities must contain a component named 'model'")
    return P["model"]


def _hessian(f: Callable[[np.ndarray], float], x: np.ndarray) -> np.ndarray:
    """``numDeriv::hessian``, falling back to small steps when it fails.

    numDeriv's first Richardson step moves each parameter by 10% of its value, which can leave
    the region where the likelihood is positive (e.g. a correlation of -0.96 stepped past -1).
    Apollo then falls back to maxLik's finite-difference Hessian; here the same Richardson
    scheme is rerun with steps of 1e-4 relative.
    """
    H = numderiv.hessian(f, x)
    if not np.all(np.isfinite(H)):
        H = numderiv.hessian(f, x, d=1e-4)
    return H
