"""Automatic differentiation backend (JAX) for the Theta-Phi MTUEM and latent class mixtures.

These functions return LOG-likelihoods written with ``jax.numpy``, so the estimation engine
can take exact gradients, individual scores and Hessians instead of numerical ones. Use them
inside ``apollo_probabilities`` and create the inputs with
``apollo.validate_inputs(..., autodiff=True)``; ``apollo.estimate`` then switches to exact
derivatives on its own.

Install the optional dependency with ``uv sync --extra jax`` (or ``pip install mtuem[jax]``).

Covered: ``mtuem_loglik`` (Theta-Phi with profiled, ``sig``, ``sig`` + ``rho`` or ``cholesky``
errors), ``class_alloc`` and ``lc``. Not covered: ``mtuem_lc_profiled_likelihood``, whose inner
fixed point is not differentiable as written.
"""

from collections.abc import Mapping

import numpy as np

try:
    import jax
    import jax.numpy as jnp
    from jax.scipy.linalg import solve_triangular
    from jax.scipy.special import logsumexp

    jax.config.update("jax_enable_x64", True)
except ImportError as err:  # pragma: no cover
    raise ImportError("mtuem.autodiff needs JAX: install it with `uv sync --extra jax`") from err

from . import apollo, likelihood
from .apollo import ClassProbs
from .utils import as_names

ESTIMATE_LIKE = ("estimate", "conditionals", "raw", "components", "output", "report")


def _vec(x) -> "jnp.ndarray":
    if x is None:
        return jnp.zeros(0)
    if isinstance(x, Mapping):
        x = list(x.values())
    if isinstance(x, (list, tuple)):
        return jnp.stack([jnp.asarray(v, dtype=float).reshape(()) for v in x]) if x else jnp.zeros(0)
    return jnp.atleast_1d(jnp.asarray(x, dtype=float))


def _to_numpy_settings(settings: Mapping) -> dict:
    """Concrete numpy copy of the settings, for the functionalities the numpy code handles."""

    def conv(v):
        if isinstance(v, Mapping):
            return {k: conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [conv(x) for x in v]
        if isinstance(v, (str, bool)) or v is None:
            return v
        a = np.asarray(v, dtype=float)
        return float(a) if a.ndim == 0 else a

    return {k: conv(v) for k, v in settings.items()}


def _residuals(s: Mapping, apollo_inputs) -> "jnp.ndarray":
    """Observed minus optimal allocations (Theta-Phi), differentiable."""
    we = s["work_elasticities"]
    Theta, Phi, thw = we["Theta"], we["Phi"], we["thw"]
    tau = s["tau"]
    Tc, Ec, w = jnp.asarray(s["Tc"], dtype=float), jnp.asarray(s["Ec"], dtype=float), jnp.asarray(s["w"], dtype=float)
    ec = Ec / (w * (tau - Tc))
    base = Theta + Phi + thw
    tpe = Phi + thw + (Theta + thw) * ec
    tw_opt = 0.5 * (tau - Tc) * (tpe + jnp.sqrt(tpe**2 - 4 * thw * ec * base)) / base

    work = as_names(s["work_times"])
    free_times = as_names(s.get("free_times"))
    free_goods = as_names(s.get("free_goods"))
    tw_used = tw_opt if s.get("optimal_tw", True) else jnp.asarray(apollo_inputs.matrix(work)[:, 0])
    cols = [tw_opt]
    th = _vec(s.get("times_elasticities"))
    if th.size:
        cols.append(jnp.outer(tau - tw_used - Tc, th / Theta))
    ph = _vec(s.get("goods_elasticities"))
    if ph.size:
        cost = _vec(s.get("goods_cost", 1))
        cols.append(jnp.outer(w * tw_used - Ec, (ph / cost) / Phi))
    opt = jnp.column_stack([c if c.ndim == 2 else c[:, None] for c in cols])
    obs = jnp.asarray(apollo_inputs.matrix(work + free_times + free_goods))
    return obs - opt


def _mvn_logpdf_chol(E, L) -> "jnp.ndarray":
    """Zero-mean normal log density with covariance L L' (L lower triangular)."""
    M = L.shape[0]
    z = solve_triangular(L, E.T, lower=True)
    logdet = 2 * jnp.sum(jnp.log(jnp.abs(jnp.diag(L))))
    return -0.5 * (jnp.sum(z**2, axis=0) + logdet + M * jnp.log(2 * jnp.pi))


def mtuem_loglik(mtuem_settings: Mapping, functionality: str = "estimate", apollo_inputs=None):
    """Per-row log-likelihood of the Theta-Phi MTUEM, the JAX counterpart of ``mtuem_likelihood``.

    Takes the same settings. With no ``sig``/``rho``/``cholesky`` the covariance is profiled
    from the residuals, as in the R package (the R code projects a non positive definite
    correlation matrix with ``nearPD``; here such a point gets log-likelihood -inf).
    """
    s = dict(mtuem_settings)
    if functionality in ("prediction", "preprocess"):
        return likelihood.mtuem_likelihood(_to_numpy_settings(s), functionality, apollo_inputs)
    N = len(apollo_inputs.database)
    if functionality == "validate":
        return jnp.zeros(N)
    if functionality == "zero_LL":
        return jnp.full(N, jnp.nan)
    if functionality not in ESTIMATE_LIKE:
        raise ValueError("Invalid value of argument 'functionality'")

    E = _residuals(s, apollo_inputs)
    M = E.shape[1]
    chol = _vec(s.get("cholesky"))
    sig = _vec(s.get("sig"))
    rho = _vec(s.get("rho"))
    if chol.size:
        if chol.size != M * (M + 1) // 2:
            raise ValueError(f"The length of cholesky ({chol.size}) must match M * (M + 1) / 2 where M is the number of equations ({M}).")
        rows, cols = np.tril_indices(M)
        L = jnp.zeros((M, M)).at[rows, cols].set(chol)
    elif sig.size == 0:
        # Profiled: standard deviations and correlations both come from the residuals.
        centred = E - E.mean(axis=0)
        L = jnp.linalg.cholesky(centred.T @ centred / (E.shape[0] - 1))
    else:
        # Explicit sds; correlations from rho (upper triangle, row-wise) or zero.
        R = jnp.eye(M)
        if rho.size:
            iu = np.triu_indices(M, 1)
            R = R.at[iu].set(rho).at[(iu[1], iu[0])].set(rho)
        L = jnp.linalg.cholesky(R * jnp.outer(sig, sig))
    ll = _mvn_logpdf_chol(E, L)
    # The same degeneracy guards as the numpy/R likelihood: zero likelihood when a standard
    # deviation is <= 1e-6 (an explicit sig must also be positive) or when the implied
    # correlation matrix has an eigenvalue <= 1e-6.
    Sigma = L @ L.T
    sd = jnp.sqrt(jnp.diag(Sigma))
    bad = ~jnp.all(jnp.isfinite(L)) | jnp.any(sd <= 1e-6)
    if sig.size and not chol.size:
        bad = bad | jnp.any(sig <= 1e-6)
    if M > 1:
        corr = Sigma / jnp.outer(sd, sd)
        bad = bad | (jnp.min(jnp.linalg.eigvalsh(corr)) <= 1e-6)
    return jnp.where(bad | ~jnp.isfinite(ll), -jnp.inf, ll)


def class_alloc(classAlloc_settings: Mapping) -> ClassProbs:
    """Class allocation probabilities (multinomial logit), the JAX counterpart of ``apollo.class_alloc``."""
    utilities = classAlloc_settings["utilities"]
    V = list(utilities.values()) if isinstance(utilities, Mapping) else list(utilities)
    V = jnp.broadcast_arrays(*[jnp.atleast_1d(jnp.asarray(v, dtype=float)) for v in V])
    V = jnp.column_stack(V)
    log_probs = V - logsumexp(V, axis=1, keepdims=True)
    S = log_probs.shape[1]
    return ClassProbs([jnp.exp(log_probs[:, s]) for s in range(S)], [log_probs[:, s] for s in range(S)])


def _safe_log(p):
    """log(p) with -inf (and zero derivative, not NaN) where p == 0."""
    p = jnp.asarray(p, dtype=float)
    positive = p > 0
    return jnp.where(positive, jnp.log(jnp.where(positive, p, 1.0)), -jnp.inf)


def lc(lc_settings: Mapping, apollo_inputs, functionality: str = "estimate"):
    """Latent class mixture in log space: log sum_s pi_s exp(sum_rows log L_s).

    ``inClassProb`` holds per-row class LOG-likelihoods (from ``mtuem_loglik``), ``classProb``
    the class probabilities (from ``class_alloc``). Returns one log-likelihood per individual.
    """
    in_class = lc_settings["inClassProb"]
    comps = [v for k, v in in_class.items() if k != "model"] if isinstance(in_class, Mapping) else list(in_class)
    class_prob = lc_settings["classProb"]
    pis = list(class_prob)
    log_pis = getattr(class_prob, "log", None) or [_safe_log(p) for p in pis]
    if functionality in ("prediction", "preprocess"):
        converted = {"inClassProb": [np.asarray(c) for c in comps], "classProb": [np.asarray(p) for p in pis]}
        return apollo.lc(converted, apollo_inputs, functionality)
    n_ind = apollo_inputs.nIndivs
    idx = jnp.asarray(apollo_inputs.indiv_index)
    first = jnp.asarray(apollo_inputs.first_rows)
    terms = []
    for ll, log_pi in zip(comps, log_pis):
        ll = jnp.asarray(ll)
        if apollo_inputs.panel:
            ll = jax.ops.segment_sum(ll, idx, num_segments=n_ind)
            log_pi = jnp.broadcast_to(jnp.asarray(log_pi, dtype=float), (apollo_inputs.nObs,))[first]
        terms.append(log_pi + ll)
    return logsumexp(jnp.stack(terms, axis=1), axis=1)
