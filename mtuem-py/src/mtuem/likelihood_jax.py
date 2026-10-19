"""JAX versions of the MTUEM likelihoods and of the latent class mixture.

These functions return per-row LOG-likelihoods written with ``jax.numpy``, so the estimation
engine can take exact gradients, individual scores and Hessians instead of numerical ones. Use
them inside the probabilities function and create the inputs with
``make_inputs(..., autodiff=True)``; ``estimate`` then switches to exact derivatives on its own.

Covered (``PORTED``): ``theta_phi``, ``alpha_beta`` and ``additive_quadratic`` (profiled,
``sig``, ``sig`` + ``rho`` or ``cholesky`` errors; the additive quadratic form also with
``omitted_times`` / ``omitted_goods`` and ``optimal_tw=False``), plus ``class_alloc`` and ``lc``.
Not covered: the profiled latent class likelihood (``lc_profiled``), whose inner fixed point is
not differentiable as written. ``require_jax_form`` tells a model which forms must fall back to
NumPy and warns (``FallbackWarning``) when they do.
"""

import warnings
from collections.abc import Mapping

import numpy as np

try:
    import jax
    import jax.numpy as jnp
    from jax.scipy.linalg import solve_triangular
    from jax.scipy.special import logsumexp

    jax.config.update("jax_enable_x64", True)
except ImportError as err:  # pragma: no cover
    raise ImportError("mtuem.likelihood_jax needs JAX, a required dependency of mtuem: reinstall mtuem") from err

from . import _engine, likelihood
from ._engine import ClassProbs, FallbackWarning
from .utils import as_names

ESTIMATE_LIKE = ("estimate", "conditionals", "raw", "components", "output", "report")

PORTED: frozenset[str] = frozenset({"theta_phi", "alpha_beta", "additive_quadratic"})

_warned: set[str] = set()


def require_jax_form(form: str) -> bool:
    """True when the likelihood ``form`` has a JAX version.

    Otherwise warns once per form with ``FallbackWarning`` and returns False: the caller then
    uses the NumPy likelihood with numerical derivatives.
    """
    if form in PORTED:
        return True
    if form not in _warned:
        _warned.add(form)
        warnings.warn(
            f"The {form} likelihood has no JAX version; falling back to NumPy with numerical derivatives.",
            FallbackWarning,
            stacklevel=2,
        )
    return False


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


def _other_functionality(numpy_fn, settings: Mapping, functionality: str, inputs, delegated: tuple[str, ...]):
    """Result of a non-estimation functionality, or None for the estimation ones.

    ``delegated`` functionalities are answered by the numpy likelihood ``numpy_fn``.
    """
    if functionality in delegated:
        return numpy_fn(_to_numpy_settings(settings), functionality, inputs)
    N = len(inputs.database)
    if functionality == "validate":
        return jnp.zeros(N)
    if functionality == "zero_LL":
        return jnp.full(N, jnp.nan)
    if functionality not in ESTIMATE_LIKE:
        raise ValueError("Invalid value of argument 'functionality'")
    return None


def _context(s: Mapping):
    """tau, Tc, Ec and w as JAX arrays."""
    return (
        s["tau"],
        jnp.asarray(s["Tc"], dtype=float),
        jnp.asarray(s["Ec"], dtype=float),
        jnp.asarray(s["w"], dtype=float),
    )


def _stack(cols: list) -> "jnp.ndarray":
    return jnp.column_stack([c if c.ndim == 2 else c[:, None] for c in cols])


def _check_cholesky(s: Mapping, M: int) -> None:
    """Validate the cholesky length early, before any functionality, as the numpy code does."""
    n = _vec(s.get("cholesky")).size
    if n and n != M * (M + 1) // 2:
        raise ValueError(f"The length of cholesky ({n}) must match M * (M + 1) / 2 where M is the number of equations ({M}).")


def _safe_sqrt(x):
    """(sqrt(x) where x >= 0 and 0 elsewhere, x >= 0), with finite derivatives everywhere.

    The argument is replaced by 1 before the square root wherever it is not positive, so the
    derivative of the discarded branch is never inf or NaN (a single ``where`` would still
    propagate NaN gradients through the unused branch). Callers mark the rows where x < 0,
    which numpy computes as NaN, with ``_mark_rows``.
    """
    positive = x > 0
    root = jnp.where(positive, jnp.sqrt(jnp.where(positive, x, 1.0)), 0.0)
    return root, x >= 0


def _mark_rows(E, ok):
    """Residuals with the rows where ``ok`` is False set to NaN (numpy's value there)."""
    return jnp.where(jnp.broadcast_to(ok, (E.shape[0],))[:, None], E, jnp.nan)


def _mvn_logpdf_chol(E, L) -> "jnp.ndarray":
    """Zero-mean normal log density with covariance L L' (L lower triangular)."""
    M = L.shape[0]
    z = solve_triangular(L, E.T, lower=True)
    logdet = 2 * jnp.sum(jnp.log(jnp.abs(jnp.diag(L))))
    return -0.5 * (jnp.sum(z**2, axis=0) + logdet + M * jnp.log(2 * jnp.pi))


def _loglik_from_errors(E, s: Mapping) -> "jnp.ndarray":
    """Per-row normal log-likelihood of the residuals ``E`` under the error structure of ``s``.

    With no ``sig``/``rho``/``cholesky`` the covariance is profiled from the residuals, as in
    the R package. The numpy/R code projects the profiled correlation matrix with ``nearPD``
    and then gives zero likelihood when an eigenvalue is <= 1e-6. A sample correlation matrix
    is positive semi-definite, so the projection only acts when an eigenvalue is
    <= 1e-6 times the largest one, and then leaves an eigenvalue of about 1e-8 times the
    largest, which the guard rejects. Here the same points are rejected directly (eigenvalue
    <= 1e-6, or <= 1e-6 times the largest when profiled), so both backends agree without the
    projection.

    Rows with a non-finite residual (for example a negative discriminant in the optimal work
    time) get -inf. The profiled covariance uses all rows, so one such row makes the whole
    profiled likelihood -inf; numpy instead profiles from the complete rows only (as R's
    ``use = "complete.obs"``) and returns NaN for the incomplete ones.

    Degenerate points get -inf with finite (zero) derivatives: the error structure is
    replaced by the identity, and bad rows by zeros, before any Cholesky factor or log is
    taken, so no NaN reaches the gradient of a latent class mixture in which only one class is
    degenerate.
    """
    M = E.shape[1]
    chol = _vec(s.get("cholesky"))
    sig = _vec(s.get("sig"))
    rho = _vec(s.get("rho"))
    row_bad = ~jnp.all(jnp.isfinite(jax.lax.stop_gradient(E)), axis=1)
    E = jnp.where(row_bad[:, None], 0.0, E)
    profiled = not chol.size and not sig.size
    if chol.size:
        rows, cols = np.tril_indices(M)
        L = jnp.zeros((M, M)).at[rows, cols].set(chol)
        Sigma = L @ L.T
    elif profiled:
        # Standard deviations and correlations both come from the residuals.
        centred = E - E.mean(axis=0)
        Sigma = centred.T @ centred / (E.shape[0] - 1)
    else:
        # Explicit sds; correlations from rho (upper triangle, row-wise) or zero.
        R = jnp.eye(M)
        if rho.size:
            iu = np.triu_indices(M, 1)
            R = R.at[iu].set(rho).at[(iu[1], iu[0])].set(rho)
        Sigma = R * jnp.outer(sig, sig)

    # The same degeneracy guards as the numpy/R likelihood: zero likelihood when a standard
    # deviation is <= 1e-6 (an explicit sig must also be positive) or when the implied
    # correlation matrix has an eigenvalue <= 1e-6. The guards only select, so they are
    # computed without derivatives.
    S = jax.lax.stop_gradient(Sigma)
    var = jnp.diag(S)
    bad = ~jnp.all(jnp.isfinite(S)) | jnp.any(var <= 1e-12)
    if profiled:
        bad = bad | jnp.any(row_bad)
    if sig.size and not chol.size:
        bad = bad | jnp.any(jax.lax.stop_gradient(sig) <= 1e-6)
    if M > 1:
        sd = jnp.sqrt(jnp.where(var > 0, var, 1.0))
        eig = jnp.linalg.eigvalsh(jnp.where(bad, jnp.eye(M), S / jnp.outer(sd, sd)))
        bad = bad | (eig[0] <= 1e-6)
        if profiled:
            bad = bad | (eig[0] <= 1e-6 * eig[-1])

    if chol.size:
        L = jnp.where(bad, jnp.eye(M), L)
    else:
        L = jnp.linalg.cholesky(jnp.where(bad, jnp.eye(M), Sigma))
    ll = _mvn_logpdf_chol(E, L)
    return jnp.where(bad | row_bad | ~jnp.isfinite(ll), -jnp.inf, ll)


# ----------------------------------------------------------------------------------------------
# Theta-Phi
# ----------------------------------------------------------------------------------------------
def _theta_phi_residuals(s: Mapping, inputs) -> "jnp.ndarray":
    """Observed minus optimal allocations (Theta-Phi), differentiable."""
    we = s["work_elasticities"]
    Theta, Phi, thw = we["Theta"], we["Phi"], we["thw"]
    tau, Tc, Ec, w = _context(s)
    ec = Ec / (w * (tau - Tc))
    base = Theta + Phi + thw
    tpe = Phi + thw + (Theta + thw) * ec
    root, ok = _safe_sqrt(tpe**2 - 4 * thw * ec * base)
    tw_opt = 0.5 * (tau - Tc) * (tpe + root) / base

    work = as_names(s["work_times"])
    free_times = as_names(s.get("free_times"))
    free_goods = as_names(s.get("free_goods"))
    tw_used = tw_opt if s.get("optimal_tw", True) else jnp.asarray(inputs.matrix(work)[:, 0])
    cols = [tw_opt]
    th = _vec(s.get("times_elasticities"))
    if th.size:
        cols.append(jnp.outer(tau - tw_used - Tc, th / Theta))
    ph = _vec(s.get("goods_elasticities"))
    if ph.size:
        cost = _vec(s.get("goods_cost", 1))
        cols.append(jnp.outer(w * tw_used - Ec, (ph / cost) / Phi))
    obs = jnp.asarray(inputs.matrix(work + free_times + free_goods))
    return _mark_rows(obs - _stack(cols), ok)


def _n_equations(s: Mapping) -> int:
    return 1 + len(as_names(s.get("free_times"))) + len(as_names(s.get("free_goods")))


def theta_phi(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Per-row log-likelihood of the Theta-Phi MTUEM, the JAX counterpart of ``likelihood.theta_phi``.

    Takes the same settings. ``prediction`` and ``preprocess`` are answered by the numpy code.
    """
    s = dict(settings)
    _check_cholesky(s, _n_equations(s))
    out = _other_functionality(likelihood.theta_phi, s, functionality, inputs, ("prediction", "preprocess"))
    if out is not None:
        return out
    return _loglik_from_errors(_theta_phi_residuals(s, inputs), s)


# ----------------------------------------------------------------------------------------------
# Alpha-beta
# ----------------------------------------------------------------------------------------------
def _alpha_beta_residuals(s: Mapping, inputs) -> "jnp.ndarray":
    """Observed minus optimal allocations (alpha-beta), differentiable."""
    we = s["work_elasticities"]
    alpha, beta = we["alpha"], we["beta"]
    tau, Tc, Ec, w = _context(s)
    ec = Ec / (w * (tau - Tc))
    betaalphaec = beta + alpha * ec
    root, ok = _safe_sqrt(betaalphaec**2 - ec * (2 * alpha + 2 * beta - 1))
    tw_opt = (tau - Tc) * (betaalphaec + root)

    work = as_names(s["work_times"])
    free_times = as_names(s.get("free_times"))
    free_goods = as_names(s.get("free_goods"))
    cols = [tw_opt]
    th = _vec(s.get("times_elasticities"))
    if th.size:
        cols.append(jnp.outer(tau - tw_opt - Tc, th / (1 - 2 * beta)))
    ph = _vec(s.get("goods_elasticities"))
    if ph.size:
        cost = _vec(s.get("goods_cost", 1))
        cols.append(jnp.outer(w * tw_opt - Ec, (ph / cost) / (1 - alpha)))
    obs = jnp.asarray(inputs.matrix(work + free_times + free_goods))
    return _mark_rows(obs - _stack(cols), ok)


def alpha_beta(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Per-row log-likelihood of the alpha-beta MTUEM, the JAX counterpart of ``likelihood.alpha_beta``.

    Takes the same settings. ``prediction``, ``preprocess`` and ``get_covar`` are answered by
    the numpy code.
    """
    s = dict(settings)
    _check_cholesky(s, _n_equations(s))
    out = _other_functionality(likelihood.alpha_beta, s, functionality, inputs, ("prediction", "preprocess", "get_covar"))
    if out is not None:
        return out
    return _loglik_from_errors(_alpha_beta_residuals(s, inputs), s)


# ----------------------------------------------------------------------------------------------
# Additive quadratic
# ----------------------------------------------------------------------------------------------
def _additive_quadratic_residuals(s: Mapping, all_cols: list[str], eq_names: list[str], inputs) -> "jnp.ndarray":
    """Observed minus optimal allocations of the equations kept (additive quadratic), differentiable."""
    theta = _vec(s["times_elasticities"])
    beta = _vec(s["times_satiations"])
    phi = _vec(s["goods_elasticities"])
    eta = _vec(s["goods_satiations"])
    P = _vec(s["goods_cost"])
    if P.size == 1:
        P = jnp.repeat(P, eta.size)
    Sth = jnp.sum(theta / beta)
    Sph = jnp.sum(P * phi / eta)
    B = jnp.sum(1 / beta)
    H = jnp.sum(P**2 / eta)

    tau, Tc, Ec, w = _context(s)
    num_1 = _vec(s["work_elasticity"])[0]
    num_2 = w * (Sph + Ec) / H
    num_3 = (Sth - (tau - Tc)) / B
    den = _vec(s["work_satiation"])[0] + w**2 / H + 1 / B
    tw_opt = (num_1 + num_2 - num_3) / den

    work = as_names(s["work_times"])
    tw_used = tw_opt if s["optimal_tw"] else jnp.asarray(inputs.matrix(work[:1])[:, 0])
    mu = (Sth - (tau - tw_used - Tc)) / B
    ti = theta / beta - jnp.outer(mu, 1 / beta)
    lam = (Sph - (w * tw_used - Ec)) / H
    xj = phi / eta - jnp.outer(lam, P / eta)
    opt = _stack([tw_opt, ti, xj])[:, [all_cols.index(c) for c in eq_names]]
    return jnp.asarray(inputs.matrix(eq_names)) - opt


def additive_quadratic(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Per-row log-likelihood of the additive quadratic MTUEM, the JAX counterpart of
    ``likelihood.additive_quadratic``.

    Takes the same settings, including ``omitted_times`` / ``omitted_goods`` and
    ``optimal_tw``. ``prediction`` and ``preprocess`` are answered by the numpy code.
    """
    s = likelihood._with_defaults(settings, likelihood._AQ_DEFAULTS)
    all_cols, eq_names = likelihood._additive_quadratic_equations(s)
    _check_cholesky(s, len(eq_names))
    out = _other_functionality(likelihood.additive_quadratic, dict(settings), functionality, inputs, ("prediction", "preprocess"))
    if out is not None:
        return out
    return _loglik_from_errors(_additive_quadratic_residuals(s, all_cols, eq_names, inputs), s)


# ----------------------------------------------------------------------------------------------
# Latent classes
# ----------------------------------------------------------------------------------------------
def class_alloc(classAlloc_settings: Mapping) -> ClassProbs:
    """Class allocation probabilities (multinomial logit), the JAX counterpart of the engine's ``class_alloc``."""
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


def lc(lc_settings: Mapping, inputs, functionality: str = "estimate"):
    """Latent class mixture in log space: log sum_s pi_s exp(sum_rows log L_s).

    ``inClassProb`` holds per-row class LOG-likelihoods (from the functions above),
    ``classProb`` the class probabilities (from ``class_alloc``). Returns one log-likelihood per
    individual.
    """
    in_class = lc_settings["inClassProb"]
    comps = [v for k, v in in_class.items() if k != "model"] if isinstance(in_class, Mapping) else list(in_class)
    class_prob = lc_settings["classProb"]
    pis = list(class_prob)
    log_pis = getattr(class_prob, "log", None) or [_safe_log(p) for p in pis]
    if functionality in ("prediction", "preprocess"):
        converted = {"inClassProb": [np.asarray(c) for c in comps], "classProb": [np.asarray(p) for p in pis]}
        return _engine.lc(converted, inputs, functionality)
    n_ind = inputs.nIndivs
    idx = jnp.asarray(inputs.indiv_index)
    first = jnp.asarray(inputs.first_rows)
    terms = []
    for ll, log_pi in zip(comps, log_pis):
        ll = jnp.asarray(ll)
        if inputs.panel:
            ll = jax.ops.segment_sum(ll, idx, num_segments=n_ind)
            log_pi = jnp.broadcast_to(jnp.asarray(log_pi, dtype=float), (inputs.nObs,))[first]
        terms.append(log_pi + ll)
    return logsumexp(jnp.stack(terms, axis=1), axis=1)
