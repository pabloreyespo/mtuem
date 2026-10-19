"""Expectation-Maximisation for latent class models."""

import time
from collections.abc import Mapping
from datetime import datetime

import numpy as np
import pandas as pd

from .estimate import Estimate, _bound_arrays, _finalise, estimate, resolve_hessian
from .inputs import Inputs
from .latent import ProbFn, _class_parts, _log_posterior
from .loglike import _LogLike
from .optimize import MESSAGES, _maximise, maximise_bounded


def _autodiff_em_objective(base: pd.Series, free: list[str], prob_fn: ProbFn, inputs: Inputs):
    """Jitted EM objective q(x, h) per individual and its exact Jacobian (autodiff backend)."""
    import jax
    import jax.numpy as jnp

    from ..likelihood_jax import _safe_log

    base = {n: float(v) for n, v in base.items()}
    idx = jnp.asarray(inputs.indiv_index)
    first = jnp.asarray(inputs.first_rows)

    def class_terms(x):
        b = dict(base)
        for i, n in enumerate(free):
            b[n] = x[i]
        pis = inputs.attach(b).pi_values
        log_pis = getattr(pis, "log", None) or [_safe_log(pi) for pi in pis]
        comps = [jnp.ravel(jnp.asarray(v)) for k, v in prob_fn(b, inputs, "estimate").items() if k != "model"]
        if inputs.panel:
            comps = [jax.ops.segment_sum(c, idx, num_segments=inputs.nIndivs) for c in comps]
            log_pis = [jnp.broadcast_to(jnp.asarray(lp), (inputs.nObs,))[first] for lp in log_pis]
        return jnp.stack(comps, axis=1) + jnp.stack([jnp.broadcast_to(lp, comps[0].shape) for lp in log_pis], axis=1)

    def q(x, h):
        return jnp.sum(jnp.where(h > 0, h * class_terms(x), 0.0), axis=1)

    return jax.jit(q), jax.jit(jax.jacfwd(q, argnums=0))


def lc_em(
    beta: Mapping[str, float],
    fixed: list[str],
    prob_fn: ProbFn,
    inputs: Inputs,
    *,
    max_iterations: int = 100,
    tolerance: float = 1e-5,
    post_em: int = 2,
    estimate_kwargs: Mapping | None = None,
    bounds: Mapping[str, tuple[float | None, float | None]] | None = None,
) -> Estimate:
    """Expectation-Maximisation for latent class models (follows the EM of Apollo).

    Each M step maximises the expected complete-data log-likelihood
    ``sum_n sum_s h_ns [log P_s(n) + log pi_s(n)]`` over all free parameters. Parameters shared
    across classes are allowed. ``post_em`` controls what follows the EM iterations: 0 returns
    the EM estimates without covariance, 1 adds the covariance matrix at the EM estimates and 2
    (default) runs a full ``estimate`` from them, with ``estimate_kwargs``.

    ``bounds`` (as in ``estimate``) are imposed in every M step with ``maximise_bounded`` and are
    passed on to the final ``estimate``.
    """
    t0 = time.perf_counter()
    start_time = datetime.now()
    est_kwargs = dict(estimate_kwargs or {})
    beta = pd.Series(beta, dtype=float)
    fixed = list(fixed or [])
    free = [n for n in beta.index if n not in fixed]
    limits = _bound_arrays(bounds, beta, free)
    silent = bool(est_kwargs.get("silent", False))

    def m_step(fun_vec, z0, scale, **kwargs):
        if limits is None:
            return _maximise(fun_vec, z0, **kwargs)
        return maximise_bounded(fun_vec, z0, limits[0] / scale, limits[1] / scale, **kwargs)

    def full_ll(b: dict) -> float:
        return float(np.sum(_log_posterior(*_class_parts(b, prob_fn, inputs, "estimate"))[1]))

    current = beta.copy()
    ll = full_ll(current.to_dict())
    ll_start = ll
    if inputs.autodiff:
        q_ad, q_jac = _autodiff_em_objective(beta, free, prob_fn, inputs)
    if not silent:
        print("Initialising EM algorithm")
    em_iter = 0
    for em_iter in range(1, int(max_iterations) + 1):
        if not silent:
            print(f"Starting iteration: {em_iter}\n- Current LL : {ll:.2f}")
        h, _ = _log_posterior(*_class_parts(current.to_dict(), prob_fn, inputs, "estimate"))

        def q(z, h=h):
            b = current.copy()
            b[free] = z
            logPz, log_piz = _class_parts(b.to_dict(), prob_fn, inputs, "estimate")
            with np.errstate(invalid="ignore"):
                return np.sum(np.where(h > 0, h * (logPz + log_piz), 0.0), axis=1)

        scale = np.maximum(np.abs(current[free].to_numpy()), 1.0)
        if inputs.autodiff:
            def q_vec(z, sc=scale, h=h):
                return np.asarray(q_ad(z * sc, h))

            def q_grad(z, sc=scale, h=h):
                S = np.asarray(q_jac(z * sc, h)) * sc
                return (S.sum(axis=0), S) if np.all(np.isfinite(S)) else (None, None)

            res = m_step(q_vec, current[free].to_numpy() / scale, scale, max_iterations=200,
                         exact_grad_and_scores=q_grad)
        else:
            res = m_step(lambda z, sc=scale: q(z * sc), current[free].to_numpy() / scale, scale, max_iterations=200)
        current[free] = res.x * scale
        new_ll = full_ll(current.to_dict())
        change = new_ll - ll
        ll = new_ll
        if not silent:
            print(f"- Updated LL: {ll:.2f}")
        if abs(change) < tolerance:
            break
    t_em = time.perf_counter() - t0
    if not silent:
        print(f"EM algorithm stopped after {em_iter} iterations, LL = {ll:.2f}")

    if post_em == 2:
        est = estimate(current.to_dict(), fixed, prob_fn, inputs, bounds=bounds, **est_kwargs)
        routine = est.routine
        est.routine = f"EM algorithm ({routine}) -> Maximum likelihood ({routine})"
        est.nIter = f"{em_iter} (EM) & {est.nIter} ({routine})"
        est.beta_start = beta
        est.LLStart = ll_start
        est.startTime = start_time
        est.timePre += t_em
        est.timeTaken += t_em
        return est

    loglike = _LogLike(current, fixed, prob_fn, inputs)
    return _finalise(
        loglike,
        current[free].to_numpy(),
        inputs,
        beta,
        fixed,
        ll_start,
        0,
        MESSAGES[0],
        "EM algorithm",
        f"{em_iter} (EM)",
        t0,
        0.0,
        t_em,
        resolve_hessian(None, inputs) if post_em == 1 else "none",
        start_time,
    )
