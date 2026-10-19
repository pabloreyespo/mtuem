"""Maximum likelihood estimation and the estimation results object."""

import math
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import pandas as pd

from ..config import Control
from .inputs import Inputs
from .latent import ProbFn, _first_row, _loglik_per_indiv
from .loglike import _LogLike, _model_component
from .optimize import _maximise, maximise_bounded


@dataclass
class Estimate:
    """Estimation results."""

    estimate: pd.Series
    beta_start: pd.Series
    fixed: list[str]
    control: Control
    maximum: float
    LLStart: float
    LL0: dict = field(default_factory=dict)
    LLC: float = float("nan")
    LLout: dict = field(default_factory=dict)
    se: pd.Series | None = None
    robse: pd.Series | None = None
    varcov: pd.DataFrame | None = None
    robvarcov: pd.DataFrame | None = None
    corrmat: pd.DataFrame | None = None
    robcorrmat: pd.DataFrame | None = None
    hessian: pd.DataFrame | None = None
    hessianEigenValue: np.ndarray | None = None
    eigValue: float = float("nan")
    code: int = 0
    message: str = ""
    successfulEstimation: bool = True
    routine: str = "bfgs"
    nIter: str | int = 0
    nObs: int = 0
    nIndivs: int = 0
    nFreeParams: int = 0
    AIC: float = float("nan")
    BIC: float = float("nan")
    startTime: datetime = field(default_factory=datetime.now)
    timePre: float = 0.0
    timeEst: float = 0.0
    timePost: float = 0.0
    timeTaken: float = 0.0
    classShares: pd.Series | None = None
    indLevelLL: np.ndarray | None = None
    scores: np.ndarray | None = None

    def __getitem__(self, key: str):
        return getattr(self, key)


def resolve_hessian(hessian: str | None, inputs: Inputs) -> str:
    """``"analytic"`` (exact, needs autodiff inputs), ``"numeric"`` or ``"none"``; None picks the best."""
    if hessian is None:
        return "analytic" if inputs.autodiff else "numeric"
    if hessian not in ("analytic", "numeric", "none"):
        raise ValueError("hessian must be 'analytic', 'numeric' or 'none'")
    if hessian == "analytic" and not inputs.autodiff:
        raise ValueError("hessian='analytic' requires autodiff inputs")
    return hessian


def _bound_arrays(
    bounds: Mapping[str, tuple[float | None, float | None]] | None, beta: pd.Series, free: list[str]
) -> tuple[np.ndarray, np.ndarray] | None:
    """Lower and upper bounds of the free parameters (+-inf where unbounded), or None if no bound applies.

    Bounds on fixed parameters are ignored. A free parameter whose starting value is not strictly
    inside its bounds raises ValueError.
    """
    if not bounds:
        return None
    unknown = [n for n in bounds if n not in beta.index]
    if unknown:
        raise ValueError(f"Bounds given for parameters not in beta: {', '.join(unknown)}")
    lower = np.full(len(free), -np.inf)
    upper = np.full(len(free), np.inf)
    for i, n in enumerate(free):
        lo, hi = bounds.get(n, (None, None))
        lower[i] = -np.inf if lo is None else float(lo)
        upper[i] = np.inf if hi is None else float(hi)
        if not lower[i] < upper[i]:
            raise ValueError(f"The lower bound of {n} is not below its upper bound")
        if not lower[i] < beta[n] < upper[i]:
            raise ValueError(f"The starting value of {n} ({beta[n]}) is not strictly inside its bounds")
    if not (np.isfinite(lower).any() or np.isfinite(upper).any()):
        return None
    return lower, upper


def _covariances(loglike: _LogLike, x: np.ndarray, hessian: str, H: np.ndarray | None = None):
    """Classical and robust covariance matrices, following ``apollo_varcov``."""
    free = loglike.free
    k = len(free)
    nan = pd.DataFrame(np.full((k, k), np.nan), index=free, columns=free)
    if hessian == "none" or k == 0:
        return dict(hessian=None, eig=None, varcov=nan, robvarcov=nan, scores=None)
    if H is None:
        H = loglike.hessian(x, hessian)
    if not np.all(np.isfinite(H)):
        return dict(hessian=None, eig=None, varcov=nan, robvarcov=nan, scores=None)
    eig = np.linalg.eigvals(H)
    eig = np.real_if_close(eig)
    try:
        varcov = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        varcov = np.full((k, k), np.nan)
    # The meat of the sandwich is the BHHH matrix, the uncentred cross-product of the
    # individual scores.
    score = loglike.scores(x)
    bhhh = score.T @ score
    robvarcov = varcov @ bhhh @ varcov
    return dict(
        hessian=pd.DataFrame(H, index=free, columns=free),
        eig=eig,
        varcov=pd.DataFrame(varcov, index=free, columns=free),
        robvarcov=pd.DataFrame(robvarcov, index=free, columns=free),
        scores=score,
    )


def _corr(vc: pd.DataFrame) -> pd.DataFrame:
    sd = np.sqrt(np.diag(vc.to_numpy()))
    with np.errstate(invalid="ignore", divide="ignore"):
        return pd.DataFrame(vc.to_numpy() / np.outer(sd, sd), index=vc.index, columns=vc.columns)


def _component_lls(beta: dict, prob_fn: ProbFn, inputs: Inputs) -> dict:
    P = prob_fn(beta, inputs, "output")
    out = {"model": float(np.sum(_loglik_per_indiv(_model_component(P), inputs)))}
    for k, v in P.items():
        if k != "model":
            out[k] = float(np.sum(_loglik_per_indiv(v, inputs)))
    return out


def _finalise(
    loglike: _LogLike,
    opt_x: np.ndarray,
    inputs: Inputs,
    beta: pd.Series,
    fixed: list[str],
    ll_start: float,
    code: int,
    message: str,
    routine: str,
    n_iter,
    t0: float,
    t_pre: float,
    t_est: float,
    hessian: str,
    start_time: datetime,
    H: np.ndarray | None = None,
) -> Estimate:
    free = loglike.free
    estimate = loglike.base.copy()
    estimate[free] = opt_x
    cov = _covariances(loglike, opt_x, hessian, H)
    varcov, robvarcov = cov["varcov"], cov["robvarcov"]
    se = pd.Series(np.nan, index=estimate.index)
    robse = pd.Series(np.nan, index=estimate.index)
    with np.errstate(invalid="ignore"):
        se[free] = np.sqrt(np.diag(varcov.to_numpy()))
        robse[free] = np.sqrt(np.diag(robvarcov.to_numpy()))
    maximum = loglike.total(opt_x)
    k = len(free)
    lls = _component_lls(estimate.to_dict(), loglike.prob_fn, inputs)
    shares = None
    if inputs.lc_pars is not None:
        pis = inputs.attach(estimate.to_dict()).pi_values
        shares = pd.Series(
            [float(np.mean(_first_row(p, inputs))) for p in pis], index=[f"class_{s + 1}" for s in range(len(pis))]
        )
    eig = cov["eig"]
    est = Estimate(
        estimate=estimate,
        beta_start=beta,
        fixed=list(fixed),
        control=inputs.control,
        maximum=maximum,
        LLStart=ll_start,
        LL0={name: float("nan") for name in lls},
        LLout=lls,
        se=se,
        robse=robse,
        varcov=varcov,
        robvarcov=robvarcov,
        corrmat=_corr(varcov),
        robcorrmat=_corr(robvarcov),
        hessian=cov["hessian"],
        hessianEigenValue=eig,
        eigValue=float(np.max(np.real(eig))) if eig is not None else float("nan"),
        code=code,
        message=message,
        successfulEstimation=code == 0,
        routine=routine,
        nIter=n_iter,
        nObs=inputs.nObs,
        nIndivs=inputs.nIndivs,
        nFreeParams=k,
        AIC=-2 * maximum + 2 * k,
        BIC=-2 * maximum + k * math.log(inputs.nObs),
        startTime=start_time,
        timePre=t_pre,
        timeEst=t_est,
        classShares=shares,
        indLevelLL=loglike.per_indiv(opt_x),
        scores=cov["scores"],
    )
    est.timeTaken = time.perf_counter() - t0
    est.timePost = est.timeTaken - t_pre - t_est
    return est


def estimate(
    beta: Mapping[str, float],
    fixed: list[str],
    prob_fn: ProbFn,
    inputs: Inputs,
    *,
    max_iterations: int = 200,
    max_function_evals: int = 100000,
    hessian: str | None = None,
    saddle_escapes: int = 3,
    silent: bool = False,
    bounds: Mapping[str, tuple[float | None, float | None]] | None = None,
) -> Estimate:
    """Maximum likelihood estimation with a BFGS line-search optimiser.

    ``hessian`` is ``"analytic"`` (exact, autodiff inputs only), ``"numeric"`` (numDeriv), or
    ``"none"``; by default the exact Hessian is used with autodiff inputs and numDeriv otherwise.
    A converged run reports code 0 and the message "successful convergence".

    ``bounds`` maps parameter names to ``(lower, upper)``, None for no bound. When a bound applies
    to a free parameter the optimiser is the adaptive barrier of ``maximise_bounded`` (R's
    ``constrOptim``, as maxLik uses for Apollo constraints). The covariance matrix is computed
    from the Hessian of the log-likelihood itself, without the barrier, also at a binding bound.
    """
    t0 = time.perf_counter()
    start_time = datetime.now()
    beta = pd.Series(beta, dtype=float)
    fixed = list(fixed or [])
    free = [n for n in beta.index if n not in fixed]
    limits = _bound_arrays(bounds, beta, free)
    hessian = resolve_hessian(hessian, inputs)
    # Parameters are scaled by |start|; flooring at 1 keeps near-zero starts (e.g. after EM)
    # from shrinking the optimiser steps to nothing.
    scale = np.maximum(np.abs(beta[free].to_numpy()), 1.0)
    loglike = _LogLike(beta, fixed, prob_fn, inputs, scale)
    x0 = beta[free].to_numpy()
    ll_start = loglike.total(x0)
    if not silent:
        print(f"Initial LL: {ll_start:.4f}")
    if not np.isfinite(ll_start):
        raise FloatingPointError("Log-likelihood calculation fails at starting values!")
    t_pre = time.perf_counter() - t0

    escapes = int(saddle_escapes) if hessian != "none" else 0
    x_start = x0 / scale
    iterations = 0
    H = None
    feasible = None
    if limits is not None:
        lower, upper = limits

        def feasible(x):
            return bool(np.all(x > lower) and np.all(x < upper))

    for attempt in range(escapes + 1):
        opt_kwargs = dict(
            max_iterations=int(max_iterations),
            max_function_evals=int(max_function_evals),
            silent=True,
            exact_grad_and_scores=loglike.grad_and_scores_scaled if loglike.autodiff else None,
        )
        if limits is None:
            res = _maximise(lambda z: loglike.per_indiv(z, scaled=True), x_start, **opt_kwargs)
        else:
            res = maximise_bounded(
                lambda z: loglike.per_indiv(z, scaled=True), x_start, lower / scale, upper / scale, **opt_kwargs
            )
        iterations += res.iterations
        if attempt == escapes or res.code != 0:
            break
        H = loglike.hessian(res.x * scale, hessian)
        x_escape = _escape_saddle(loglike, res.x * scale, res.f, H, feasible)
        if x_escape is None:
            break
        if not silent:
            print("Converged to a saddle point; stepping along the direction of positive curvature and re-estimating.")
        H = None
        x_start = x_escape / scale
    t_est = time.perf_counter() - t0 - t_pre
    if not silent:
        how = "exact derivatives (JAX)" if loglike.autodiff else "numerical methods (numDeriv)"
        print(f"Estimation finished: {res.message}. Computing covariance matrix using {how}.")
    return _finalise(
        loglike,
        res.x * scale,
        inputs,
        beta,
        fixed,
        ll_start,
        res.code,
        res.message,
        "bfgs",
        iterations,
        t0,
        t_pre,
        t_est,
        hessian,
        start_time,
        H,
    )


def _escape_saddle(loglike: _LogLike, x: np.ndarray, f: float, H: np.ndarray, feasible=None) -> np.ndarray | None:
    """A point with a higher likelihood along the most positive-curvature direction, if any.

    BFGS stops wherever the gradient vanishes, including saddle points such as the symmetric
    solution of a latent class model started from identical classes. Numerical noise usually
    drifts off such points; this does it deliberately. With ``feasible`` only points it accepts
    (strictly inside the bounds) are candidates.
    """
    if not np.all(np.isfinite(H)):
        return None
    eigval, eigvec = np.linalg.eigh((H + H.T) / 2)
    if eigval[-1] <= 0:
        return None
    v = eigvec[:, -1]
    size = max(1.0, float(np.max(np.abs(x))))
    for t in (0.1, 0.03, 0.01, 0.003, 0.001, 0.0003):
        for sign in (1.0, -1.0):
            xn = x + sign * t * size * v
            if feasible is not None and not feasible(xn):
                continue
            if loglike.total(xn) > f + 1e-6:
                return xn
    return None
