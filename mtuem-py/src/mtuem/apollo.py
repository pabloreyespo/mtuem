"""A small, Apollo-like estimation engine for continuous-density models such as the MTUEM.

The R version of mtuem only supplies likelihood functions and relies on the Apollo package for
everything else. This module provides the pieces of Apollo those likelihoods need, keeping
Apollo's workflow and the layout of its printed and saved results:

* ``validate_inputs``  ~ ``apollo_validateInputs``
* ``ApolloInputs.attach`` ~ ``apollo_attach`` (returns a namespace instead of attaching)
* ``prepare_prob``     ~ ``apollo_prepareProb``
* ``estimate``         ~ ``apollo_estimate``
* ``lc_em``            ~ ``apollo_lcEM``
* ``model_output``     ~ ``apollo_modelOutput``
* ``save_output``      ~ ``apollo_saveOutput``
* ``delta_method``     ~ ``apollo_deltaMethod`` (expressions are Python callables)
* ``prediction``       ~ ``apollo_prediction``
* ``class_alloc``, ``lc``, ``lc_conditionals`` ~ their Apollo counterparts

Standard errors use a port of ``numDeriv::hessian`` and Apollo's sandwich estimator, so they
match the R output closely.
"""

import getpass
import io
import math
import platform
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from . import numderiv

try:
    __version__ = version("mtuem")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "dev"

ProbFn = Callable[..., dict]


# ==============================================================================================
# Control and inputs
# ==============================================================================================
@dataclass
class ApolloControl:
    """Model-level settings, the equivalent of ``apollo_control``."""

    modelName: str = "model"
    modelDescr: str = ""
    indivID: str | None = None
    outputDirectory: str = "."
    noValidation: bool = False
    noDiagnostics: bool = False
    nCores: int = 1
    # Apollo's option, default FALSE there. mtuem-py defaults to TRUE: likelihood functions
    # return log-likelihoods, so tiny densities never underflow to zero.
    workInLogs: bool = True

    @classmethod
    def from_any(cls, control) -> "ApolloControl":
        if isinstance(control, ApolloControl):
            return control
        known = {k: v for k, v in dict(control).items() if k in cls.__dataclass_fields__}
        return cls(**known)


class ApolloInputs:
    """The equivalent of ``apollo_inputs``: database, control, parameters and lcPars."""

    def __init__(
        self,
        database: pd.DataFrame,
        apollo_beta: Mapping[str, float],
        apollo_fixed: list[str],
        apollo_control: ApolloControl,
        apollo_lcPars: Callable | None = None,
        autodiff: bool = False,
    ):
        self.database = database.reset_index(drop=True)
        # True when apollo_probabilities is written with mtuem.autodiff (exact derivatives)
        self.autodiff = autodiff
        if autodiff and not apollo_control.workInLogs:
            raise ValueError("The autodiff backend always works in logs: do not set workInLogs = False")
        # True when likelihood functions return log-likelihoods (apollo_control$workInLogs)
        self.work_in_logs = bool(apollo_control.workInLogs) or autodiff
        self.apollo_beta = pd.Series(apollo_beta, dtype=float)
        self.apollo_fixed = list(apollo_fixed)
        self.apollo_control = apollo_control
        self.apollo_lcPars = apollo_lcPars
        self.obs_matrix = None
        self._columns = {}
        for col in self.database.columns:
            s = self.database[col]
            if pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s):
                self._columns[col] = s.to_numpy(dtype=float)
            else:
                self._columns[col] = s.to_numpy()
        ids = self.database[apollo_control.indivID].to_numpy() if apollo_control.indivID else np.arange(len(database))
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

    def attach(self, apollo_beta: Mapping[str, float]) -> SimpleNamespace:
        """Database columns, parameter values and (for latent class) lcPars in one namespace."""
        beta = {k: (float(v) if isinstance(v, (int, float, np.number)) else v) for k, v in dict(apollo_beta).items()}
        ns = SimpleNamespace(**self._columns, **beta)
        if self.apollo_lcPars is not None:
            for k, v in self.apollo_lcPars(ns, self).items():
                setattr(ns, k, v)
        return ns


def validate_inputs(
    database: pd.DataFrame,
    apollo_beta: Mapping[str, float],
    apollo_fixed: list[str] | None = None,
    apollo_control: ApolloControl | Mapping | None = None,
    apollo_lcPars: Callable | None = None,
    silent: bool = False,
    autodiff: bool = False,
) -> ApolloInputs:
    """Check and bundle the model inputs, as ``apollo_validateInputs``.

    With ``autodiff=True`` the model is written with ``mtuem.autodiff``: ``apollo_probabilities``
    returns log-likelihoods and estimation uses exact (JAX) gradients, scores and Hessians.
    """
    if autodiff:
        from . import autodiff as _ad  # noqa: F401  (fails early with a clear message without JAX)
    control = ApolloControl.from_any(apollo_control or {})
    fixed = list(apollo_fixed or [])
    missing = [f for f in fixed if f not in apollo_beta]
    if missing:
        raise ValueError(f"apollo_fixed contains names not in apollo_beta: {', '.join(missing)}")
    if control.indivID and control.indivID not in database.columns:
        raise ValueError(f"indivID '{control.indivID}' is not a column of the database")
    for name, value in dict(apollo_beta).items():
        if not np.isfinite(float(value)):
            raise ValueError(f"Starting value of '{name}' is not finite")
    inputs = ApolloInputs(database, apollo_beta, fixed, control, apollo_lcPars, autodiff)
    if not silent:
        print("All checks on apollo_control completed.")
        print("All checks on database completed.")
    return inputs


def prepare_prob(P: dict, apollo_inputs: ApolloInputs, functionality: str = "estimate") -> dict:
    """Validate the list of model components, as ``apollo_prepareProb``."""
    if "model" not in P:
        raise ValueError("The list of probabilities must contain a component named 'model'")
    return P


# ==============================================================================================
# Latent class helpers
# ==============================================================================================
class ClassProbs(list):
    """Class probabilities (a list, as ``apollo_classAlloc`` returns) that also carry exact logs.

    ``lc`` and the posterior computations use ``.log``, so a class whose probability underflows
    to 0 keeps a finite log-probability.
    """

    def __init__(self, probs, log):
        super().__init__(probs)
        self.log = log


def class_alloc(classAlloc_settings: Mapping) -> ClassProbs:
    """Class allocation probabilities from utilities (a multinomial logit), as ``apollo_classAlloc``."""
    utilities = classAlloc_settings["utilities"]
    V = list(utilities.values()) if isinstance(utilities, Mapping) else list(utilities)
    V = np.broadcast_arrays(*[np.atleast_1d(np.asarray(v, dtype=float)) for v in V])
    V = np.column_stack(V)
    V = V - V.max(axis=1, keepdims=True)
    log_probs = V - np.log(np.exp(V).sum(axis=1, keepdims=True))
    S = log_probs.shape[1]
    return ClassProbs([np.exp(log_probs[:, s]) for s in range(S)], [log_probs[:, s] for s in range(S)])


def _log_class_probs(class_prob) -> list[np.ndarray]:
    """Exact log class probabilities when available (``ClassProbs``), else log of the levels."""
    logs = getattr(class_prob, "log", None)
    if logs is not None:
        return [np.asarray(lp, dtype=float) for lp in logs]
    with np.errstate(divide="ignore"):
        return [np.log(np.asarray(p, dtype=float)) for p in class_prob]


def _panel_prod(x: np.ndarray, inputs: ApolloInputs) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == inputs.nIndivs and not inputs.panel:
        return x
    with np.errstate(divide="ignore"):
        return np.exp(np.bincount(inputs.indiv_index, weights=np.log(x), minlength=inputs.nIndivs))


def _first_row(x: np.ndarray, inputs: ApolloInputs) -> np.ndarray:
    x = np.broadcast_to(np.asarray(x, dtype=float), (inputs.nObs,))
    return x[inputs.first_rows] if inputs.panel else np.asarray(x)


def _class_components(in_class: Any) -> list:
    if isinstance(in_class, Mapping):
        return [v for k, v in in_class.items() if k != "model"]
    return list(in_class)


def lc(lc_settings: Mapping, apollo_inputs: ApolloInputs, functionality: str = "estimate"):
    """Latent class mixture of in-class likelihoods, as ``apollo_lc``.

    ``lc_settings`` holds ``inClassProb`` (dict or list of class likelihoods) and ``classProb``
    (list of class allocation probabilities). With panel data each class likelihood is
    multiplied over the rows of an individual before mixing.
    """
    in_class = _class_components(lc_settings["inClassProb"])
    class_prob = lc_settings["classProb"]
    pis = list(class_prob)
    if len(in_class) != len(pis):
        raise ValueError("inClassProb and classProb must have the same number of classes")
    n = apollo_inputs.nObs
    if functionality == "preprocess":
        return {"componentName": lc_settings.get("componentName", "model")}
    logs = getattr(apollo_inputs, "work_in_logs", False)
    if functionality == "validate":
        return np.zeros(n) if logs else np.ones(n)
    if functionality == "zero_LL":
        return np.full(n, np.nan)
    if functionality == "prediction":
        frames = []
        for s, (p, pi) in enumerate(zip(in_class, pis), start=1):
            df = pd.DataFrame(p).add_suffix(f"_class{s}")
            df[f"prob_class{s}"] = np.broadcast_to(np.asarray(pi, dtype=float), (n,))
            frames.append(df)
        return pd.concat(frames, axis=1)
    if functionality == "gradient":
        return {"like": None, "grad": None}
    if logs:
        # Mixture in logs: log sum_s exp(log pi_s + sum_rows log L_s), one value per individual.
        terms = [_first_row(lp, apollo_inputs) + _loglik_per_indiv(p, apollo_inputs)
                 for p, lp in zip(in_class, _log_class_probs(class_prob))]
        return _log_posterior(np.column_stack(terms), np.zeros((apollo_inputs.nIndivs, len(terms))))[1]
    total = np.zeros(apollo_inputs.nIndivs)
    for p, pi in zip(in_class, pis):
        total += _first_row(pi, apollo_inputs) * _panel_prod(p, apollo_inputs)
    return total


def _loglik_per_indiv(L, inputs: ApolloInputs) -> np.ndarray:
    L = np.asarray(L, dtype=float).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        logL = L if getattr(inputs, "work_in_logs", False) else np.log(L)
    if L.size == inputs.nObs and inputs.panel:
        return np.bincount(inputs.indiv_index, weights=logL, minlength=inputs.nIndivs)
    return logL


def _class_parts(beta: Mapping, prob_fn: ProbFn, inputs: ApolloInputs, functionality: str = "conditionals"):
    """Per-individual class LOG-likelihoods and log class allocation probabilities at ``beta``.

    Working in logs keeps posteriors and EM weights accurate when class likelihoods are tiny
    (for example with a participation penalty), for both the numpy and the autodiff backend.
    """
    if inputs.apollo_lcPars is None:
        raise ValueError("apollo_lcPars is required for latent class functions")
    pis = inputs.attach(beta).pi_values
    P = prob_fn(dict(beta), inputs, functionality)
    comps = [v for k, v in P.items() if k != "model"]
    if len(comps) != len(pis):
        raise ValueError(
            "Could not find one component per class in the list returned by apollo_probabilities. "
            "Latent class models estimated with mtuem_lc_profiled_likelihood have no class components."
        )
    logPs = np.column_stack([_loglik_per_indiv(c, inputs) for c in comps])
    log_pis = np.column_stack([_first_row(lp, inputs) for lp in _log_class_probs(pis)])
    return logPs, log_pis


def _log_posterior(logPs: np.ndarray, log_pis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Posterior class probabilities and per-individual mixture log-likelihood, stably."""
    a = logPs + log_pis
    m = np.max(a, axis=1, keepdims=True)
    with np.errstate(invalid="ignore"):
        e = np.exp(a - m)
        tot = e.sum(axis=1, keepdims=True)
    return e / tot, (m + np.log(tot)).ravel()


def lc_conditionals(model: "ApolloModel", apollo_probabilities: ProbFn, apollo_inputs: ApolloInputs) -> pd.DataFrame:
    """Posterior class membership probabilities, as ``apollo_lcConditionals``."""
    post, _ = _log_posterior(*_class_parts(model.estimate.to_dict(), apollo_probabilities, apollo_inputs))
    ids = apollo_inputs.ids[apollo_inputs.first_rows] if apollo_inputs.panel else apollo_inputs.ids
    out = pd.DataFrame(post, columns=[f"X{s + 1}" for s in range(post.shape[1])])
    out.insert(0, "ID", ids)
    return out


# ==============================================================================================
# Optimiser
# ==============================================================================================
BGW_MESSAGES = {
    3: "X-convergence (favorable)",
    4: "Relative function convergence (favorable)",
    7: "Singular convergence (unfavorable)",
    8: "False convergence (unfavorable)",
    9: "Function evaluation limit reached without convergence (unfavorable)",
    10: "Iteration limit reached without convergence (unfavorable)",
}


@dataclass
class _OptResult:
    x: np.ndarray
    f: float
    code: int
    message: str
    iterations: int
    function_evals: int


def _maximise(
    fun_vec: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    max_iterations: int = 500,
    max_function_evals: int = 100000,
    rel_f_tol: float = 1e-10,
    rel_grad_tol: float = 1e-6,
    silent: bool = True,
    exact_grad_and_scores: Callable | None = None,
) -> _OptResult:
    """Maximise a sum of per-individual log-likelihoods, ``sum(fun_vec(x))``.

    Quasi-Newton (BFGS) with a backtracking line search. As in BGW, the curvature model is
    seeded with the BHHH matrix (outer product of the individual scores), which the central
    differences used for the gradient deliver at no extra cost. Non-finite values are treated
    as an over-long step, so the search backs away from regions where the likelihood is zero.
    Convergence codes follow BGW: 4 is relative function convergence (the gain predicted by
    the quadratic model is below ``rel_f_tol`` times |f|).
    """
    evals = 0  # function evaluations outside the gradient, as BGW counts them

    def fvec(x):
        v = np.asarray(fun_vec(x), dtype=float)
        return v if np.all(np.isfinite(v)) else None

    def f(x):
        nonlocal evals
        evals += 1
        v = fvec(x)
        return float(v.sum()) if v is not None else -np.inf

    def grad_and_scores(x):
        n = x.size
        S = None
        for i in range(n):
            h = 1e-5 * max(1.0, abs(x[i]))
            e = np.zeros(n)
            e[i] = h
            up, down = fvec(x + e), fvec(x - e)
            if up is None or down is None:
                return None, None
            col = (up - down) / (2 * h)
            if S is None:
                S = np.zeros((col.size, n))
            S[:, i] = col
        return S.sum(axis=0), S

    if exact_grad_and_scores is not None:
        grad_and_scores = exact_grad_and_scores  # noqa: F811  (autodiff replaces finite differences)

    def bhhh_inverse(S):
        B = S.T @ S
        B += 1e-8 * np.trace(B) / B.shape[0] * np.eye(B.shape[0])
        try:
            return np.linalg.inv(B)
        except np.linalg.LinAlgError:
            return np.eye(B.shape[0])

    x = np.asarray(x0, dtype=float).copy()
    fx = f(x)
    if not np.isfinite(fx):
        raise FloatingPointError("Log-likelihood is not finite at the starting values")
    g, S = grad_and_scores(x)
    if g is None:
        raise FloatingPointError("Log-likelihood is not finite next to the starting values")
    n = x.size
    Hinv = bhhh_inverse(S)
    restarted = True
    code = 10
    it = 0
    for it in range(1, max_iterations + 1):
        direction = Hinv @ g
        if direction @ g <= 0:
            Hinv = bhhh_inverse(S)
            restarted = True
            direction = Hinv @ g
        slope = direction @ g
        step = 1.0
        accepted = False
        while step > 1e-14:
            xn = x + step * direction
            fn = f(xn)
            if np.isfinite(fn) and fn >= fx + 1e-4 * step * slope:
                accepted = True
                break
            step *= 0.5 if np.isfinite(fn) else 0.1
        gn = Sn = None
        if accepted:
            gn, Sn = grad_and_scores(xn)
            accepted = gn is not None
        if not accepted:
            if not restarted:
                Hinv = bhhh_inverse(S)
                restarted = True
                continue
            # No step improves the function: converged if the predicted gain is negligible
            # (numerical noise), otherwise false convergence, as BGW reports it.
            code = 4 if 0.5 * slope / max(abs(fx), 1.0) < 1e-7 else 8
            break
        s = xn - x
        y = -(gn - g)  # change in the gradient of -f, the function BFGS minimises
        sy = s @ y
        if sy > 1e-12 * np.linalg.norm(s) * np.linalg.norm(y):
            rho = 1.0 / sy
            V = np.eye(n) - rho * np.outer(s, y)
            Hinv = V @ Hinv @ V.T + rho * np.outer(s, s)
        restarted = False
        rel_change = abs(fn - fx) / max(abs(fn), 1.0)
        x, fx, g, S = xn, fn, gn, Sn
        pred = 0.5 * float(g @ Hinv @ g) / max(abs(fx), 1.0)
        rel_grad = float(np.max(np.abs(g) * np.maximum(np.abs(x), 1.0))) / max(abs(fx), 1.0)
        if not silent:
            print(f"  it {it:4d}   LL {fx: .6f}   predicted rel. gain {pred:.2e}   rel. gradient {rel_grad:.2e}")
        if pred < rel_f_tol and rel_change < 1e3 * rel_f_tol:
            if rel_grad < 1e2 * rel_grad_tol:
                code = 4
                break
            # The quasi-Newton matrix has lost curvature along directions where the gradient
            # is still large: re-seed it with BHHH.
            Hinv = bhhh_inverse(S)
            restarted = True
        if evals >= max_function_evals:
            code = 9
            break
    return _OptResult(x, fx, code, BGW_MESSAGES[code], it, evals)


# ==============================================================================================
# Model object
# ==============================================================================================
@dataclass
class ApolloModel:
    """Estimation results. Field names follow the list returned by ``apollo_estimate``."""

    estimate: pd.Series
    apollo_beta: pd.Series
    apollo_fixed: list[str]
    apollo_control: ApolloControl
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
    estimationRoutine: str = "bgw"
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


# ==============================================================================================
# Estimation
# ==============================================================================================
DEFAULT_ESTIMATE_SETTINGS = {
    "estimationRoutine": "bgw",
    "maxIterations": 200,
    "writeIter": False,
    "hessianRoutine": "numDeriv",
    "silent": False,
    "scaling": None,
    "bgw_settings": {},
}


class _LogLike:
    """Log-likelihood wrapper over the free parameters, optionally scaled."""

    def __init__(self, apollo_beta: pd.Series, apollo_fixed: list[str], prob_fn: ProbFn, inputs: ApolloInputs, scale=None):
        self.names = list(apollo_beta.index)
        self.free = [n for n in self.names if n not in apollo_fixed]
        self.base = apollo_beta.astype(float).copy()
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
            ll = jnp.ravel(jnp.asarray(prob_fn(b, inputs, "estimate")["model"]))
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

    def hessian(self, x_free: np.ndarray) -> np.ndarray:
        """Hessian of the total log-likelihood: exact with autodiff, else numDeriv with fallback."""
        if self.autodiff:
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
        return _loglik_per_indiv(P["model"], self.inputs)

    def total(self, x_free: np.ndarray, scaled: bool = False) -> float:
        ll = self.per_indiv(x_free, scaled)
        return float(np.sum(ll)) if np.all(np.isfinite(ll)) else -np.inf


def _hessian(f: Callable[[np.ndarray], float], x: np.ndarray) -> np.ndarray:
    """``numDeriv::hessian``, falling back to small steps when it fails, as Apollo does.

    numDeriv's first Richardson step moves each parameter by 10% of its value, which can leave
    the region where the likelihood is positive (e.g. a correlation of -0.96 stepped past -1).
    Apollo then falls back to maxLik's finite-difference Hessian; here the same Richardson
    scheme is rerun with steps of 1e-4 relative.
    """
    H = numderiv.hessian(f, x)
    if not np.all(np.isfinite(H)):
        H = numderiv.hessian(f, x, d=1e-4)
    return H


def _covariances(loglike: _LogLike, x: np.ndarray, hessian_routine: str, H: np.ndarray | None = None):
    """Classical and robust covariance matrices, following ``apollo_varcov``."""
    free = loglike.free
    k = len(free)
    nan = pd.DataFrame(np.full((k, k), np.nan), index=free, columns=free)
    if hessian_routine == "none" or k == 0:
        return dict(hessian=None, eig=None, varcov=nan, robvarcov=nan, scores=None)
    if H is None:
        H = loglike.hessian(x)
    if not np.all(np.isfinite(H)):
        return dict(hessian=None, eig=None, varcov=nan, robvarcov=nan, scores=None)
    eig = np.linalg.eigvals(H)
    eig = np.real_if_close(eig)
    try:
        varcov = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        varcov = np.full((k, k), np.nan)
    # Apollo's BGW route takes the meat from BGW's BHHH matrix, the uncentred cross-product
    # of the individual scores.
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


def _component_lls(beta: dict, prob_fn: ProbFn, inputs: ApolloInputs) -> dict:
    P = prob_fn(beta, inputs, "output")
    out = {"model": float(np.sum(_loglik_per_indiv(P["model"], inputs)))}
    for k, v in P.items():
        if k != "model":
            out[k] = float(np.sum(_loglik_per_indiv(v, inputs)))
    return out


def _finalise_model(
    loglike: _LogLike,
    opt_x: np.ndarray,
    inputs: ApolloInputs,
    apollo_beta: pd.Series,
    apollo_fixed: list[str],
    ll_start: float,
    code: int,
    message: str,
    routine: str,
    n_iter,
    t0: float,
    t_pre: float,
    t_est: float,
    hessian_routine: str,
    start_time: datetime,
    H: np.ndarray | None = None,
) -> ApolloModel:
    free = loglike.free
    estimate = loglike.base.copy()
    estimate[free] = opt_x
    cov = _covariances(loglike, opt_x, hessian_routine, H)
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
    if inputs.apollo_lcPars is not None:
        pis = inputs.attach(estimate.to_dict()).pi_values
        shares = pd.Series(
            [float(np.mean(_first_row(p, inputs))) for p in pis], index=[f"class_{s + 1}" for s in range(len(pis))]
        )
    eig = cov["eig"]
    model = ApolloModel(
        estimate=estimate,
        apollo_beta=apollo_beta,
        apollo_fixed=list(apollo_fixed),
        apollo_control=inputs.apollo_control,
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
        successfulEstimation=code in (0, 3, 4, 5, 6),
        estimationRoutine=routine,
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
    model.timeTaken = time.perf_counter() - t0
    model.timePost = model.timeTaken - t_pre - t_est
    return model


def estimate(
    apollo_beta: Mapping[str, float],
    apollo_fixed: list[str],
    apollo_probabilities: ProbFn,
    apollo_inputs: ApolloInputs,
    estimate_settings: Mapping | None = None,
) -> ApolloModel:
    """Maximum likelihood estimation, the equivalent of ``apollo_estimate``.

    ``estimationRoutine`` may be ``"bgw"`` (default) or ``"bfgs"``; both run the same BFGS
    line-search optimiser. A converged ``"bgw"`` run reports code 4 (relative function
    convergence) and a converged ``"bfgs"`` run reports code 0, the codes Apollo's two
    routines use, so code that checks ``model.code`` behaves as in R.
    """
    t0 = time.perf_counter()
    start_time = datetime.now()
    settings = dict(DEFAULT_ESTIMATE_SETTINGS)
    settings.update(estimate_settings or {})
    beta = pd.Series(apollo_beta, dtype=float)
    fixed = list(apollo_fixed or [])
    free = [n for n in beta.index if n not in fixed]
    scale = settings.get("scaling")
    if scale is None:
        # Apollo scales by |start|; flooring at 1 keeps near-zero starts (e.g. after EM) from
        # shrinking the optimiser steps to nothing.
        scale = np.maximum(np.abs(beta[free].to_numpy()), 1.0)
    else:
        scale = np.asarray([dict(scale).get(n, 1.0) for n in free], dtype=float)
    loglike = _LogLike(beta, fixed, apollo_probabilities, apollo_inputs, scale)
    silent = bool(settings.get("silent"))
    x0 = beta[free].to_numpy()
    ll_start = loglike.total(x0)
    if not silent:
        print(f"Initial LL: {ll_start:.4f}")
    if not np.isfinite(ll_start):
        raise FloatingPointError("Log-likelihood calculation fails at starting values!")
    t_pre = time.perf_counter() - t0

    routine = settings["estimationRoutine"].lower()
    if routine not in ("bgw", "bfgs"):
        raise ValueError("estimationRoutine must be 'bgw' or 'bfgs'")
    max_evals = dict(settings.get("bgw_settings") or {}).get("maxFunctionEvals", 100000)
    hessian_routine = settings.get("hessianRoutine", "numDeriv")
    escapes = int(settings.get("saddleEscapes", 3)) if hessian_routine != "none" else 0
    x_start = x0 / scale
    iterations = 0
    H = None
    for attempt in range(escapes + 1):
        res = _maximise(
            lambda z: loglike.per_indiv(z, scaled=True),
            x_start,
            max_iterations=int(settings["maxIterations"]),
            max_function_evals=int(max_evals),
            silent=silent or not settings.get("writeIter", False),
            exact_grad_and_scores=loglike.grad_and_scores_scaled if loglike.autodiff else None,
        )
        iterations += res.iterations
        if attempt == escapes or res.code != 4:
            break
        H = loglike.hessian(res.x * scale)
        x_escape = _escape_saddle(loglike, res.x * scale, res.f, H)
        if x_escape is None:
            break
        if not silent:
            print("Converged to a saddle point; stepping along the direction of positive curvature and re-estimating.")
        H = None
        x_start = x_escape / scale
    code = res.code
    if routine == "bfgs":
        code = 0 if res.code == 4 else 1
    t_est = time.perf_counter() - t0 - t_pre
    if not silent:
        how = "exact derivatives (JAX)" if loglike.autodiff else "numerical methods (numDeriv)"
        print(f"Estimation finished: {res.message}. Computing covariance matrix using {how}.")
    return _finalise_model(
        loglike,
        res.x * scale,
        apollo_inputs,
        beta,
        fixed,
        ll_start,
        code,
        res.message,
        routine,
        iterations,
        t0,
        t_pre,
        t_est,
        hessian_routine,
        start_time,
        H,
    )


def _escape_saddle(loglike: _LogLike, x: np.ndarray, f: float, H: np.ndarray) -> np.ndarray | None:
    """A point with a higher likelihood along the most positive-curvature direction, if any.

    BFGS stops wherever the gradient vanishes, including saddle points such as the symmetric
    solution of a latent class model started from identical classes. Apollo usually drifts off
    such points by numerical noise; this does it deliberately.
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
            if loglike.total(xn) > f + 1e-6:
                return xn
    return None


def _autodiff_em_objective(base: pd.Series, free: list[str], prob_fn: ProbFn, inputs: ApolloInputs):
    """Jitted EM objective q(x, h) per individual and its exact Jacobian (autodiff backend)."""
    import jax
    import jax.numpy as jnp

    from .autodiff import _safe_log

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
    apollo_beta: Mapping[str, float],
    apollo_fixed: list[str],
    apollo_probabilities: ProbFn,
    apollo_inputs: ApolloInputs,
    lcEM_settings: Mapping | None = None,
    estimate_settings: Mapping | None = None,
) -> ApolloModel:
    """Expectation-Maximisation for latent class models, the equivalent of ``apollo_lcEM``.

    Each M step maximises the expected complete-data log-likelihood
    ``sum_n sum_s h_ns [log P_s(n) + log pi_s(n)]`` over all free parameters. Unlike Apollo,
    parameters shared across classes are allowed. ``postEM`` controls what follows the EM
    iterations: 0 returns the EM estimates without covariance, 1 adds the covariance matrix at the
    EM estimates and 2 (default) runs a full ``estimate`` from them.
    """
    t0 = time.perf_counter()
    start_time = datetime.now()
    s = {"EMstoppingCriterion": 1e-5, "EMmaxIterations": 100, "postEM": 2, "silent": False}
    s.update(lcEM_settings or {})
    est_set = dict(estimate_settings or {})
    beta = pd.Series(apollo_beta, dtype=float)
    fixed = list(apollo_fixed or [])
    free = [n for n in beta.index if n not in fixed]
    inputs = apollo_inputs
    silent = bool(s["silent"])

    def full_ll(b: dict) -> float:
        return float(np.sum(_log_posterior(*_class_parts(b, apollo_probabilities, inputs, "estimate"))[1]))

    current = beta.copy()
    ll = full_ll(current.to_dict())
    ll_start = ll
    if inputs.autodiff:
        q_ad, q_jac = _autodiff_em_objective(beta, free, apollo_probabilities, inputs)
    if not silent:
        print("Initialising EM algorithm")
    em_iter = 0
    for em_iter in range(1, int(s["EMmaxIterations"]) + 1):
        if not silent:
            print(f"Starting iteration: {em_iter}\n- Current LL : {ll:.2f}")
        h, _ = _log_posterior(*_class_parts(current.to_dict(), apollo_probabilities, inputs, "estimate"))

        def q(z, h=h):
            b = current.copy()
            b[free] = z
            logPz, log_piz = _class_parts(b.to_dict(), apollo_probabilities, inputs, "estimate")
            with np.errstate(invalid="ignore"):
                return np.sum(np.where(h > 0, h * (logPz + log_piz), 0.0), axis=1)

        scale = np.maximum(np.abs(current[free].to_numpy()), 1.0)
        if inputs.autodiff:
            def q_vec(z, sc=scale, h=h):
                return np.asarray(q_ad(z * sc, h))

            def q_grad(z, sc=scale, h=h):
                S = np.asarray(q_jac(z * sc, h)) * sc
                return (S.sum(axis=0), S) if np.all(np.isfinite(S)) else (None, None)

            res = _maximise(q_vec, current[free].to_numpy() / scale, max_iterations=200, exact_grad_and_scores=q_grad)
        else:
            res = _maximise(lambda z, sc=scale: q(z * sc), current[free].to_numpy() / scale, max_iterations=200)
        current[free] = res.x * scale
        new_ll = full_ll(current.to_dict())
        change = new_ll - ll
        ll = new_ll
        if not silent:
            print(f"- Updated LL: {ll:.2f}")
        if abs(change) < s["EMstoppingCriterion"]:
            break
    t_em = time.perf_counter() - t0
    if not silent:
        print(f"EM algorithm stopped after {em_iter} iterations, LL = {ll:.2f}")

    if s["postEM"] == 2:
        est_set.setdefault("estimationRoutine", "bgw")
        model = estimate(current.to_dict(), fixed, apollo_probabilities, inputs, est_set)
        routine = model.estimationRoutine
        model.estimationRoutine = f"EM algorithm ({routine}) -> Maximum likelihood ({routine})"
        model.nIter = f"{em_iter} (EM) & {model.nIter} ({routine})"
        model.apollo_beta = beta
        model.LLStart = ll_start
        model.startTime = start_time
        model.timePre += t_em
        model.timeTaken += t_em
        return model

    loglike = _LogLike(current, fixed, apollo_probabilities, inputs)
    model = _finalise_model(
        loglike,
        current[free].to_numpy(),
        inputs,
        beta,
        fixed,
        ll_start,
        4,
        BGW_MESSAGES[4],
        "EM algorithm",
        f"{em_iter} (EM)",
        t0,
        0.0,
        t_em,
        "numDeriv" if s["postEM"] == 1 else "none",
        start_time,
    )
    return model


# ==============================================================================================
# Output
# ==============================================================================================
def _r_num(x, digits: int = 7) -> str:
    """Format a number the way R's ``cat`` prints it (up to 7 significant digits)."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, (float, np.floating)) and math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    s = f"{x:.{digits}g}"
    if "e" in s:
        mant, exp = s.split("e")
        s = f"{mant}e{exp[0]}{exp[1:].zfill(2)}"
    return s


def _format_column(values: np.ndarray, n_signif: int = 4, width: int = 11) -> list[str]:
    """Port of the column formatting in ``apollo_print``."""
    out, scie, ndec = [], [], []
    for v in values:
        if v is None or not np.isfinite(v):
            m = "NA".rjust(width) if (v is None or np.isnan(v)) else ("Inf" if v > 0 else "-Inf").rjust(width)
        else:
            m = f"%#{width}.{n_signif}g" % v
        is_sci = "e" in m
        nd = len(m) - m.index(".") - 1 if "." in m else 0
        if not is_sci and nd > 6 and np.isfinite(v):
            m = f"%{width}.{n_signif}e" % v
            is_sci = True
        out.append(m)
        scie.append(is_sci)
        ndec.append(nd)
    non_sci = [d for d, s in zip(ndec, scie) if not s]
    obj = min(max(non_sci), 6) if non_sci else 6
    for i, v in enumerate(values):
        if not scie[i] and ndec[i] < obj and v is not None and np.isfinite(v):
            out[i] = f"%{width}.{obj}f" % v
    w = max(len(m) for m in out)
    return [m.rjust(w) for m in out]


def apollo_print_matrix(M: pd.DataFrame, n_signif: int = 4, width: int = 11) -> str:
    """Render a numeric table exactly as ``apollo_print`` does (via R's data.frame print)."""
    cols = [_format_column(M[c].to_numpy(dtype=float), n_signif, width) for c in M.columns]
    rn = [str(r) for r in M.index]
    rw = max(len(r) for r in rn) if rn else 0
    widths = [max(len(str(c)), len(col[0]) if col else 0) for c, col in zip(M.columns, cols)]
    lines = [" " * rw + "".join(" " + str(c).rjust(w) for c, w in zip(M.columns, widths))]
    for i, r in enumerate(rn):
        lines.append(r.ljust(rw) + "".join(" " + col[i].rjust(w) for col, w in zip(cols, widths)))
    return "\n".join(lines)


def _hms(t: float) -> str:
    h = int(t // 3600)
    m = int((t - h * 3600) // 60)
    s = round(t - h * 3600 - m * 60, 2)
    return f"{h:02d}:{m:02d}:{_r_num(s)}"


def _estimates_table(model: ApolloModel) -> pd.DataFrame:
    with np.errstate(invalid="ignore", divide="ignore"):
        return pd.DataFrame(
            {
                "Estimate": model.estimate,
                "s.e.": model.se,
                "t.rat.(0)": model.estimate / model.se,
                "Rob.s.e.": model.robse,
                "Rob.t.rat.(0)": model.estimate / model.robse,
            }
        )


def _signif(x: float, digits: int) -> float:
    if x == 0 or not np.isfinite(x):
        return x
    return round(x, digits - int(math.floor(math.log10(abs(x)))) - 1)


def model_output_text(model: ApolloModel, printCovar: bool = False, printCorr: bool = False) -> str:
    """The text that ``model_output`` prints."""
    buf = io.StringIO()
    w = buf.write
    ctrl = model.apollo_control
    w(f"Model run by {getpass.getuser()} using mtuem-py {__version__} on Python "
      f"{platform.python_version()} for {platform.system()}.\n")
    w("Output layout follows Apollo; please acknowledge Apollo by citing Hess & Palma (2019)\n")
    w("  DOI 10.1016/j.jocm.2019.100170\n")
    w("  www.ApolloChoiceModelling.com\n\n")
    w(f"Model name                                  : {ctrl.modelName}\n")
    w(f"Model description                           : {ctrl.modelDescr}\n")
    w(f"Model run at                                : {model.startTime}\n")
    w(f"Estimation method                           : {model.estimationRoutine}\n")
    w(f"Estimation diagnosis                        : {model.message}\n")
    eig = model.hessianEigenValue
    if eig is None:
        w("Optimisation diagnosis                      : Unknown (no Hessian requested by user)\n")
    else:
        if np.iscomplexobj(eig):
            txt = ["Non-symmetrical hessian", "Complex eigenvalues", "NA"]
        else:
            if np.all(eig < 0):
                txt = ["Maximum found", "Negative definite"]
            elif np.any(eig == 0):
                txt = ["Inconclusive test", "Some eigenvalues are zero"]
            elif np.all(eig > 0):
                txt = ["Minimum found", "Positive definite"]
            else:
                txt = ["Saddle point found", "Some eigenvalues are positive and others negative"]
            txt.append(_r_num(round(float(np.max(eig)), 6), 15))
        if not np.iscomplexobj(eig) and np.all(eig < 0):
            txt.append(_r_num(_signif(float(np.max(eig) / np.min(eig)), 6), 15))
        else:
            txt.append("not calculated (Hessian is not negative definite)")
        w(f"Optimisation diagnosis                      : {txt[0]}\n")
        w(f"     hessian properties                     : {txt[1]}\n")
        w(f"     maximum eigenvalue                     : {txt[2]}\n")
        w(f"     reciprocal of condition number         : {txt[3]}\n")
    w(f"Number of individuals                       : {model.nIndivs}\n")
    w(f"Number of rows in database                  : {model.nObs}\n")
    w(f"Number of modelled outcomes                 : {model.nObs}\n\n")
    w(f"Number of cores used                        :  {ctrl.nCores} \n")
    w("Model without mixing\n\n")
    multi = len(model.LLout) > 1
    w(f"LL(start)                                   : {_r_num(round(model.LLStart, 2))}\n")
    if multi:
        w("LL (whole model) at equal shares, LL(0)     : NA\n")
        w("LL (whole model) at observed shares, LL(C)  : NA\n")
        w(f"LL(final, whole model)                      : {_r_num(round(model.maximum, 2))}\n")
    else:
        w("LL at equal shares, LL(0)                   : NA\n")
        w("LL at observed shares, LL(C)                : NA\n")
        w(f"LL(final)                                   : {_r_num(round(model.maximum, 2))}\n")
    for label in ("Rho-squared vs equal shares", "Adj.Rho-squared vs equal shares",
                  "Rho-squared vs observed shares", "Adj.Rho-squared vs observed shares"):
        w(f"{label.ljust(45)}:  Not applicable \n")
    w(f"AIC                                         :  {_r_num(round(model.AIC, 2))} \n")
    w(f"BIC                                         :  {_r_num(round(model.BIC, 2))} \n")
    if multi:
        w("\n")
        for nam, val in list(model.LLout.items())[1:]:
            sp1 = " " * (33 - len(nam) - 6) if 6 + len(nam) < 33 else ""
            sp2 = " " * (33 - len(nam) - 10) if 9 + len(nam) < 33 else ""
            w(f"LL(0,{nam}){sp1}: Not applicable\n")
            w(f"LL(final,{nam}){sp2}: {_r_num(round(val, 2))}\n")
    w("\n")
    w(f"Estimated parameters                        : {model.nFreeParams}\n")
    w(f"Time taken (hh:mm:ss)                       :  {_hms(model.timeTaken)} \n")
    w(f"     pre-estimation                         :  {_hms(model.timePre)} \n")
    w(f"     estimation                             :  {_hms(model.timeEst)} \n")
    w(f"     post-estimation                        :  {_hms(model.timePost)} \n")
    tmp = "" if model.successfulEstimation else f"({model.message})"
    w(f"Iterations                                  :  {model.nIter} {tmp} \n\n")
    w("Unconstrained optimisation.\n\n")
    w("Estimates:\n")
    w(apollo_print_matrix(_estimates_table(model)) + "\n\n")
    if model.classShares is not None:
        w("\nSummary of class allocation for model component :\n")
        shares = pd.DataFrame({"Mean prob.": model.classShares.round(4)})
        w(shares.to_string() + "\n\n")
    if printCovar:
        w("\nClassical covariance matrix:\n" + apollo_print_matrix(model.varcov) + "\n")
        w("\nRobust covariance matrix:\n" + apollo_print_matrix(model.robvarcov) + "\n")
    if printCorr:
        w("\nClassical correlation matrix:\n" + apollo_print_matrix(model.corrmat) + "\n")
        w("\nRobust correlation matrix:\n" + apollo_print_matrix(model.robcorrmat) + "\n")
    return buf.getvalue()


def model_output(model: ApolloModel, modelOutput_settings: Mapping | None = None) -> None:
    """Print estimation results in Apollo's layout, as ``apollo_modelOutput``."""
    s = dict(modelOutput_settings or {})
    print(model_output_text(model, s.get("printCovar", False), s.get("printCorr", False)), end="")


def save_output(model: ApolloModel, saveOutput_settings: Mapping | None = None) -> dict[str, Path]:
    """Write ``<modelName>_output.txt`` and ``<modelName>_estimates.csv``, as ``apollo_saveOutput``.

    With ``saveCov`` / ``saveCorr`` (default True) the classical and robust covariance and
    correlation matrices are also written to csv files, as Apollo does.
    """
    s = {"printCovar": True, "printCorr": True, "saveCov": True, "saveCorr": True}
    s.update(saveOutput_settings or {})
    outdir = Path(model.apollo_control.outputDirectory or ".")
    outdir.mkdir(parents=True, exist_ok=True)
    name = model.apollo_control.modelName
    files = {}
    files["output"] = outdir / f"{name}_output.txt"
    files["output"].write_text(model_output_text(model, s["printCovar"], s["printCorr"]), encoding="utf-8")
    est = _estimates_table(model)
    est.columns = ["Estimate", "Std.err.", "t-ratio(0)", "Rob.std.err.", "Rob.t-ratio(0)"]
    files["estimates"] = outdir / f"{name}_estimates.csv"
    est.to_csv(files["estimates"], na_rep="NA")
    if s["saveCov"]:
        for key, mat in (("covar", model.varcov), ("robcovar", model.robvarcov)):
            files[key] = outdir / f"{name}_{key}.csv"
            mat.to_csv(files[key], na_rep="NA")
    if s["saveCorr"]:
        for key, mat in (("corr", model.corrmat), ("robcorr", model.robcorrmat)):
            files[key] = outdir / f"{name}_{key}.csv"
            mat.to_csv(files[key], na_rep="NA")
    for f in files.values():
        print(f"Saved {f}")
    return files


# ==============================================================================================
# Post-estimation
# ==============================================================================================
def delta_method(model: ApolloModel, deltaMethod_settings: Mapping, silent: bool = False) -> pd.DataFrame:
    """Delta method standard errors for functions of the estimates, as ``apollo_deltaMethod``.

    ``deltaMethod_settings["expression"]`` maps a label to a callable that receives the full
    parameter vector as a dict (fixed parameters included) and returns a number. ``varcov`` is
    ``"robust"`` (default, as in Apollo) or ``"classical"``. Fixed parameters carry no variance.
    """
    expressions = deltaMethod_settings["expression"]
    which = deltaMethod_settings.get("varcov", "robust")
    vc = model.robvarcov if which == "robust" else model.varcov
    free = list(vc.index)
    b0 = model.estimate.copy()
    x0 = b0[free].to_numpy(dtype=float)
    rows = []
    for label, fn in expressions.items():
        def f(z, fn=fn):
            b = b0.copy()
            b[free] = z
            return float(fn(b.to_dict()))

        value = f(x0)
        g = numderiv.gradient(f, x0)
        se = float(np.sqrt(g @ vc.to_numpy() @ g))
        rows.append({"Expression": label, "Value": round(value, 4), "s.e.": round(se, 4),
                     "t-ratio (0)": round(value / se, 2) if se > 0 else float("nan")})
    out = pd.DataFrame(rows)
    if not silent:
        print(f"\nRunning Delta method computation for user-defined function using {which} standard errors\n")
        print(out.to_string(index=False))
    return out


def prediction(
    model: ApolloModel,
    apollo_probabilities: ProbFn,
    apollo_inputs: ApolloInputs,
    prediction_settings: Mapping | None = None,
) -> pd.DataFrame:
    """Model predictions at the estimates, as ``apollo_prediction``."""
    s = {"silent": False}
    s.update(prediction_settings or {})
    P = apollo_probabilities(model.estimate.to_dict(), apollo_inputs, "prediction")
    pred = pd.DataFrame(P["model"]).reset_index(drop=True)
    if not s["silent"]:
        print("Running predictions from model using estimated parameters...\n")
        print("Prediction at user provided parameters")
        summ = pd.DataFrame([pred.sum(), pred.mean()], index=["Aggregate", "Average"]).round(2)
        print(summ.to_string())
        print()
    ids = apollo_inputs.ids
    obs = apollo_inputs.database.groupby(apollo_inputs.indiv_index).cumcount().to_numpy() + 1
    pred.insert(0, "Observation", obs)
    pred.insert(0, "ID", ids)
    return pred
