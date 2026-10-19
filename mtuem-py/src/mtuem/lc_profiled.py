"""Latent class MTUEM with the class-specific error covariance profiled out of the likelihood."""

from collections.abc import Mapping

import numpy as np
import pandas as pd

from .likelihood import mtuem_likelihood, mtuem_residuals, works_in_logs
from .utils import cov_complete, cov_to_sig_rho

# Covariances, posterior weights and class shares of the last successful evaluation, the
# equivalent of the package-local environment used by the R version.
_PROFILE_STASH: dict = {}


def _mvn_logdensity(E: np.ndarray, Sigma: np.ndarray) -> np.ndarray | None:
    """Zero-mean multivariate normal log density through a Cholesky factor (None if degenerate)."""
    try:
        R = np.linalg.cholesky(Sigma)
    except np.linalg.LinAlgError:
        return None
    dg = np.diag(R)
    if np.any(~np.isfinite(dg)) or np.any(dg <= 1e-8):
        return None
    M = E.shape[1]
    z = np.linalg.solve(R, E.T)
    q = np.sum(z**2, axis=0)
    logdet = 2 * np.sum(np.log(dg))
    return -0.5 * (q + logdet + M * np.log(2 * np.pi))


def _mvn_density(E: np.ndarray, Sigma: np.ndarray) -> np.ndarray | None:
    """Zero-mean multivariate normal density (None if degenerate)."""
    ld = _mvn_logdensity(E, Sigma)
    return None if ld is None else np.exp(ld)


def _posterior(log_num: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Posterior weights and log mixture density from log(pi_s) + log f_s, without underflow."""
    m = np.max(log_num, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        e = np.exp(log_num - m)
        tot = e.sum(axis=1, keepdims=True)
        return e / tot, (m + np.log(tot)).ravel()


def _profile_core(class_settings, classProb, apollo_inputs, em_settings, logs: bool = False) -> dict | None:
    """Fixed point for the class covariances, computed in logs.

    Returns the mixture likelihood as ``logP`` and ``P``. In density mode (``logs=False``) a
    mixture density that underflows to zero makes the whole likelihood zero, as in R.
    """
    S = len(class_settings)
    N = len(apollo_inputs.database)
    E = [mtuem_residuals(cs, apollo_inputs) for cs in class_settings]
    if any(not np.all(np.isfinite(e)) for e in E):
        return None
    M = E[0].shape[1]
    min_weight = em_settings.get("minWeight")
    if min_weight is None:
        min_weight = 2 * M
    Pi = np.column_stack([np.broadcast_to(np.asarray(p, dtype=float), (N,)) for p in classProb])
    class_logs = getattr(classProb, "log", None)
    with np.errstate(divide="ignore"):
        log_Pi = (np.column_stack([np.broadcast_to(np.asarray(lp, dtype=float), (N,)) for lp in class_logs])
                  if class_logs is not None else np.log(Pi))

    # Deterministic unweighted start, recomputed on every call so the objective is path
    # independent (warm starts would break numerical gradients).
    Sigma = [cov_complete(e) for e in E]
    logdens = np.empty((N, S))

    def mixture():
        for s in range(S):
            ld = _mvn_logdensity(E[s], Sigma[s])
            if ld is None:
                return None
            logdens[:, s] = ld
        W, logP = _posterior(log_Pi + logdens)
        if np.any(~np.isfinite(logP)) or (not logs and np.any(np.exp(logP) <= 0)):
            return None
        return W, logP

    for _ in range(int(em_settings["maxIterations"])):
        out = mixture()
        if out is None:
            return None
        W, _ = out
        if np.any(W.sum(axis=0) < min_weight):
            return None
        Sigma_new = []
        for s in range(S):
            ws = W[:, s]
            Es = E[s] * np.sqrt(ws)[:, None]
            Sigma_new.append(Es.T @ Es / ws.sum())
        delta = max(np.max(np.abs(a - b)) for a, b in zip(Sigma_new, Sigma))
        Sigma = Sigma_new
        if delta < em_settings["tol"]:
            break

    out = mixture()
    if out is None:
        return None
    W, logP = out
    return {"P": np.exp(logP), "logP": logP, "Sigma": Sigma, "W": W, "Pi": Pi}


def _stash(core: dict) -> None:
    _PROFILE_STASH.update(Sigma=core["Sigma"], W=core["W"], Pi=core["Pi"])


def mtuem_lc_profiled_likelihood(lc_settings: Mapping, functionality: str = "estimate", apollo_inputs=None):
    """Latent class Theta-Phi likelihood with the covariance concentrated out (inner fixed point).

    ``lc_settings`` holds ``class_settings`` (one ``mtuem_settings`` per class, without ``sig``,
    ``rho`` or ``cholesky``), ``classProb`` (from the engine's ``class_alloc``) and optionally
    ``em_settings`` with ``tol`` (1e-10), ``maxIterations`` (50) and ``minWeight`` (2 * M).
    The result is already the class mixture, so do not pass it through the engine's ``lc``.
    """
    em = {"tol": 1e-10, "maxIterations": 50, "minWeight": None}
    em.update(lc_settings.get("em_settings") or {})
    class_settings = lc_settings["class_settings"]
    classProb = lc_settings["classProb"]
    S = len(class_settings)

    if functionality == "preprocess":
        return {"componentName": lc_settings.get("componentName") or "..", "gradient": False}
    N = len(apollo_inputs.database)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(apollo_inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    logs = works_in_logs(apollo_inputs)
    if functionality in ("estimate", "conditionals", "raw", "components"):
        core = _profile_core(class_settings, classProb, apollo_inputs, em, logs)
        if core is None:
            return np.full(N, -np.inf) if logs else np.zeros(N)
        _stash(core)
        return core["logP"] if logs else core["P"]

    if functionality in ("output", "report"):
        ans = mtuem_lc_profiled_likelihood(lc_settings, "estimate", apollo_inputs)
        if functionality == "report":
            from .utils import get_tw_thph, get_values_of_time_thph

            for s, cs in enumerate(class_settings, start=1):
                tw = get_tw_thph(cs["work_elasticities"], cs["tau"], cs["Tc"], cs["Ec"], cs["w"])
                vot = get_values_of_time_thph(tw, cs["work_elasticities"], cs["Tc"], cs["Ec"], cs["w"]).mean(axis=0)
                print(f"Class {s} VoL: {vot[0]:g} \nClass {s} VTAW: {vot[1]:g} ")
        return ans

    if functionality == "prediction":
        core = _profile_core(class_settings, classProb, apollo_inputs, em, works_in_logs(apollo_inputs))
        if core is None:
            return np.zeros(N)
        _stash(core)
        frames = []
        for s in range(S):
            opt = mtuem_likelihood(class_settings[s], "prediction", apollo_inputs)
            frames.append(opt.add_suffix(f"_class{s + 1}"))
        post = pd.DataFrame(core["W"], columns=[f"post_class{s + 1}" for s in range(S)])
        return pd.concat(frames + [post], axis=1)

    raise ValueError("Invalid value of argument 'functionality'")


def mtuem_get_profiled_sigma() -> dict:
    """Covariances, sds, correlations, posteriors and shares from the LAST profiled evaluation.

    Re-evaluate ``apollo_probabilities`` at ``model.estimate`` first to read them at the optimum.
    """
    if "Sigma" not in _PROFILE_STASH:
        raise RuntimeError(
            "No profiled covariance is available yet. Call apollo_probabilities (which in turn calls "
            "mtuem_lc_profiled_likelihood) at least once before calling mtuem_get_profiled_sigma()."
        )
    Sigma = _PROFILE_STASH["Sigma"]
    sig, rho = zip(*(cov_to_sig_rho(s) for s in Sigma))
    return {
        "Sigma": list(Sigma),
        "sig": list(sig),
        "rho": list(rho),
        "posterior": _PROFILE_STASH["W"],
        "shares": _PROFILE_STASH["Pi"].mean(axis=0),
    }


def mtuem_profile_ic(model, nClass: int, nEq: int) -> dict:
    """AIC/BIC counting the profiled-out covariance parameters, nClass * nEq * (nEq + 1) / 2."""
    k = len(model.estimate) - len(model.fixed) + nClass * nEq * (nEq + 1) / 2
    ll = model.maximum
    n = model.nObs
    return {"k": k, "AIC": -2 * ll + 2 * k, "BIC": -2 * ll + k * np.log(n)}
