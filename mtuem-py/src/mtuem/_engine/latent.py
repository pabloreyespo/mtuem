"""Latent class helpers: class allocation, the class mixture and posterior class probabilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .inputs import Inputs

if TYPE_CHECKING:
    from .estimate import Estimate

ProbFn = Callable[..., dict]


class ClassProbs(list):
    """Class probabilities (a list) that also carry exact logs.

    ``lc`` and the posterior computations use ``.log``, so a class whose probability underflows
    to 0 keeps a finite log-probability.
    """

    def __init__(self, probs, log):
        super().__init__(probs)
        self.log = log


def class_alloc(classAlloc_settings: Mapping) -> ClassProbs:
    """Class allocation probabilities from utilities (a multinomial logit)."""
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


def _first_row(x: np.ndarray, inputs: Inputs) -> np.ndarray:
    x = np.broadcast_to(np.asarray(x, dtype=float), (inputs.nObs,))
    return x[inputs.first_rows] if inputs.panel else np.asarray(x)


def _class_components(in_class: Any) -> list:
    if isinstance(in_class, Mapping):
        return [v for k, v in in_class.items() if k != "model"]
    return list(in_class)


def lc(lc_settings: Mapping, inputs: Inputs, functionality: str = "estimate"):
    """Latent class mixture of in-class log-likelihoods (follows ``apollo_lc``).

    ``lc_settings`` holds ``inClassProb`` (dict or list of class log-likelihoods) and
    ``classProb`` (list of class allocation probabilities). With panel data each class
    log-likelihood is summed over the rows of an individual before mixing.
    """
    in_class = _class_components(lc_settings["inClassProb"])
    class_prob = lc_settings["classProb"]
    pis = list(class_prob)
    if len(in_class) != len(pis):
        raise ValueError("inClassProb and classProb must have the same number of classes")
    n = inputs.nObs
    if functionality == "preprocess":
        return {"componentName": lc_settings.get("componentName", "model")}
    if functionality == "validate":
        return np.zeros(n)
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
    # Mixture in logs: log sum_s exp(log pi_s + sum_rows log L_s), one value per individual.
    terms = [_first_row(lp, inputs) + _loglik_per_indiv(p, inputs)
             for p, lp in zip(in_class, _log_class_probs(class_prob))]
    return _log_posterior(np.column_stack(terms), np.zeros((inputs.nIndivs, len(terms))))[1]


def _loglik_per_indiv(L, inputs: Inputs) -> np.ndarray:
    """Log-likelihood per individual: the rows of a panel are summed."""
    logL = np.asarray(L, dtype=float).ravel()
    if logL.size == inputs.nObs and inputs.panel:
        return np.bincount(inputs.indiv_index, weights=logL, minlength=inputs.nIndivs)
    return logL


def _class_parts(beta: Mapping, prob_fn: ProbFn, inputs: Inputs, functionality: str = "conditionals"):
    """Per-individual class LOG-likelihoods and log class allocation probabilities at ``beta``.

    Working in logs keeps posteriors and EM weights accurate when class likelihoods are tiny
    (for example with a participation penalty), for both the numpy and the autodiff backend.
    """
    if inputs.lc_pars is None:
        raise ValueError("lc_pars is required for latent class functions")
    pis = inputs.attach(beta).pi_values
    P = prob_fn(dict(beta), inputs, functionality)
    comps = [v for k, v in P.items() if k != "model"]
    if len(comps) != len(pis):
        raise ValueError(
            "Could not find one component per class in the dict returned by the probabilities function. "
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


def lc_conditionals(est: Estimate, prob_fn: ProbFn, inputs: Inputs) -> pd.DataFrame:
    """Posterior class membership probabilities (follows ``apollo_lcConditionals``)."""
    post, _ = _log_posterior(*_class_parts(est.estimate.to_dict(), prob_fn, inputs))
    ids = inputs.ids[inputs.first_rows] if inputs.panel else inputs.ids
    out = pd.DataFrame(post, columns=[f"X{s + 1}" for s in range(post.shape[1])])
    out.insert(0, "ID", ids)
    return out
