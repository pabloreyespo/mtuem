"""Post-estimation helpers: error correlations, Cholesky decoding and values of time."""

from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy import stats

from . import apollo
from .lc_profiled import mtuem_get_profiled_sigma
from .utils import cov_complete, get_tw_thph


# ----------------------------------------------------------------------------------------------
# Error covariance
# ----------------------------------------------------------------------------------------------
def mtuem_get_corr(model, apollo_probabilities, apollo_inputs) -> dict:
    """Covariance, correlation and sds of observed minus predicted allocations."""
    pred = apollo.prediction(model, apollo_probabilities, apollo_inputs, {"silent": True})
    drop = ["ID", "Observation", "Tfi", "Xfj", "VoL", "VTAW"]
    opt = pred.drop(columns=[c for c in drop if c in pred.columns])
    err = apollo_inputs.database[list(opt.columns)].to_numpy(dtype=float) - opt.to_numpy(dtype=float)
    cols = list(opt.columns)
    cov = cov_complete(err)
    sd = np.sqrt(np.diag(cov))
    return {
        "covar": pd.DataFrame(cov, index=cols, columns=cols),
        "corr": pd.DataFrame(cov / np.outer(sd, sd), index=cols, columns=cols),
        "sigma": pd.Series(sd, index=cols),
    }


def cholesky_to_cov(estimates: Mapping, cholesky_names: list[str]) -> np.ndarray:
    """Rebuild Sigma = L L' from the named Cholesky entries (row-wise lower triangle)."""
    K = len(cholesky_names)
    M = int(round((-1 + np.sqrt(1 + 8 * K)) / 2))
    L = np.zeros((M, M))
    L[np.tril_indices(M)] = [float(estimates[n]) for n in cholesky_names]
    return L @ L.T


def mtuem_cholesky_deltaMethod(
    model, cholesky_names: list[str], eq_names: list[str] | None = None, deltaMethod_settings_extra: Mapping | None = None
) -> pd.DataFrame:
    """Standard deviations and correlations implied by the Cholesky factor, with delta method s.e."""
    K = len(cholesky_names)
    M = int(round((-1 + np.sqrt(1 + 8 * K)) / 2))
    if M * (M + 1) // 2 != K:
        raise ValueError("The number of Cholesky parameter names must be of the form M * (M + 1) / 2.")

    def lab(i):
        return eq_names[i] if eq_names is not None else str(i + 1)

    def sigma_of(b):
        return cholesky_to_cov(b, cholesky_names)

    expressions = {}
    for i in range(M):
        expressions[f"sig_{lab(i)}"] = lambda b, i=i: float(np.sqrt(sigma_of(b)[i, i]))
    for i in range(1, M):
        for j in range(i):
            def rho(b, i=i, j=j):
                S = sigma_of(b)
                return float(S[i, j] / np.sqrt(S[i, i] * S[j, j]))

            expressions[f"rho_{lab(j)}_{lab(i)}"] = rho
    settings = {"expression": expressions}
    settings.update(deltaMethod_settings_extra or {})
    return apollo.delta_method(model, settings)


# ----------------------------------------------------------------------------------------------
# Values of time
# ----------------------------------------------------------------------------------------------
def _vot_par(b: pd.Series, name: str) -> float:
    if name not in b.index:
        raise KeyError(
            f"Parameter '{name}' is not among the model estimates. Check nClass and the parameter naming: "
            "mtuem_values_of_time expects the Theta-Phi parameterization, with a '_s' suffix per class when nClass > 1."
        )
    return float(b[name])


def _wtd_se(x: np.ndarray, wt: np.ndarray) -> float:
    ok = ~np.isnan(x) & ~np.isnan(wt)
    x, wt = x[ok], wt[ok]
    m = np.average(x, weights=wt)
    v = np.sum(wt * (x - m) ** 2) / (np.sum(wt) - 1)
    with np.errstate(invalid="ignore"):
        return float(np.sqrt(v / np.sum(wt)))


def _wmean(x: np.ndarray, wt: np.ndarray) -> float:
    ok = ~np.isnan(x)
    return float(np.sum(x[ok] * wt[ok]) / np.sum(wt[ok]))


def _vot_weights(weights, model, apollo_probabilities, apollo_inputs, nClass: int, N: int) -> np.ndarray:
    if weights is None:
        if nClass == 1:
            return np.ones((N, 1))
        weights = _vot_posterior(model, apollo_probabilities, apollo_inputs, nClass)
    if isinstance(weights, pd.DataFrame):
        weights = weights.drop(columns=[c for c in ("ID", "Observation") if c in weights.columns]).to_numpy(float)
    elif isinstance(weights, (list, tuple)):
        weights = np.column_stack([np.asarray(w, dtype=float) for w in weights])
    weights = np.asarray(weights, dtype=float)
    if weights.ndim == 1:
        weights = np.repeat(weights[:, None], nClass, axis=1)
    if weights.shape != (N, nClass):
        raise ValueError(f"weights must have {N} rows and {nClass} columns, but has {weights.shape[0]} x {weights.shape[1]}.")
    return weights


def _vot_posterior(model, apollo_probabilities, apollo_inputs, nClass: int) -> np.ndarray:
    post = None
    if apollo_probabilities is not None:
        try:
            post = apollo.lc_conditionals(model, apollo_probabilities, apollo_inputs)
        except (ValueError, AttributeError, KeyError):
            post = None
    if post is None:
        try:
            if apollo_probabilities is not None:
                apollo_probabilities(model.estimate.to_dict(), apollo_inputs, "estimate")
            prof = mtuem_get_profiled_sigma()
            if prof["posterior"].shape[1] == nClass:
                post = prof["posterior"]
        except RuntimeError:
            post = None
    if post is None:
        raise ValueError(
            "Could not recover the posterior class probabilities. Pass apollo_probabilities so they can be "
            "computed with lc_conditionals, or supply them through the 'weights' argument."
        )
    return post


def _vot_report(summ: pd.DataFrame, nClass: int, conf_level: float) -> None:
    labels = {"VoL": "VoL              ", "VTAW": "VTAW             ", "wage": "observed wage    "}
    print(f"\nValues of time, {100 * conf_level:.0f}% confidence intervals")
    for s in range(1, nClass + 1):
        if nClass > 1:
            print(f"\nCLASS {s}:")
        for _, r in summ[summ["class"] == s].iterrows():
            print(f"  {labels.get(r['quantity'], r['quantity'])} {r['mean']:9.3f}  [ {r['lower']:9.3f} ; {r['upper']:9.3f} ]")
    print()


def mtuem_values_of_time(
    model,
    apollo_inputs,
    apollo_probabilities=None,
    nClass: int = 1,
    weights=None,
    tau: float = 168,
    Tc: str = "Tc",
    Ec: str = "Ec",
    w: str = "w",
    class_specific_budgets: bool = False,
    standard_errors: bool = True,
    conf_level: float = 0.95,
    silent: bool = False,
) -> dict:
    """Values of leisure (VoL) and of time assigned to work (VTAW) for a Theta-Phi MTUEM.

    The leading fraction of each value of time is averaged over the sample (weighted by class
    membership in latent class models) and the delta method is applied to the parameter ratio,
    exactly as ``mtuem::mtuem_values_of_time``.
    """
    db = apollo_inputs.database
    N = len(db)
    b = model.estimate
    W = _vot_weights(weights, model, apollo_probabilities, apollo_inputs, nClass, N)
    z = stats.norm.ppf(1 - (1 - conf_level) / 2)
    wage = db[w].to_numpy(dtype=float)
    summ, indiv, deltas = [], [], []
    for s in range(1, nClass + 1):
        suffix = f"_{s}" if nClass > 1 else ""
        names = {k: f"{k}{suffix}" for k in ("Theta", "Phi", "thw")}
        we = {k: _vot_par(b, n) for k, n in names.items()}
        Tc_s = db[f"{Tc}_{s}" if class_specific_budgets else Tc].to_numpy(dtype=float)
        Ec_s = db[f"{Ec}_{s}" if class_specific_budgets else Ec].to_numpy(dtype=float)
        Tw = get_tw_thph(we, tau, Tc_s, Ec_s, wage).ravel()
        cte_vol_i = (wage * Tw - Ec_s) / (tau - Tw - Tc_s)
        cte_vtaw_i = (wage * Tw - Ec_s) / Tw
        coef_vol = we["Theta"] / we["Phi"]
        coef_vtaw = we["thw"] / we["Phi"]
        wt = W[:, s - 1]
        cte_vol = _wmean(cte_vol_i, wt)
        cte_vtaw = _wmean(cte_vtaw_i, wt)
        indiv.append(pd.DataFrame({"class": s, "Tw": Tw, "VoL": coef_vol * cte_vol_i,
                                   "VTAW": coef_vtaw * cte_vtaw_i, "weight": wt}))
        if standard_errors:
            delta = apollo.delta_method(model, {"expression": {
                "VoL": lambda p, c=cte_vol, n=names: c * p[n["Theta"]] / p[n["Phi"]],
                "VTAW": lambda p, c=cte_vtaw, n=names: c * p[n["thw"]] / p[n["Phi"]],
            }}, silent=silent)
            deltas.append(delta)
            vol_se, vtaw_se = delta["s.e."].to_numpy(dtype=float)
        else:
            vol_se = vtaw_se = np.nan
        summ.append(pd.DataFrame({
            "class": s,
            "quantity": ["VoL", "VTAW", "wage"],
            "mean": [coef_vol * cte_vol, coef_vtaw * cte_vtaw, _wmean(wage, wt)],
            "s.e.": [vol_se, vtaw_se, _wtd_se(wage, wt)],
        }))
    summ = pd.concat(summ, ignore_index=True)
    summ["lower"] = summ["mean"] - z * summ["s.e."]
    summ["upper"] = summ["mean"] + z * summ["s.e."]
    if not silent:
        _vot_report(summ, nClass, conf_level)
    return {
        "summary": summ,
        "individual": pd.concat(indiv, ignore_index=True),
        "delta": deltas if standard_errors else None,
    }
