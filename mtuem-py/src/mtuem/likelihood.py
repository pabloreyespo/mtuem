"""MTUEM likelihoods (Theta-Phi, alpha-beta and additive quadratic) and their residuals.

Each likelihood takes the same ``mtuem_settings`` dictionary as the R function of the same
name, a ``functionality`` string, and the ``apollo_inputs`` object (in R it is read from the
calling frame; here it is passed explicitly).
"""

from collections.abc import Mapping

import numpy as np
import pandas as pd

from .utils import (
    as_names,
    as_vector,
    cholesky_to_sigma,
    cor_complete,
    cov_complete,
    cov_to_sig_rho,
    get_additive_quadratic_essentials,
    get_ti_additive_quadratic,
    get_ti_albe,
    get_ti_thph,
    get_tw_additive_quadratic,
    get_tw_albe,
    get_tw_thph,
    get_values_of_time_additive_quadratic,
    get_values_of_time_albe,
    get_values_of_time_thph,
    get_xi_additive_quadratic,
    get_xi_albe,
    get_xi_thph,
    mvn_logpdf_cond,
    near_pd,
    rho_from_upper,
)

ESTIMATE_LIKE = ("estimate", "conditionals", "raw", "components")

_DEFAULTS = {
    "free_times": [],
    "free_goods": [],
    "goods_cost": 1,
    "times_elasticities": [],
    "goods_elasticities": [],
    "sig": [],
    "rho": [],
    "cholesky": [],
    "optimal_tw": True,
}


def _with_defaults(settings: Mapping, extra: Mapping | None = None) -> dict:
    out = dict(_DEFAULTS)
    if extra:
        out.update(extra)
    out.update(settings)
    return out


def works_in_logs(apollo_inputs) -> bool:
    """True when the model works in logs (``apollo_control$workInLogs``, default in mtuem-py).

    Then likelihood functions return log-likelihoods, which never underflow. Inputs without
    the flag (for example a bare namespace with a ``database``) get densities, as in R.
    """
    return bool(getattr(apollo_inputs, "work_in_logs", False))


def _zero(N: int, logs: bool) -> np.ndarray:
    """A zero likelihood: -inf in logs."""
    return np.full(N, -np.inf) if logs else np.zeros(N)


def _n_rows(apollo_inputs) -> int:
    return len(apollo_inputs.database)


def _observed(apollo_inputs, columns: list[str]) -> np.ndarray:
    obs_matrix = getattr(apollo_inputs, "obs_matrix", None)
    if obs_matrix is not None:
        return np.asarray(obs_matrix, dtype=float)
    if hasattr(apollo_inputs, "matrix"):
        return apollo_inputs.matrix(columns)
    return apollo_inputs.database[columns].to_numpy(dtype=float)


def _preprocess(settings: Mapping) -> dict:
    name = settings.get("componentName") or settings.get("componentName2") or ".."
    return {"componentName": name, "gradient": False}


def _error_structure(s: dict, M: int):
    """Resolve sig/rho/cholesky into (sig, rho, estimate_sig, estimate_rho)."""
    sig = as_vector(s["sig"])
    rho_vals = as_vector(s["rho"])
    chol = as_vector(s["cholesky"])
    estimate_sig = sig.size > 0
    estimate_rho = rho_vals.size > 0
    rho = None
    if chol.size > 0:
        sig, rho = cov_to_sig_rho(cholesky_to_sigma(chol, M))
        estimate_sig = estimate_rho = True
    elif estimate_rho:
        rho = rho_from_upper(rho_vals, sig.size)
    return sig, rho, estimate_sig, estimate_rho


def _likelihood_from_errors(err: np.ndarray, s: dict, single_equation: bool, M: int, logs: bool = False) -> np.ndarray:
    """Normal likelihood (or log-likelihood when ``logs``) of the residuals under the error structure."""
    N = err.shape[0]
    sig, rho, estimate_sig, estimate_rho = _error_structure(s, M)
    if not estimate_sig:
        sig = np.sqrt(np.diag(cov_complete(err)))

    if single_equation:
        if np.any(~np.isfinite(sig)) or np.any(sig <= 1e-6):
            return _zero(N, logs)
        mu = err / sig
        ll = (-0.5 * mu**2 - np.log(sig) - 0.5 * np.log(2 * np.pi)).sum(axis=1)
        return ll if logs else np.exp(ll)

    if not estimate_rho:
        if estimate_sig:
            rho = np.eye(sig.size)
        else:
            rho = cor_complete(err)
            if np.any(~np.isfinite(rho)):
                return _zero(N, logs)
            rho = near_pd(rho, corr=True, keep_diag=True)
    ll = mvn_logpdf_cond(err, sig, rho)
    if ll is None:
        return _zero(N, logs)
    return ll if logs else np.exp(ll)


def _allocation_frame(opt: dict[str, np.ndarray], vot: np.ndarray) -> pd.DataFrame:
    df = pd.DataFrame({k: np.asarray(v, dtype=float).ravel() for k, v in opt.items()})
    df["VoL"] = vot[:, 0]
    df["VTAW"] = vot[:, 1]
    return df


# ----------------------------------------------------------------------------------------------
# Residuals (Theta-Phi)
# ----------------------------------------------------------------------------------------------
def mtuem_residuals(mtuem_settings: Mapping, apollo_inputs) -> np.ndarray:
    """Observed minus optimal allocations, columns ``work_times + free_times + free_goods``."""
    s = _with_defaults(mtuem_settings)
    work_times = as_names(s["work_times"])
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    we = s["work_elasticities"]
    tw_opt = get_tw_thph(we, s["tau"], s["Tc"], s["Ec"], s["w"])
    tw_used = tw_opt if s["optimal_tw"] else _observed(apollo_inputs, work_times[:1]).ravel()
    ti_opt = get_ti_thph(s["times_elasticities"], we["Theta"], tw_used, s["tau"], s["Tc"])
    xj_opt = get_xi_thph(s["goods_elasticities"], s["goods_cost"], we["Phi"], tw_used, s["Ec"], s["w"])
    opt = np.column_stack([tw_opt, ti_opt, xj_opt])
    cols = work_times + free_times + free_goods
    return _observed(apollo_inputs, cols) - opt


# ----------------------------------------------------------------------------------------------
# Theta-Phi likelihood
# ----------------------------------------------------------------------------------------------
def mtuem_likelihood(mtuem_settings: Mapping, functionality: str = "estimate", apollo_inputs=None):
    """Likelihood of the Theta-Phi MTUEM (Jara-Diaz et al., 2008).

    Mirrors ``mtuem::mtuem_likelihood``. ``work_elasticities`` must hold ``Theta``, ``Phi`` and
    ``thw``. The error structure is profiled (no ``sig``/``rho``), uses explicit ``sig`` with
    profiled correlations, explicit ``sig`` and ``rho``, or a ``cholesky`` factor.

    When the model works in logs (``workInLogs``, the default for inputs created by
    ``apollo.validate_inputs``) the estimation functionalities return log-likelihoods; otherwise
    densities, as in R. The same holds for ``mtuem_ab_likelihood`` and
    ``additive_quadratic_mtuem_likelihood``.
    """
    s = _with_defaults(mtuem_settings)
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    work_times = as_names(s["work_times"])
    flag_times = len(as_vector(s["times_elasticities"])) > 0
    flag_goods = len(as_vector(s["goods_elasticities"])) > 0
    M = 1 + len(free_times) + len(free_goods)
    if len(as_vector(s["cholesky"])) > 0:
        cholesky_to_sigma(s["cholesky"], M)  # validates the length early, as R does

    if functionality == "preprocess":
        return _preprocess(mtuem_settings)
    N = _n_rows(apollo_inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(apollo_inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    we = s["work_elasticities"]
    tau, Tc, Ec, w = s["tau"], s["Tc"], s["Ec"], s["w"]

    if functionality in ESTIMATE_LIKE:
        err = mtuem_residuals(s, apollo_inputs)
        return _likelihood_from_errors(err, s, not (flag_times or flag_goods), M, works_in_logs(apollo_inputs))

    if functionality in ("output", "report"):
        ans = mtuem_likelihood(mtuem_settings, "estimate", apollo_inputs)
        if functionality == "report":
            tw_opt = get_tw_thph(we, tau, Tc, Ec, w)
            vot = get_values_of_time_thph(tw_opt, we, Tc, Ec, w).mean(axis=0)
            print(f"VoL: {vot[0]:g} \nVTAW: {vot[1]:g} ")
        return ans

    if functionality == "prediction":
        tw_opt = get_tw_thph(we, tau, Tc, Ec, w)
        opt = {work_times[0]: tw_opt}
        if flag_times:
            ti = get_ti_thph(s["times_elasticities"], we["Theta"], tw_opt, tau, Tc)
            opt.update({n: ti[:, k] for k, n in enumerate(free_times)})
            opt["Tfi"] = tau - np.asarray(Tc) - ti.sum(axis=1) - tw_opt.ravel()
        else:
            opt["T"] = tau - np.asarray(Tc) - tw_opt.ravel()
        if flag_goods:
            xj = get_xi_thph(s["goods_elasticities"], s["goods_cost"], we["Phi"], tw_opt, Ec, w)
            opt.update({n: xj[:, k] for k, n in enumerate(free_goods)})
            opt["Xfj"] = np.asarray(w) * tw_opt.ravel() - np.asarray(Ec) - xj.sum(axis=1)
        else:
            opt["X"] = np.asarray(w) * tw_opt.ravel() - np.asarray(Ec)
        return _allocation_frame(opt, get_values_of_time_thph(tw_opt, we, Tc, Ec, w))

    if functionality == "gradient":
        return {"like": None, "grad": None}
    raise ValueError("Invalid value of argument 'functionality'")


# ----------------------------------------------------------------------------------------------
# Alpha-Beta likelihood
# ----------------------------------------------------------------------------------------------
def mtuem_ab_likelihood(mtuem_settings: Mapping, functionality: str = "estimate", apollo_inputs=None):
    """Likelihood of the alpha-beta MTUEM. Mirrors ``mtuem::mtuem_ab_likelihood``.

    Besides the usual functionalities it accepts ``"get_covar"``, returning the covariance,
    correlation and standard deviations of the residuals.
    """
    s = _with_defaults(mtuem_settings)
    s.pop("optimal_tw", None)
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    work_times = as_names(s["work_times"])
    flag_times = len(as_vector(s["times_elasticities"])) > 0
    flag_goods = len(as_vector(s["goods_elasticities"])) > 0
    M = 1 + len(free_times) + len(free_goods)
    if len(as_vector(s["cholesky"])) > 0:
        cholesky_to_sigma(s["cholesky"], M)

    if functionality == "preprocess":
        return _preprocess(mtuem_settings)
    N = _n_rows(apollo_inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(apollo_inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    we = s["work_elasticities"]
    tau, Tc, Ec, w = s["tau"], s["Tc"], s["Ec"], s["w"]

    if functionality in ESTIMATE_LIKE + ("get_covar",):
        tw_opt = get_tw_albe(we, tau, Tc, Ec, w)
        ti = get_ti_albe(s["times_elasticities"], we["beta"], tw_opt, tau, Tc)
        xj = get_xi_albe(s["goods_elasticities"], s["goods_cost"], we["alpha"], tw_opt, Ec, w)
        opt = np.column_stack([tw_opt, ti, xj])
        cols = work_times + free_times + free_goods
        err = _observed(apollo_inputs, cols) - opt
        if functionality == "get_covar":
            cov = cov_complete(err)
            return {
                "covar": pd.DataFrame(cov, index=cols, columns=cols),
                "corr": pd.DataFrame(cor_complete(err), index=cols, columns=cols),
                "sigma": pd.Series(np.sqrt(np.diag(cov)), index=cols),
            }
        return _likelihood_from_errors(err, s, not (flag_times or flag_goods), M, works_in_logs(apollo_inputs))

    if functionality in ("output", "report"):
        ans = mtuem_ab_likelihood(mtuem_settings, "estimate", apollo_inputs)
        if functionality == "report":
            tw_opt = get_tw_albe(we, tau, Tc, Ec, w)
            vot = get_values_of_time_albe(tw_opt, we, Tc, Ec, w).mean(axis=0)
            print(f"VoL: {vot[0]:g} \nVTAW: {vot[1]:g} ")
        return ans

    if functionality == "prediction":
        tw_opt = get_tw_albe(we, tau, Tc, Ec, w)
        opt = {work_times[0]: tw_opt}
        if flag_times:
            ti = get_ti_albe(s["times_elasticities"], we["beta"], tw_opt, tau, Tc)
            opt.update({n: ti[:, k] for k, n in enumerate(free_times)})
            opt["Tfi"] = tau - np.asarray(Tc) - ti.sum(axis=1) - tw_opt.ravel()
        else:
            opt["T"] = tau - np.asarray(Tc) - tw_opt.ravel()
        if flag_goods:
            xj = get_xi_albe(s["goods_elasticities"], s["goods_cost"], we["alpha"], tw_opt, Ec, w)
            opt.update({n: xj[:, k] for k, n in enumerate(free_goods)})
            opt["Xfj"] = np.asarray(w) * tw_opt.ravel() - np.asarray(Ec) - xj.sum(axis=1)
        else:
            opt["X"] = np.asarray(w) * tw_opt.ravel() - np.asarray(Ec)
        return _allocation_frame(opt, get_values_of_time_albe(tw_opt, we, Tc, Ec, w))

    if functionality == "gradient":
        return {"like": None, "grad": None}
    raise ValueError("Invalid value of argument 'functionality'")


# ----------------------------------------------------------------------------------------------
# Additive quadratic likelihood
# ----------------------------------------------------------------------------------------------
def additive_quadratic_mtuem_likelihood(
    mtuem_settings: Mapping, functionality: str = "estimate", apollo_inputs=None
):
    """Likelihood of the additive quadratic MTUEM. Mirrors the R function of the same name.

    ``omitted_times`` / ``omitted_goods`` drop equations from the likelihood while their
    parameters stay in the aggregates S_theta, B, S_phi and H.
    """
    s = _with_defaults(
        mtuem_settings,
        {"times_satiations": [], "goods_satiations": [], "omitted_times": [], "omitted_goods": []},
    )
    work_times = as_names(s["work_times"])
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    omitted = as_names(s["omitted_times"]) + as_names(s["omitted_goods"])
    bad = [o for o in omitted if o not in free_times + free_goods]
    if bad:
        raise ValueError(
            "omitted_times / omitted_goods name categories that are not in free_times / free_goods: "
            + ", ".join(bad)
            + "."
        )
    all_cols = work_times + free_times + free_goods
    eq_names = [c for c in all_cols if c not in omitted]
    if not eq_names:
        raise ValueError("Every equation has been omitted; nothing left to estimate.")
    flag_times = len(as_vector(s["times_elasticities"])) > 0
    flag_goods = len(as_vector(s["goods_elasticities"])) > 0
    M = len(eq_names)
    if len(as_vector(s["cholesky"])) > 0:
        cholesky_to_sigma(s["cholesky"], M)

    agg = get_additive_quadratic_essentials(
        s["times_elasticities"], s["times_satiations"], s["goods_elasticities"], s["goods_satiations"], s["goods_cost"]
    )
    Sth, Sph, B, H = agg["Sth"], agg["Sph"], agg["B"], agg["H"]

    if functionality == "preprocess":
        return _preprocess(mtuem_settings)
    N = _n_rows(apollo_inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(apollo_inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    tau, Tc, Ec, w = s["tau"], s["Tc"], s["Ec"], s["w"]
    tw_opt = get_tw_additive_quadratic(s["work_elasticity"], s["work_satiation"], Sth, Sph, B, H, tau, Tc, Ec, w)

    if functionality in ESTIMATE_LIKE:
        tw_used = tw_opt if s["optimal_tw"] else apollo_inputs.database[work_times[0]].to_numpy(dtype=float)
        ti = get_ti_additive_quadratic(s["times_elasticities"], s["times_satiations"], Sth, B, tau, tw_used, Tc)
        xj = get_xi_additive_quadratic(
            s["goods_elasticities"], s["goods_satiations"], s["goods_cost"], Sph, H, tw_used, Ec, w
        )
        opt = pd.DataFrame(np.column_stack([tw_opt, ti, xj]), columns=all_cols)[eq_names].to_numpy()
        err = _observed(apollo_inputs, eq_names) - opt
        return _likelihood_from_errors(err, s, M == 1, M, works_in_logs(apollo_inputs))

    if functionality in ("output", "report"):
        ans = additive_quadratic_mtuem_likelihood(mtuem_settings, "estimate", apollo_inputs)
        if functionality == "report":
            vot = get_values_of_time_additive_quadratic(Sth, Sph, B, H, tw_opt, tau, Tc, Ec, w).mean(axis=0)
            print(f"VoL: {vot[0]:g} \nVTAW: {vot[1]:g} ")
        return ans

    if functionality == "prediction":
        opt = {work_times[0]: tw_opt}
        if flag_times:
            ti = get_ti_additive_quadratic(s["times_elasticities"], s["times_satiations"], Sth, B, tau, tw_opt, Tc)
            opt.update({n: ti[:, k] for k, n in enumerate(free_times)})
        else:
            opt["T"] = tau - np.asarray(Tc) - tw_opt.ravel()
        if flag_goods:
            xj = get_xi_additive_quadratic(
                s["goods_elasticities"], s["goods_satiations"], s["goods_cost"], Sph, H, tw_opt, Ec, w
            )
            opt.update({n: xj[:, k] for k, n in enumerate(free_goods)})
        else:
            opt["X"] = np.asarray(w) * tw_opt.ravel() - np.asarray(Ec)
        vot = get_values_of_time_additive_quadratic(Sth, Sph, B, H, tw_opt, tau, Tc, Ec, w)
        return _allocation_frame(opt, vot)

    if functionality == "gradient":
        return {"like": None, "grad": None}
    raise ValueError("Invalid value of argument 'functionality'")
