"""MTUEM likelihoods (Theta-Phi, alpha-beta and additive quadratic) and their residuals.

Each likelihood takes the same ``settings`` dictionary as the corresponding R function
(``theta_phi``: ``mtuem_likelihood``, ``alpha_beta``: ``mtuem_ab_likelihood``,
``additive_quadratic``: ``additive_quadratic_mtuem_likelihood``), a ``functionality`` string,
and the model ``inputs`` (in R they are read from the calling frame; here they are passed
explicitly). ``mtuem.likelihood_jax`` holds the JAX versions.

For latent class models, ``class_alloc`` (class allocation probabilities), ``lc`` (the mixture
of class log-likelihoods) and ``lc_profiled`` (the profiled latent class likelihood) are the
NumPy counterparts of the functions of the same name in ``mtuem.likelihood_jax``.
"""

from collections.abc import Mapping

import numpy as np
import pandas as pd

from ._engine import class_alloc, lc
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

__all__ = [
    "additive_quadratic",
    "alpha_beta",
    "class_alloc",
    "lc",
    "lc_profiled",
    "residuals",
    "theta_phi",
    "works_in_logs",
]

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
_AQ_DEFAULTS = {"times_satiations": [], "goods_satiations": [], "omitted_times": [], "omitted_goods": []}


def _with_defaults(settings: Mapping, extra: Mapping | None = None) -> dict:
    out = dict(_DEFAULTS)
    if extra:
        out.update(extra)
    out.update(settings)
    return out


def works_in_logs(inputs) -> bool:
    """True when the model works in logs (always the case for inputs made by the engine).

    Then likelihood functions return log-likelihoods, which never underflow. Inputs without
    the flag (for example a bare namespace with a ``database``) get densities, as in R.
    """
    return bool(getattr(inputs, "work_in_logs", False))


def _zero(N: int, logs: bool) -> np.ndarray:
    """A zero likelihood: -inf in logs."""
    return np.full(N, -np.inf) if logs else np.zeros(N)


def _n_rows(inputs) -> int:
    return len(inputs.database)


def _observed(inputs, columns: list[str]) -> np.ndarray:
    obs_matrix = getattr(inputs, "obs_matrix", None)
    if obs_matrix is not None:
        return np.asarray(obs_matrix, dtype=float)
    if hasattr(inputs, "matrix"):
        return inputs.matrix(columns)
    return inputs.database[columns].to_numpy(dtype=float)


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
def residuals(settings: Mapping, inputs) -> np.ndarray:
    """Observed minus optimal allocations, columns ``work_times + free_times + free_goods``."""
    s = _with_defaults(settings)
    work_times = as_names(s["work_times"])
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    we = s["work_elasticities"]
    tw_opt = get_tw_thph(we, s["tau"], s["Tc"], s["Ec"], s["w"])
    tw_used = tw_opt if s["optimal_tw"] else _observed(inputs, work_times[:1]).ravel()
    ti_opt = get_ti_thph(s["times_elasticities"], we["Theta"], tw_used, s["tau"], s["Tc"])
    xj_opt = get_xi_thph(s["goods_elasticities"], s["goods_cost"], we["Phi"], tw_used, s["Ec"], s["w"])
    opt = np.column_stack([tw_opt, ti_opt, xj_opt])
    cols = work_times + free_times + free_goods
    return _observed(inputs, cols) - opt


# ----------------------------------------------------------------------------------------------
# Theta-Phi likelihood
# ----------------------------------------------------------------------------------------------
def theta_phi(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Likelihood of the Theta-Phi MTUEM (Jara-Diaz et al., 2008).

    Mirrors ``mtuem::mtuem_likelihood``. ``work_elasticities`` must hold ``Theta``, ``Phi`` and
    ``thw``. The error structure is profiled (no ``sig``/``rho``), uses explicit ``sig`` with
    profiled correlations, explicit ``sig`` and ``rho``, or a ``cholesky`` factor.

    When the model works in logs (always the case for inputs created by the engine's
    ``make_inputs``) the estimation functionalities return log-likelihoods; otherwise
    densities, as in R. The same holds for ``alpha_beta`` and ``additive_quadratic``.
    """
    s = _with_defaults(settings)
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    work_times = as_names(s["work_times"])
    flag_times = len(as_vector(s["times_elasticities"])) > 0
    flag_goods = len(as_vector(s["goods_elasticities"])) > 0
    M = 1 + len(free_times) + len(free_goods)
    if len(as_vector(s["cholesky"])) > 0:
        cholesky_to_sigma(s["cholesky"], M)  # validates the length early, as R does

    if functionality == "preprocess":
        return _preprocess(settings)
    N = _n_rows(inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    we = s["work_elasticities"]
    tau, Tc, Ec, w = s["tau"], s["Tc"], s["Ec"], s["w"]

    if functionality in ESTIMATE_LIKE:
        err = residuals(s, inputs)
        return _likelihood_from_errors(err, s, not (flag_times or flag_goods), M, works_in_logs(inputs))

    if functionality == "output":
        return theta_phi(settings, "estimate", inputs)

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

    raise ValueError("Invalid value of argument 'functionality'")


# ----------------------------------------------------------------------------------------------
# Alpha-Beta likelihood
# ----------------------------------------------------------------------------------------------
def alpha_beta(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Likelihood of the alpha-beta MTUEM. Mirrors ``mtuem::mtuem_ab_likelihood``.

    Besides the usual functionalities it accepts ``"get_covar"``, returning the covariance,
    correlation and standard deviations of the residuals.
    """
    s = _with_defaults(settings)
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
        return _preprocess(settings)
    N = _n_rows(inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(inputs) else np.ones(N)
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
        err = _observed(inputs, cols) - opt
        if functionality == "get_covar":
            cov = cov_complete(err)
            return {
                "covar": pd.DataFrame(cov, index=cols, columns=cols),
                "corr": pd.DataFrame(cor_complete(err), index=cols, columns=cols),
                "sigma": pd.Series(np.sqrt(np.diag(cov)), index=cols),
            }
        return _likelihood_from_errors(err, s, not (flag_times or flag_goods), M, works_in_logs(inputs))

    if functionality == "output":
        return alpha_beta(settings, "estimate", inputs)

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

    raise ValueError("Invalid value of argument 'functionality'")


# ----------------------------------------------------------------------------------------------
# Additive quadratic likelihood
# ----------------------------------------------------------------------------------------------
def _additive_quadratic_equations(s: Mapping) -> tuple[list[str], list[str]]:
    """All allocation columns and the ones kept as equations (``omitted_*`` removed), checked."""
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
    return all_cols, eq_names


def additive_quadratic(settings: Mapping, functionality: str = "estimate", inputs=None):
    """Likelihood of the additive quadratic MTUEM. Mirrors ``mtuem::additive_quadratic_mtuem_likelihood``.

    ``omitted_times`` / ``omitted_goods`` drop equations from the likelihood while their
    parameters stay in the aggregates S_theta, B, S_phi and H.
    """
    s = _with_defaults(settings, _AQ_DEFAULTS)
    work_times = as_names(s["work_times"])
    free_times = as_names(s["free_times"])
    free_goods = as_names(s["free_goods"])
    all_cols, eq_names = _additive_quadratic_equations(s)
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
        return _preprocess(settings)
    N = _n_rows(inputs)
    if functionality == "validate":
        return np.zeros(N) if works_in_logs(inputs) else np.ones(N)
    if functionality == "zero_LL":
        return np.full(N, np.nan)

    tau, Tc, Ec, w = s["tau"], s["Tc"], s["Ec"], s["w"]
    tw_opt = get_tw_additive_quadratic(s["work_elasticity"], s["work_satiation"], Sth, Sph, B, H, tau, Tc, Ec, w)

    if functionality in ESTIMATE_LIKE:
        tw_used = tw_opt if s["optimal_tw"] else inputs.database[work_times[0]].to_numpy(dtype=float)
        ti = get_ti_additive_quadratic(s["times_elasticities"], s["times_satiations"], Sth, B, tau, tw_used, Tc)
        xj = get_xi_additive_quadratic(
            s["goods_elasticities"], s["goods_satiations"], s["goods_cost"], Sph, H, tw_used, Ec, w
        )
        opt = pd.DataFrame(np.column_stack([tw_opt, ti, xj]), columns=all_cols)[eq_names].to_numpy()
        err = _observed(inputs, eq_names) - opt
        return _likelihood_from_errors(err, s, M == 1, M, works_in_logs(inputs))

    if functionality == "output":
        return additive_quadratic(settings, "estimate", inputs)

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

    raise ValueError("Invalid value of argument 'functionality'")


# ----------------------------------------------------------------------------------------------
# Profiled latent class likelihood
# ----------------------------------------------------------------------------------------------
def lc_profiled(lc_settings: Mapping, functionality: str = "estimate", inputs=None):
    """Latent class Theta-Phi likelihood with the class covariances profiled out.

    The same function as ``mtuem.lc_profiled.lc_profiled`` (see there for the settings), here
    beside the other NumPy likelihoods; ``mtuem.likelihood_jax.lc_profiled`` is its JAX version.
    """
    # Imported here: mtuem.lc_profiled builds on the residuals of this module.
    from .lc_profiled import lc_profiled as profiled

    return profiled(lc_settings, functionality, inputs)
