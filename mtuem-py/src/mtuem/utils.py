"""Closed-form optimal allocations, values of time and small linear algebra helpers.

Every function mirrors the exported helper of the same name in the R package (R/utils.R).
Elasticity arguments accept a dict, a list or a scalar, like the R lists they replace.
"""

from collections.abc import Mapping

import numpy as np


def as_vector(x) -> np.ndarray:
    """Flatten a dict, list, scalar or array of parameters, like R's ``unlist``."""
    if x is None:
        return np.zeros(0)
    if isinstance(x, Mapping):
        x = list(x.values())
    if isinstance(x, (list, tuple)):
        if len(x) == 0:
            return np.zeros(0)
        return np.concatenate([np.atleast_1d(np.asarray(v, dtype=float)).ravel() for v in x])
    return np.atleast_1d(np.asarray(x, dtype=float)).ravel()


def as_names(x) -> list[str]:
    """Flatten a column-name specification (None, str, list) into a list of strings."""
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return [str(v) for v in x]


def _col(x) -> np.ndarray:
    return np.asarray(x, dtype=float).reshape(-1, 1)


# ----------------------------------------------------------------------------------------------
# Theta-Phi parameterisation
# ----------------------------------------------------------------------------------------------
def get_tw_thph(work_elasticities: Mapping, tau, Tc, Ec, w) -> np.ndarray:
    """Optimal work time, Theta-Phi formulation. Returns an (N, 1) array."""
    Theta = work_elasticities["Theta"]
    Phi = work_elasticities["Phi"]
    thw = work_elasticities["thw"]
    Tc = np.asarray(Tc, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ec = np.asarray(Ec, dtype=float) / (np.asarray(w, dtype=float) * (tau - Tc))
    base = Theta + Phi + thw
    thetaphiec = Phi + thw + (Theta + thw) * ec
    with np.errstate(invalid="ignore"):
        sqrt_term = np.sqrt(thetaphiec**2 - 4 * thw * ec * base)
    tw_opt = 0.5 * (tau - Tc) * (thetaphiec + sqrt_term) / base
    return _col(tw_opt)


def get_ti_thph(times_elasticities, Theta, Tw, tau, Tc) -> np.ndarray:
    """Optimal free times, Theta-Phi formulation. Returns an (N, K) array."""
    th = as_vector(times_elasticities)
    left = np.asarray(tau - np.asarray(Tw, dtype=float).ravel() - np.asarray(Tc, dtype=float), dtype=float)
    return np.outer(left, th / Theta) if th.size else np.zeros((left.size, 0))


def get_xi_thph(goods_elasticities, goods_cost, Phi, Tw, Ec, w) -> np.ndarray:
    """Optimal free expenditures, Theta-Phi formulation. Returns an (N, J) array."""
    ph = as_vector(goods_elasticities)
    if ph.size == 0:
        return np.zeros((np.asarray(Tw).size, 0))
    cost = as_vector(goods_cost)
    income = np.asarray(w, dtype=float) * np.asarray(Tw, dtype=float).ravel() - np.asarray(Ec, dtype=float)
    return np.outer(income, (ph / cost) / Phi)


def get_values_of_time_thph(tw, work_elasticities: Mapping, Tc, Ec, w) -> np.ndarray:
    """VoL and VTAW per observation, Theta-Phi formulation. Returns (N, 2).

    As in the R package the leisure denominator uses a fixed 168-hour week.
    """
    tw = np.asarray(tw, dtype=float).ravel()
    w = np.asarray(w, dtype=float)
    cte_vol = (w * tw - np.asarray(Ec, dtype=float)) / (168 - tw - np.asarray(Tc, dtype=float))
    vol = cte_vol * work_elasticities["Theta"] / work_elasticities["Phi"]
    return np.column_stack([vol, vol - w])


# ----------------------------------------------------------------------------------------------
# Alpha-Beta parameterisation
# ----------------------------------------------------------------------------------------------
def get_tw_albe(work_elasticities: Mapping, tau, Tc, Ec, w) -> np.ndarray:
    """Optimal work time, alpha-beta formulation. Returns an (N, 1) array."""
    alpha = work_elasticities["alpha"]
    beta = work_elasticities["beta"]
    Tc = np.asarray(Tc, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ec = np.asarray(Ec, dtype=float) / (np.asarray(w, dtype=float) * (tau - Tc))
    betaalphaec = beta + alpha * ec
    with np.errstate(invalid="ignore"):
        sqrt_term = np.sqrt(betaalphaec**2 - ec * (2 * alpha + 2 * beta - 1))
    return _col((tau - Tc) * (betaalphaec + sqrt_term))


def get_ti_albe(times_elasticities, beta, Tw, tau, Tc) -> np.ndarray:
    th = as_vector(times_elasticities)
    left = np.asarray(tau - np.asarray(Tw, dtype=float).ravel() - np.asarray(Tc, dtype=float), dtype=float)
    return np.outer(left, th / (1 - 2 * beta)) if th.size else np.zeros((left.size, 0))


def get_xi_albe(goods_elasticities, goods_cost, alpha, Tw, Ec, w) -> np.ndarray:
    ph = as_vector(goods_elasticities)
    if ph.size == 0:
        return np.zeros((np.asarray(Tw).size, 0))
    cost = as_vector(goods_cost)
    income = np.asarray(w, dtype=float) * np.asarray(Tw, dtype=float).ravel() - np.asarray(Ec, dtype=float)
    return np.outer(income, (ph / cost) / (1 - alpha))


def get_values_of_time_albe(tw, work_elasticities: Mapping, Tc, Ec, w) -> np.ndarray:
    tw = np.asarray(tw, dtype=float).ravel()
    w = np.asarray(w, dtype=float)
    cte_vol = (w * tw - np.asarray(Ec, dtype=float)) / (168 - tw - np.asarray(Tc, dtype=float))
    coef = (1 - 2 * work_elasticities["beta"]) / (1 - 2 * work_elasticities["alpha"])
    vol = cte_vol * coef
    return np.column_stack([vol, vol - w])


# ----------------------------------------------------------------------------------------------
# Additive quadratic parameterisation
# ----------------------------------------------------------------------------------------------
def get_additive_quadratic_essentials(
    times_elasticities, times_satiations, goods_elasticities, goods_satiations, goods_cost=1
) -> dict[str, float]:
    """Aggregates S_theta, S_phi, B and H of the additive quadratic MTUEM."""
    theta = as_vector(times_elasticities)
    beta = as_vector(times_satiations)
    phi = as_vector(goods_elasticities)
    eta = as_vector(goods_satiations)
    P = as_vector(goods_cost)
    if P.size == 1:
        P = np.repeat(P, eta.size)
    return {
        "Sth": float(np.sum(theta / beta)),
        "Sph": float(np.sum(P * phi / eta)),
        "B": float(np.sum(1 / beta)),
        "H": float(np.sum(P**2 / eta)),
    }


def get_tw_additive_quadratic(work_elasticity, work_satiation, Sth, Sph, B, H, tau, Tc, Ec, w) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    num_1 = as_vector(work_elasticity)[0]
    num_2 = w * (Sph + np.asarray(Ec, dtype=float)) / H
    num_3 = (Sth - (tau - np.asarray(Tc, dtype=float))) / B
    den = as_vector(work_satiation)[0] + w**2 / H + 1 / B
    return _col((num_1 + num_2 - num_3) / den)


def get_ti_additive_quadratic(times_elasticities, times_satiations, Sth, B, tau, Tw, Tc) -> np.ndarray:
    theta = as_vector(times_elasticities)
    beta = as_vector(times_satiations)
    Tw = np.asarray(Tw, dtype=float).ravel()
    mu = (Sth - (tau - Tw - np.asarray(Tc, dtype=float))) / B
    return np.outer(np.ones(mu.size), theta / beta) - np.outer(mu, 1 / beta)


def get_xi_additive_quadratic(goods_elasticities, goods_satiations, goods_cost, Sph, H, Tw, Ec, w) -> np.ndarray:
    phi = as_vector(goods_elasticities)
    eta = as_vector(goods_satiations)
    P = as_vector(goods_cost)
    if P.size == 1:
        P = np.repeat(P, eta.size)
    Tw = np.asarray(Tw, dtype=float).ravel()
    lam = (Sph - (np.asarray(w, dtype=float) * Tw - np.asarray(Ec, dtype=float))) / H
    return np.outer(np.ones(lam.size), phi / eta) - np.outer(lam, P / eta)


def get_values_of_time_additive_quadratic(Sth, Sph, B, H, Tw, tau, Tc, Ec, w) -> np.ndarray:
    Tw = np.asarray(Tw, dtype=float).ravel()
    w = np.asarray(w, dtype=float)
    mu = (Sth - (tau - Tw - np.asarray(Tc, dtype=float))) / B
    lam = (Sph - (w * Tw - np.asarray(Ec, dtype=float))) / H
    vol = mu / lam
    return np.column_stack([vol, vol - w])


# ----------------------------------------------------------------------------------------------
# Error structure helpers
# ----------------------------------------------------------------------------------------------
def get_cond_err(mu: np.ndarray, rho: np.ndarray) -> dict[str, np.ndarray]:
    """Sequential conditional decomposition of standardised normal errors.

    Returns ``cond_mu`` (the standardised conditional residuals) and ``cond_sd`` (the
    conditional variances), exactly as ``get_cond_err`` in the R package.
    """
    mu = np.asarray(mu, dtype=float)
    neq = mu.shape[1]
    cmu = np.zeros_like(mu)
    csd = np.ones(neq)
    for j in range(1, neq):
        if j == 1:
            cmu[:, 1] = rho[1, 0] * mu[:, 0]
            csd[1] = 1 - rho[1, 0] ** 2
        elif j == 2:
            cmu[:, 2] = (
                (rho[1, 2] - rho[0, 2] * rho[0, 1]) * mu[:, 1] + (rho[0, 2] - rho[1, 2] * rho[0, 1]) * mu[:, 0]
            ) / csd[1]
            csd[2] = 1 - (rho[1, 2] ** 2 - 2 * rho[1, 0] * rho[1, 2] * rho[0, 2] + rho[0, 2] ** 2) / csd[1]
        else:
            idx = list(range(j - 1, -1, -1))
            inv_rho = np.linalg.pinv(rho[np.ix_(idx, idx)])
            cmu[:, j] = rho[j, idx] @ inv_rho @ mu[:, idx].T
            csd[j] = rho[j, j] - rho[j, idx] @ inv_rho @ rho[idx, j]
    with np.errstate(invalid="ignore", divide="ignore"):
        cmu = (mu - cmu) / np.sqrt(csd)
    return {"cond_mu": cmu, "cond_sd": csd}


def cholesky_to_sigma(cholesky, M: int) -> np.ndarray:
    """Covariance L L' from the row-wise lower triangle (L11, L21, L22, L31, ...)."""
    vals = as_vector(cholesky)
    if vals.size != M * (M + 1) // 2:
        raise ValueError(
            f"The length of cholesky ({vals.size}) must match M * (M + 1) / 2 where M is the "
            f"number of equations ({M})."
        )
    L = np.zeros((M, M))
    L[np.tril_indices(M)] = vals
    return L @ L.T


def cov_to_sig_rho(Sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sig = np.sqrt(np.diag(Sigma))
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = Sigma / np.outer(sig, sig)
    np.fill_diagonal(rho, 1.0)
    return sig, rho


def rho_from_upper(rho_vals, M: int) -> np.ndarray:
    """Correlation matrix from its upper triangle read row-wise (rho_12, rho_13, rho_23...)."""
    rho = np.eye(M)
    rho[np.triu_indices(M, 1)] = as_vector(rho_vals)
    return np.triu(rho) + np.triu(rho, 1).T


def complete_rows(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return x[np.all(np.isfinite(x), axis=1)]


def cov_complete(x: np.ndarray) -> np.ndarray:
    """``stats::cov(x, use = "complete.obs")``; all NaN when fewer than two rows are complete."""
    rows = complete_rows(x)
    if rows.shape[0] < 2:
        return np.full((rows.shape[1], rows.shape[1]), np.nan)
    return np.atleast_2d(np.cov(rows, rowvar=False, ddof=1))


def cor_complete(x: np.ndarray) -> np.ndarray:
    """``stats::cor(x, use = "complete.obs")``; all NaN when fewer than two rows are complete."""
    rows = complete_rows(x)
    if rows.shape[0] < 2:
        return np.full((rows.shape[1], rows.shape[1]), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.atleast_2d(np.corrcoef(rows, rowvar=False))


def near_pd(
    x: np.ndarray,
    corr: bool = False,
    keep_diag: bool = False,
    eig_tol: float = 1e-6,
    conv_tol: float = 1e-7,
    posd_tol: float = 1e-8,
    maxit: int = 100,
) -> np.ndarray:
    """Nearest positive definite matrix, a port of ``Matrix::nearPD`` (Higham 2002)."""
    x = np.asarray(x, dtype=float)
    x = (x + x.T) / 2
    n = x.shape[0]
    diag0 = np.diag(x).copy()
    D_S = np.zeros_like(x)
    X = x.copy()
    it = 0
    converged = False
    while it < maxit and not converged:
        Y = X
        R = Y - D_S
        d, Q = np.linalg.eigh(R)
        d, Q = d[::-1], Q[:, ::-1]
        p = d > eig_tol * d[0]
        if not np.any(p):
            raise np.linalg.LinAlgError("Matrix seems negative semi-definite")
        Q = Q[:, p]
        X = (Q * d[p]) @ Q.T
        D_S = X - R
        if corr:
            np.fill_diagonal(X, 1.0)
        elif keep_diag:
            np.fill_diagonal(X, diag0)
        conv = np.abs(Y - X).sum(axis=1).max() / np.abs(Y).sum(axis=1).max()
        it += 1
        converged = conv <= conv_tol
    d, Q = np.linalg.eigh(X)
    d, Q = d[::-1], Q[:, ::-1]
    eps = posd_tol * abs(d[0])
    if d[n - 1] < eps:
        d = np.where(d < eps, eps, d)
        o_diag = np.diag(X).copy()
        X = Q @ (d[:, None] * Q.T)
        D = np.sqrt(np.maximum(eps, o_diag) / np.diag(X))
        X = D[:, None] * X * D[None, :]
    if corr:
        np.fill_diagonal(X, 1.0)
    elif keep_diag:
        np.fill_diagonal(X, diag0)
    return X


def mvn_logpdf_cond(err: np.ndarray, sig: np.ndarray, rho: np.ndarray) -> np.ndarray | None:
    """Zero-mean multivariate normal log density through ``get_cond_err``.

    Returns None when the covariance is degenerate, signalling a zero likelihood.
    """
    if np.any(~np.isfinite(sig)) or np.any(sig <= 1e-6) or np.any(~np.isfinite(rho)):
        return None
    eigs = np.linalg.eigvalsh(rho)
    if np.any(~np.isfinite(eigs)) or np.any(eigs <= 1e-6):
        return None
    mu = err / sig
    cond = get_cond_err(mu, rho)
    cmu, csd = cond["cond_mu"], cond["cond_sd"]
    const = np.sum(np.log(sig)) + 0.5 * np.sum(np.log(csd)) + 0.5 * cmu.shape[1] * np.log(2 * np.pi)
    return -0.5 * np.sum(cmu**2, axis=1) - const
