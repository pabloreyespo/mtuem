"""Shared fixtures for the tests ported from the R testthat suite (tests/testthat/helper-mtuem.R)."""

from types import SimpleNamespace

import numpy as np
import pandas as pd

import mtuem
from mtuem import _engine as engine
from mtuem.config import Control


def inputs_of(database: pd.DataFrame, **extra) -> SimpleNamespace:
    """Minimal engine inputs: any object with a ``database`` (and optionally ``obs_matrix``)."""
    return SimpleNamespace(database=database, **extra)


def mtuem_fixture_thph(
    n=200,
    work_elasticities=None,
    times_elasticities=None,
    goods_elasticities=None,
    tau=168,
    seed=1,
    sd=2,
):
    """Theta-Phi fixture: observed allocations are the optimal ones plus normal noise."""
    work_elasticities = {"Theta": 1, "Phi": 1, "thw": 0} if work_elasticities is None else work_elasticities
    times_elasticities = {"th1": 0.3, "th2": 0.2} if times_elasticities is None else times_elasticities
    goods_elasticities = {} if goods_elasticities is None else goods_elasticities
    rng = np.random.default_rng(seed)
    Tc = rng.uniform(40, 80, n)
    Ec = rng.uniform(50, 150, n)
    w = rng.uniform(10, 25, n)

    tw_opt = mtuem.utils.get_tw_thph(work_elasticities, tau, Tc, Ec, w)
    opt = tw_opt
    nms = ["Tw"]
    if len(times_elasticities) > 0:
        ti_opt = mtuem.utils.get_ti_thph(times_elasticities, work_elasticities["Theta"], tw_opt, tau, Tc)
        opt = np.column_stack([opt, ti_opt])
        nms += [f"Tf{i + 1}" for i in range(len(times_elasticities))]
    if len(goods_elasticities) > 0:
        xj_opt = mtuem.utils.get_xi_thph(goods_elasticities, 1, work_elasticities["Phi"], tw_opt, Ec, w)
        opt = np.column_stack([opt, xj_opt])
        nms += [f"Ef{i + 1}" for i in range(len(goods_elasticities))]
    M = opt.shape[1]
    obs = opt + rng.normal(0, sd, (n, M)) if sd > 0 else opt.copy()

    database = pd.DataFrame({"PeID": np.arange(1, n + 1), "Tc": Tc, "EcI": Ec, "w": w})
    database = pd.concat([database, pd.DataFrame(obs, columns=nms)], axis=1)

    nt, ng = len(times_elasticities), len(goods_elasticities)
    settings = {
        "work_times": "Tw",
        "free_times": nms[1 : 1 + nt] if nt else [],
        "free_goods": nms[len(nms) - ng :] if ng else [],
        "work_elasticities": work_elasticities,
        "times_elasticities": times_elasticities,
        "goods_elasticities": goods_elasticities,
        "Tc": Tc, "Ec": Ec, "w": w, "tau": tau,
    }
    return SimpleNamespace(
        database=database, inputs=inputs_of(database), mtuem_settings=settings,
        opt=opt, obs=obs, names=nms, Tc=Tc, Ec=Ec, w=w, tau=tau, n=n,
    )


def mtuem_dmvnorm_ref(E, Sigma) -> np.ndarray:
    """Zero-mean multivariate normal density through inv() and slogdet()."""
    E = np.asarray(E, dtype=float)
    if E.ndim == 1:
        E = E[:, None]
    M = E.shape[1]
    Sinv = np.linalg.inv(Sigma)
    logdet = np.linalg.slogdet(Sigma)[1]
    q = np.sum((E @ Sinv) * E, axis=1)
    return np.exp(-0.5 * (q + logdet + M * np.log(2 * np.pi)))


def mtuem_chol_vec(Sigma) -> list[float]:
    """Lower Cholesky factor flattened row-major (L11, L21, L22, L31, ...)."""
    L = np.linalg.cholesky(np.asarray(Sigma, dtype=float))
    return [float(L[r, c]) for r in range(L.shape[0]) for c in range(r + 1)]


def cov2cor(S) -> np.ndarray:
    S = np.asarray(S, dtype=float)
    d = np.sqrt(np.diag(S))
    return S / np.outer(d, d)


def mtuem_fake_model(estimate: dict, varcov) -> engine.Estimate:
    """Minimal estimate object that engine.delta_method accepts."""
    est = pd.Series(estimate, dtype=float)
    vc = pd.DataFrame(np.asarray(varcov, dtype=float), index=est.index, columns=est.index)
    return engine.Estimate(
        estimate=est, beta_start=est.copy(), fixed=[], control=Control(),
        maximum=0.0, LLStart=0.0, varcov=vc, robvarcov=vc.copy(),
    )
