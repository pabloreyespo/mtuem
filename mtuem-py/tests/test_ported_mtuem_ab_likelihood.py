"""Port of tests/testthat/test_mtuem_ab_likelihood.R."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy import stats

import mtuem
from helper_mtuem import inputs_of, mtuem_chol_vec, mtuem_dmvnorm_ref


def ab_fixture(n=60, times_elasticities=None, goods_elasticities=None, seed=1, sd=2):
    """Alpha-beta fixture: observed allocations are the optimal ones plus noise."""
    times_elasticities = {"th1": 0.1} if times_elasticities is None else times_elasticities
    goods_elasticities = {} if goods_elasticities is None else goods_elasticities
    rng = np.random.default_rng(seed)
    tau = 168
    Tc = rng.uniform(40, 80, n)
    Ec = rng.uniform(50, 150, n)
    w = rng.uniform(10, 25, n)
    we = {"alpha": 0.3, "beta": 0.25}

    tw = mtuem.utils.get_tw_albe(we, tau, Tc, Ec, w)
    opt, nms = tw, ["Tw"]
    if times_elasticities:
        opt = np.column_stack([opt, mtuem.utils.get_ti_albe(times_elasticities, we["beta"], tw, tau, Tc)])
        nms += [f"Tf{i + 1}" for i in range(len(times_elasticities))]
    if goods_elasticities:
        opt = np.column_stack([opt, mtuem.utils.get_xi_albe(goods_elasticities, 1, we["alpha"], tw, Ec, w)])
        nms += [f"Ef{i + 1}" for i in range(len(goods_elasticities))]
    obs = opt + rng.normal(0, sd, opt.shape)
    database = pd.concat(
        [pd.DataFrame({"PeID": np.arange(1, n + 1), "Tc": Tc, "EcI": Ec, "w": w}), pd.DataFrame(obs, columns=nms)],
        axis=1,
    )
    nt, ng = len(times_elasticities), len(goods_elasticities)
    return SimpleNamespace(
        inputs=inputs_of(database), database=database, opt=opt, obs=obs, n=n, tau=tau, Tc=Tc, Ec=Ec, w=w,
        mtuem_settings={
            "work_times": "Tw",
            "free_times": nms[1 : 1 + nt] if nt else [],
            "free_goods": nms[len(nms) - ng :] if ng else [],
            "work_elasticities": we,
            "times_elasticities": times_elasticities,
            "goods_elasticities": goods_elasticities,
            "Tc": Tc, "Ec": Ec, "w": w, "tau": tau,
        },
    )


def test_alpha_beta_single_equation_likelihood_equals_the_univariate_normal_density():
    fx = ab_fixture(n=40, times_elasticities={}, seed=51)
    settings = dict(fx.mtuem_settings, sig=[2.4])

    L = mtuem.likelihood.alpha_beta(settings, "estimate", fx.inputs)

    err = fx.obs[:, 0] - fx.opt[:, 0]
    assert_allclose(L, stats.norm.pdf(err, scale=2.4), rtol=1e-12)


def test_alpha_beta_cholesky_specification_reproduces_the_multivariate_normal_density():
    fx = ab_fixture(n=50, times_elasticities={"th1": 0.1}, seed=52)
    Sigma = np.array([[4.0, 1.1], [1.1, 2.25]])
    settings = dict(fx.mtuem_settings, cholesky=mtuem_chol_vec(Sigma))

    L = mtuem.likelihood.alpha_beta(settings, "estimate", fx.inputs)

    assert_allclose(L, mtuem_dmvnorm_ref(fx.obs - fx.opt, Sigma), rtol=1e-10)


def test_alpha_beta_rejects_a_cholesky_vector_of_the_wrong_length():
    fx = ab_fixture(n=10, seed=53)
    settings = dict(fx.mtuem_settings, cholesky=[1, 0.2])  # 2 entries for a 2-equation model

    with pytest.raises(ValueError, match="must match M"):
        mtuem.likelihood.alpha_beta(settings, "estimate", fx.inputs)


def test_alpha_beta_guards_against_a_degenerate_sigma():
    fx = ab_fixture(n=14, times_elasticities={}, seed=54)
    settings = dict(fx.mtuem_settings, sig=[0])

    assert_allclose(mtuem.likelihood.alpha_beta(settings, "estimate", fx.inputs), np.zeros(14))


def test_alpha_beta_get_covar_returns_the_residual_covariance_correlation_and_sigmas():
    fx = ab_fixture(n=80, times_elasticities={"th1": 0.1}, seed=55)

    cv = mtuem.likelihood.alpha_beta(fx.mtuem_settings, "get_covar", fx.inputs)

    err = fx.obs - fx.opt
    assert_allclose(cv["covar"].to_numpy(), np.cov(err, rowvar=False), rtol=1e-10)
    assert_allclose(cv["corr"].to_numpy(), np.corrcoef(err, rowvar=False), rtol=1e-10)
    assert_allclose(cv["sigma"].to_numpy(), np.sqrt(np.diag(np.cov(err, rowvar=False))), rtol=1e-10)


def test_alpha_beta_handles_the_engine_bookkeeping_functionalities():
    fx = ab_fixture(n=13, seed=56)
    named = dict(fx.mtuem_settings, componentName="ab_component")
    inp = fx.inputs

    assert_allclose(mtuem.likelihood.alpha_beta(fx.mtuem_settings, "validate", inp), np.ones(13))
    assert np.all(np.isnan(mtuem.likelihood.alpha_beta(fx.mtuem_settings, "zero_LL", inp)))
    assert mtuem.likelihood.alpha_beta(named, "preprocess", inp)["componentName"] == "ab_component"
    assert mtuem.likelihood.alpha_beta(fx.mtuem_settings, "gradient", inp) == {"like": None, "grad": None}
    with pytest.raises(ValueError, match="Invalid value of argument"):
        mtuem.likelihood.alpha_beta(fx.mtuem_settings, "nonsense", inp)


def test_alpha_beta_prediction_closes_both_budgets_and_reports_values_of_time():
    fx = ab_fixture(n=25, times_elasticities={"th1": 0.1}, goods_elasticities={"ph1": 0.2}, seed=57)

    pred = mtuem.likelihood.alpha_beta(fx.mtuem_settings, "prediction", fx.inputs)

    assert list(pred.columns) == ["Tw", "Tf1", "Tfi", "Ef1", "Xfj", "VoL", "VTAW"]
    assert_allclose(pred["Tw"] + pred["Tf1"] + pred["Tfi"], fx.tau - fx.Tc, rtol=1e-9)
    assert_allclose(pred["Ef1"] + pred["Xfj"], fx.w * pred["Tw"] - fx.Ec, rtol=1e-9)
    assert_allclose(pred["VoL"], fx.w + pred["VTAW"], rtol=1e-9)


def test_alpha_beta_output_matches_estimate_and_report_prints_the_values_of_time(capsys):
    fx = ab_fixture(n=20, seed=58)
    settings = dict(fx.mtuem_settings, sig=[2, 2])

    out = mtuem.likelihood.alpha_beta(settings, "output", fx.inputs)
    est = mtuem.likelihood.alpha_beta(settings, "estimate", fx.inputs)
    assert_allclose(out, est)

    mtuem.likelihood.alpha_beta(settings, "report", fx.inputs)
    assert "VTAW" in capsys.readouterr().out
