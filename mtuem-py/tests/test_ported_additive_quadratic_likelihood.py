"""Port of tests/testthat/test_additive_quadratic_likelihood.R."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose
from scipy import stats

import mtuem
from helper_mtuem import inputs_of, mtuem_chol_vec, mtuem_dmvnorm_ref


def aq_fixture(n=60, seed=1, sd=1, goods_cost=None):
    """Additive quadratic fixture: observed allocations are the optimal ones plus noise."""
    goods_cost = {"P1": 1, "P2": 1} if goods_cost is None else goods_cost
    rng = np.random.default_rng(seed)
    tau = 168
    Tc = rng.uniform(40, 80, n)
    Ec = rng.uniform(50, 150, n)
    w = rng.uniform(10, 25, n)

    times_elasticities = {"theta_1": 800, "theta_2": 3000}
    times_satiations = {"beta_1": 20, "beta_2": 300}
    goods_elasticities = {"phi_1": 100, "phi_2": 80}
    goods_satiations = {"eta_1": 1, "eta_2": 0.5}
    work_elasticity, work_satiation = 0, 10

    ess = mtuem.utils.get_additive_quadratic_essentials(
        times_elasticities, times_satiations, goods_elasticities, goods_satiations, goods_cost
    )
    tw = mtuem.utils.get_tw_additive_quadratic(
        work_elasticity, work_satiation, ess["Sth"], ess["Sph"], ess["B"], ess["H"], tau, Tc, Ec, w
    )
    ti = mtuem.utils.get_ti_additive_quadratic(times_elasticities, times_satiations, ess["Sth"], ess["B"], tau, tw, Tc)
    xj = mtuem.utils.get_xi_additive_quadratic(
        goods_elasticities, goods_satiations, goods_cost, ess["Sph"], ess["H"], tw, Ec, w
    )
    nms = ["Tw", "Tf1", "Tf2", "Ef1", "Ef2"]
    opt = np.column_stack([tw, ti, xj])
    obs = opt + rng.normal(0, sd, opt.shape)
    database = pd.concat(
        [pd.DataFrame({"PeID": np.arange(1, n + 1), "Tc": Tc, "EcI": Ec, "w": w}), pd.DataFrame(obs, columns=nms)],
        axis=1,
    )
    return SimpleNamespace(
        inputs=inputs_of(database), database=database, opt=opt, obs=obs, ess=ess, n=n, tau=tau,
        Tc=Tc, Ec=Ec, w=w, cols={c: i for i, c in enumerate(nms)},
        mtuem_settings={
            "work_times": "Tw", "free_times": ["Tf1", "Tf2"], "free_goods": ["Ef1", "Ef2"],
            "goods_cost": goods_cost, "work_elasticity": work_elasticity, "work_satiation": work_satiation,
            "times_elasticities": times_elasticities, "times_satiations": times_satiations,
            "goods_elasticities": goods_elasticities, "goods_satiations": goods_satiations,
            "Tc": Tc, "Ec": Ec, "w": w, "tau": tau,
        },
    )


def test_additive_quadratic_likelihood_reproduces_the_multivariate_normal_density():
    # With all five equations the residual covariance is singular (the budgets make
    # them collinear), so a Cholesky factor is supplied for a reduced set instead.
    fx = aq_fixture(n=50, seed=61)
    Sigma = np.array([[1.5, 0.3, -0.2], [0.3, 1.0, 0.1], [-0.2, 0.1, 2.0]])
    settings = dict(fx.mtuem_settings, omitted_times="Tf2", omitted_goods="Ef2", cholesky=mtuem_chol_vec(Sigma))

    L = mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)

    idx = [fx.cols[c] for c in ("Tw", "Tf1", "Ef1")]
    assert_allclose(L, mtuem_dmvnorm_ref((fx.obs - fx.opt)[:, idx], Sigma), rtol=1e-10)


def test_omitting_all_but_one_equation_falls_back_to_the_univariate_normal():
    fx = aq_fixture(n=30, seed=62)
    settings = dict(fx.mtuem_settings, omitted_times=["Tf1", "Tf2"], omitted_goods=["Ef1", "Ef2"], sig=[1.3])

    L = mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)

    err = fx.obs[:, 0] - fx.opt[:, 0]
    assert_allclose(L, stats.norm.pdf(err, scale=1.3), rtol=1e-12)


def test_omitting_an_equation_keeps_its_parameters_in_the_aggregates():
    # Dropping Tf2 from the likelihood must not change the predicted Tw, because
    # theta_2 and beta_2 still enter S_theta and B.
    fx = aq_fixture(n=30, seed=63)
    settings = dict(fx.mtuem_settings, omitted_times="Tf2", omitted_goods="Ef2")

    pred_full = mtuem.likelihood.additive_quadratic(fx.mtuem_settings, "prediction", fx.inputs)
    pred_omit = mtuem.likelihood.additive_quadratic(settings, "prediction", fx.inputs)

    pd.testing.assert_frame_equal(pred_omit, pred_full, rtol=1e-12)


def test_omitted_categories_must_name_existing_free_categories():
    fx = aq_fixture(n=10, seed=64)
    settings = dict(fx.mtuem_settings, omitted_times="Tf9")

    with pytest.raises(ValueError, match="not in free_times / free_goods"):
        mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)


def test_omitting_every_equation_is_rejected():
    fx = aq_fixture(n=10, seed=65)
    settings = dict(fx.mtuem_settings, work_times=[], omitted_times=["Tf1", "Tf2"], omitted_goods=["Ef1", "Ef2"])

    with pytest.raises(ValueError, match="nothing left to estimate"):
        mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)


def test_the_cholesky_length_is_checked_against_the_surviving_equations():
    # Three equations are omitted, so three Cholesky entries are needed for the two
    # survivors (Tw, Ef1), not fifteen for all five.
    fx = aq_fixture(n=10, seed=66)
    settings = dict(fx.mtuem_settings, omitted_times=["Tf1", "Tf2"], omitted_goods="Ef2", cholesky=[1.0] * 15)

    with pytest.raises(ValueError, match="must match M"):
        mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)


def test_additive_quadratic_prediction_closes_both_budgets_and_reports_values_of_time():
    fx = aq_fixture(n=25, seed=67, goods_cost={"P1": 1, "P2": 2})

    pred = mtuem.likelihood.additive_quadratic(fx.mtuem_settings, "prediction", fx.inputs)

    # No residual time or goods category is reported, because the parameterised
    # categories already exhaust both budgets.
    assert list(pred.columns) == ["Tw", "Tf1", "Tf2", "Ef1", "Ef2", "VoL", "VTAW"]
    assert_allclose(pred["Tw"] + pred["Tf1"] + pred["Tf2"], fx.tau - fx.Tc, rtol=1e-9)
    assert_allclose(1 * pred["Ef1"] + 2 * pred["Ef2"], fx.w * pred["Tw"] - fx.Ec, rtol=1e-9)
    assert_allclose(pred["VoL"], fx.w + pred["VTAW"], rtol=1e-9)


def test_additive_quadratic_handles_the_engine_bookkeeping_functionalities():
    fx = aq_fixture(n=11, seed=68)
    named = dict(fx.mtuem_settings, componentName="aq_component")
    inp = fx.inputs
    f = mtuem.likelihood.additive_quadratic

    assert_allclose(f(fx.mtuem_settings, "validate", inp), np.ones(11))
    assert np.all(np.isnan(f(fx.mtuem_settings, "zero_LL", inp)))
    assert f(named, "preprocess", inp)["componentName"] == "aq_component"
    with pytest.raises(ValueError, match="Invalid value of argument"):
        f(fx.mtuem_settings, "nonsense", inp)


def test_additive_quadratic_output_matches_estimate_and_report_prints_the_values_of_time(capsys):
    fx = aq_fixture(n=20, seed=69)
    settings = dict(fx.mtuem_settings, omitted_times="Tf2", omitted_goods="Ef2", sig=[1, 1, 1])
    f = mtuem.likelihood.additive_quadratic

    assert_allclose(f(settings, "output", fx.inputs), f(settings, "estimate", fx.inputs))
    f(settings, "report", fx.inputs)
    assert "VoL" in capsys.readouterr().out


def test_optimal_tw_false_builds_the_other_equations_from_observed_work_time():
    fx = aq_fixture(n=25, seed=70)
    settings = dict(
        fx.mtuem_settings, omitted_times="Tf2", omitted_goods=["Ef1", "Ef2"], optimal_tw=False, sig=[1, 1]
    )

    L = mtuem.likelihood.additive_quadratic(settings, "estimate", fx.inputs)

    ti_from_obs = mtuem.utils.get_ti_additive_quadratic(
        settings["times_elasticities"], settings["times_satiations"], fx.ess["Sth"], fx.ess["B"], fx.tau,
        fx.database["Tw"].to_numpy(), fx.Tc,
    )
    err = np.column_stack([fx.obs[:, 0] - fx.opt[:, 0], fx.obs[:, 1] - ti_from_obs[:, 0]])
    assert_allclose(L, mtuem_dmvnorm_ref(err, np.eye(2)), rtol=1e-10)
