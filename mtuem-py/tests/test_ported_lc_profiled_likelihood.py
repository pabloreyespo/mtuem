"""Port of tests/testthat/test_lc_profiled_likelihood.R."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import mtuem
from helper_mtuem import cov2cor, inputs_of, mtuem_dmvnorm_ref


def lc_fixture(n=400, seed=1):
    """Two-class fixture with a work equation plus one free time equation."""
    rng = np.random.default_rng(seed)
    tau = 168
    Tc = rng.uniform(40, 90, n)
    Ec = rng.uniform(50, 200, n)
    w = rng.uniform(8, 25, n)

    we1 = {"Theta": 1, "Phi": 1, "thw": 0}
    we2 = {"Theta": 1, "Phi": 2.5, "thw": 0.5}
    te1, te2 = {"th1": 0.4}, {"th1": 0.15}

    tw1 = mtuem.get_tw_thph(we1, tau, Tc, Ec, w)
    tw2 = mtuem.get_tw_thph(we2, tau, Tc, Ec, w)
    opt1 = np.column_stack([tw1, mtuem.get_ti_thph(te1, we1["Theta"], tw1, tau, Tc)])
    opt2 = np.column_stack([tw2, mtuem.get_ti_thph(te2, we2["Theta"], tw2, tau, Tc)])

    S1 = np.array([[4.0, 1.0], [1.0, 2.25]])
    S2 = np.array([[9.0, -1.5], [-1.5, 3.24]])
    pi1 = 0.6
    z = rng.binomial(1, 1 - pi1, n) + 1

    e1 = rng.normal(size=(n, 2)) @ np.linalg.cholesky(S1).T
    e2 = rng.normal(size=(n, 2)) @ np.linalg.cholesky(S2).T
    obs = np.where((z == 1)[:, None], opt1 + e1, opt2 + e2)

    database = pd.DataFrame({"PeID": np.arange(1, n + 1), "Tw": obs[:, 0], "Tf1": obs[:, 1], "Tc": Tc, "EcI": Ec, "w": w})
    class_settings = [
        {"work_times": "Tw", "free_times": "Tf1", "work_elasticities": we1, "times_elasticities": te1,
         "Tc": Tc, "Ec": Ec, "w": w, "tau": tau},
        {"work_times": "Tw", "free_times": "Tf1", "work_elasticities": we2, "times_elasticities": te2,
         "Tc": Tc, "Ec": Ec, "w": w, "tau": tau},
    ]
    return SimpleNamespace(
        apollo_inputs=inputs_of(database), database=database, n=n,
        lc_settings={"class_settings": class_settings, "classProb": [np.full(n, pi1), np.full(n, 1 - pi1)]},
        class_settings=class_settings, pi1=pi1, Sigma_true=[S1, S2],
    )


def test_the_profiled_likelihood_handles_the_apollo_bookkeeping_functionalities():
    fx = lc_fixture(n=60, seed=71)
    named = dict(fx.lc_settings, componentName="lc_component")
    f = mtuem.mtuem_lc_profiled_likelihood
    inp = fx.apollo_inputs

    assert_allclose(f(fx.lc_settings, "validate", inp), np.ones(60))
    assert np.all(np.isnan(f(fx.lc_settings, "zero_LL", inp)))
    assert f(named, "preprocess", inp)["componentName"] == "lc_component"
    assert f(fx.lc_settings, "preprocess", inp)["componentName"] == ".."
    with pytest.raises(ValueError, match="Invalid value of argument"):
        f(fx.lc_settings, "nonsense", inp)


def test_the_profiled_likelihood_is_a_proper_mixture_of_the_class_densities():
    fx = lc_fixture(n=300, seed=72)

    P = mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    prof = mtuem.mtuem_get_profiled_sigma()

    E1 = mtuem.mtuem_residuals(fx.class_settings[0], fx.apollo_inputs)
    E2 = mtuem.mtuem_residuals(fx.class_settings[1], fx.apollo_inputs)
    P_ref = fx.pi1 * mtuem_dmvnorm_ref(E1, prof["Sigma"][0]) + (1 - fx.pi1) * mtuem_dmvnorm_ref(E2, prof["Sigma"][1])

    assert_allclose(P, P_ref, rtol=1e-8)
    assert np.all(P > 0)


def test_the_profiled_likelihood_is_path_independent_across_repeated_calls():
    # The inner fixed point must restart from the deterministic unweighted covariance
    # every time, otherwise numerical gradients would be evaluated on a moving objective.
    fx = lc_fixture(n=200, seed=73)

    P_first = mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    P_third = mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)

    assert_array_equal(P_first, P_third)


def test_the_fixed_point_converges_on_a_stationary_covariance():
    fx = lc_fixture(n=300, seed=74)

    mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    prof = mtuem.mtuem_get_profiled_sigma()

    # One more weighted second-moment update leaves Sigma unchanged.
    E = [mtuem.mtuem_residuals(cs, fx.apollo_inputs) for cs in fx.class_settings]
    W = prof["posterior"]
    for s in range(2):
        ws = W[:, s]
        Es = E[s] * np.sqrt(ws)[:, None]
        Sigma_next = Es.T @ Es / ws.sum()
        assert_allclose(Sigma_next, prof["Sigma"][s], rtol=1e-7)


def test_mtuem_get_profiled_sigma_reports_sigmas_correlations_posteriors_and_shares():
    fx = lc_fixture(n=300, seed=75)

    mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    prof = mtuem.mtuem_get_profiled_sigma()

    for s in range(2):
        assert_allclose(prof["sig"][s], np.sqrt(np.diag(prof["Sigma"][s])), rtol=1e-12)
        assert_allclose(prof["rho"][s], cov2cor(prof["Sigma"][s]), rtol=1e-10)
        assert_allclose(np.diag(prof["rho"][s]), np.ones(2), rtol=1e-12)
    assert prof["posterior"].shape == (300, 2)
    assert_allclose(prof["posterior"].sum(axis=1), np.ones(300), rtol=1e-12)
    assert_allclose(prof["shares"], [fx.pi1, 1 - fx.pi1], rtol=1e-12)


def test_a_single_class_specification_reduces_to_the_plain_multivariate_normal():
    fx = lc_fixture(n=200, seed=76)
    one_class = {"class_settings": fx.class_settings[:1], "classProb": [np.ones(200)]}

    P = mtuem.mtuem_lc_profiled_likelihood(one_class, "estimate", fx.apollo_inputs)

    # With a single class the posterior weights are all one, so the profiled covariance
    # is the (uncentred) residual second moment.
    E = mtuem.mtuem_residuals(fx.class_settings[0], fx.apollo_inputs)
    Sigma = E.T @ E / E.shape[0]
    assert_allclose(P, mtuem_dmvnorm_ref(E, Sigma), rtol=1e-8)


def test_prediction_returns_per_class_allocations_and_posterior_weights():
    fx = lc_fixture(n=100, seed=77)

    pred = mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "prediction", fx.apollo_inputs)

    assert len(pred) == 100
    assert {"post_class1", "post_class2"} <= set(pred.columns)
    assert all(c.endswith("_class1") for c in pred.columns[:6])
    assert_allclose(pred[["post_class1", "post_class2"]].sum(axis=1), np.ones(100), rtol=1e-12)


def test_output_matches_estimate_and_report_prints_per_class_values_of_time(capsys):
    fx = lc_fixture(n=80, seed=78)
    f = mtuem.mtuem_lc_profiled_likelihood

    assert_allclose(f(fx.lc_settings, "output", fx.apollo_inputs), f(fx.lc_settings, "estimate", fx.apollo_inputs))
    capsys.readouterr()
    f(fx.lc_settings, "report", fx.apollo_inputs)
    assert "Class 2 VoL" in capsys.readouterr().out


def test_em_settings_control_the_fixed_point_and_change_the_answer_when_truncated():
    fx = lc_fixture(n=200, seed=79)
    truncated = dict(fx.lc_settings, em_settings={"maxIterations": 1, "tol": 0})

    P_conv = mtuem.mtuem_lc_profiled_likelihood(fx.lc_settings, "estimate", fx.apollo_inputs)
    P_trunc = mtuem.mtuem_lc_profiled_likelihood(truncated, "estimate", fx.apollo_inputs)

    assert np.all(np.isfinite(P_trunc))
    assert not np.allclose(P_conv, P_trunc, rtol=1.5e-8, atol=0)


def _ic_model(n_est, fixed, maximum, nObs):
    return SimpleNamespace(
        estimate=pd.Series(np.zeros(n_est), index=[f"b{i}" for i in range(1, n_est + 1)]),
        apollo_fixed=fixed, maximum=maximum, nObs=nObs,
    )


def test_mtuem_profile_ic_counts_the_profiled_covariance_parameters():
    # 10 estimated entries of which 3 are fixed, 2 classes, 2 equations per class: k = 7 + 2 * 3 = 13.
    model = _ic_model(10, ["b1", "b2", "b3"], -1234.5, 500)

    ic = mtuem.mtuem_profile_ic(model, nClass=2, nEq=2)

    assert ic["k"] == 13
    assert_allclose(ic["AIC"], -2 * -1234.5 + 2 * 13)
    assert_allclose(ic["BIC"], -2 * -1234.5 + 13 * np.log(500))


def test_mtuem_profile_ic_penalises_a_larger_covariance_more_heavily():
    model = _ic_model(6, [], -100, 200)

    small = mtuem.mtuem_profile_ic(model, nClass=2, nEq=2)
    large = mtuem.mtuem_profile_ic(model, nClass=2, nEq=3)

    assert large["k"] - small["k"] == 2 * (6 - 3)
    assert large["BIC"] > small["BIC"]
