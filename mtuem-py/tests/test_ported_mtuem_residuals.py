"""Port of tests/testthat/test_mtuem_residuals.R."""

import numpy as np
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import inputs_of, mtuem_fixture_thph


def test_mtuem_residuals_returns_observed_minus_optimal_for_every_equation():
    fx = mtuem_fixture_thph(n=25, seed=12)

    err = mtuem.likelihood.residuals(fx.mtuem_settings, fx.apollo_inputs)

    assert err.shape == fx.opt.shape
    # Columns are Tw, Tf1, Tf2 by position (numpy arrays carry no names).
    assert fx.names == ["Tw", "Tf1", "Tf2"]
    assert_allclose(err, fx.obs - fx.opt, atol=1e-10)


def test_mtuem_residuals_fills_in_optional_settings_with_their_defaults():
    fx = mtuem_fixture_thph(n=15, times_elasticities={}, seed=3)
    minimal = {
        "work_times": "Tw",
        "work_elasticities": fx.mtuem_settings["work_elasticities"],
        "Tc": fx.Tc, "Ec": fx.Ec, "w": fx.w, "tau": fx.tau,
    }

    err = mtuem.likelihood.residuals(minimal, fx.apollo_inputs)

    assert err.shape == (15, 1)
    assert_allclose(err.ravel(), fx.obs[:, 0] - fx.opt[:, 0], atol=1e-10)


def test_mtuem_residuals_covers_free_goods_as_well_as_free_times():
    fx = mtuem_fixture_thph(n=30, times_elasticities={"th1": 0.3}, goods_elasticities={"ph1": 0.4}, seed=6)

    err = mtuem.likelihood.residuals(fx.mtuem_settings, fx.apollo_inputs)

    assert fx.names == ["Tw", "Tf1", "Ef1"]
    assert_allclose(err, fx.obs - fx.opt, atol=1e-10)


def test_mtuem_residuals_uses_obs_matrix_in_preference_to_the_database():
    fx = mtuem_fixture_thph(n=20, seed=8)
    inputs_with_matrix = inputs_of(fx.database, obs_matrix=fx.obs + 5)

    err_db = mtuem.likelihood.residuals(fx.mtuem_settings, fx.apollo_inputs)
    err_mat = mtuem.likelihood.residuals(fx.mtuem_settings, inputs_with_matrix)

    assert_allclose(err_mat - err_db, np.full((20, 3), 5.0), atol=1e-10)


def test_optimal_tw_false_derives_the_other_equations_from_observed_work_time():
    fx = mtuem_fixture_thph(n=20, times_elasticities={"th1": 0.3}, seed=15)
    settings = dict(fx.mtuem_settings, optimal_tw=False)

    err = mtuem.likelihood.residuals(settings, fx.apollo_inputs)

    # The free-time equation is built from the observed Tw, while the work
    # equation still compares against the optimal Tw.
    ti_from_obs = mtuem.get_ti_thph(
        settings["times_elasticities"], settings["work_elasticities"]["Theta"],
        fx.database["Tw"].to_numpy(), fx.tau, fx.Tc,
    )
    assert_allclose(err[:, 1], fx.database["Tf1"].to_numpy() - ti_from_obs[:, 0], atol=1e-10)
    assert_allclose(err[:, 0], fx.obs[:, 0] - fx.opt[:, 0], atol=1e-10)


def test_mtuem_residuals_is_exactly_zero_when_observations_equal_the_optimum():
    fx = mtuem_fixture_thph(n=12, seed=2, sd=0)

    err = mtuem.likelihood.residuals(fx.mtuem_settings, fx.apollo_inputs)

    assert_allclose(err, np.zeros((12, 3)), atol=1e-10)
