"""Port of tests/testthat/test_additive_quadratic.R."""

import numpy as np
from numpy.testing import assert_allclose

import mtuem


def aq_params():
    return {
        "times_elasticities": {"theta_1": 800, "theta_i": 3000},
        "times_satiations": {"beta_1": 20, "beta_i": 300},
        "goods_elasticities": {"phi_1": 100, "phi_i": 80},
        "goods_satiations": {"eta_1": 1, "eta_i": 0.5},
        "goods_cost": {"P1": 1, "Pi": 2},
    }


def essentials(p, goods_cost=None):
    return mtuem.utils.get_additive_quadratic_essentials(
        p["times_elasticities"], p["times_satiations"], p["goods_elasticities"], p["goods_satiations"],
        p["goods_cost"] if goods_cost is None else goods_cost,
    )


def vals(d):
    return np.array(list(d.values()), dtype=float)


def test_get_additive_quadratic_essentials_computes_the_four_sufficient_statistics():
    p = aq_params()
    theta, beta = vals(p["times_elasticities"]), vals(p["times_satiations"])
    phi, eta, P = vals(p["goods_elasticities"]), vals(p["goods_satiations"]), vals(p["goods_cost"])

    ess = essentials(p)

    assert_allclose(ess["Sth"], np.sum(theta / beta), rtol=1e-12)
    assert_allclose(ess["B"], np.sum(1 / beta), rtol=1e-12)
    assert_allclose(ess["Sph"], np.sum(P * phi / eta), rtol=1e-12)
    assert_allclose(ess["H"], np.sum(P**2 / eta), rtol=1e-12)


def test_get_additive_quadratic_essentials_recycles_a_scalar_price_across_all_goods():
    p = aq_params()

    scalar = essentials(p, goods_cost=1)
    vector = essentials(p, goods_cost=[1, 1])

    assert scalar.keys() == vector.keys()
    for k in scalar:
        assert_allclose(scalar[k], vector[k])


def test_get_tw_additive_quadratic_matches_the_analytical_work_time_expression():
    p = aq_params()
    tau, Tc, Ec, w = 168, 80, 200, 15
    theta_w, beta_w = 0, 10
    ess = essentials(p)

    tw = mtuem.utils.get_tw_additive_quadratic(theta_w, beta_w, ess["Sth"], ess["Sph"], ess["B"], ess["H"], tau, Tc, Ec, w)

    expected = (theta_w + w * (ess["Sph"] + Ec) / ess["H"] - (ess["Sth"] - (tau - Tc)) / ess["B"]) / (
        beta_w + w**2 / ess["H"] + 1 / ess["B"]
    )
    assert tw.shape == (1, 1)
    assert_allclose(tw.item(), expected, rtol=1e-9)


def test_optimal_work_time_zeroes_the_work_first_order_condition():
    p = aq_params()
    tau = 168
    Tc = np.array([70, 80, 60.0])
    Ec = np.array([150, 200, 120.0])
    w = np.array([12, 15, 20.0])
    theta_w, beta_w = 5, 10
    ess = essentials(p)

    tw = mtuem.utils.get_tw_additive_quadratic(
        theta_w, beta_w, ess["Sth"], ess["Sph"], ess["B"], ess["H"], tau, Tc, Ec, w
    ).ravel()

    mu = (ess["Sth"] - (tau - tw - Tc)) / ess["B"]
    lam = (ess["Sph"] - (w * tw - Ec)) / ess["H"]
    assert_allclose(theta_w - beta_w * tw - mu + w * lam, np.zeros(3), atol=1e-9)


def test_optimal_times_and_goods_satisfy_the_time_and_money_budgets_exactly():
    p = aq_params()
    tau = 168
    Tc = np.array([70, 65.0])
    Ec = np.array([150, 180.0])
    w = np.array([12, 16.0])
    theta_w, beta_w = 0, 8
    ess = essentials(p)
    tw = mtuem.utils.get_tw_additive_quadratic(theta_w, beta_w, ess["Sth"], ess["Sph"], ess["B"], ess["H"], tau, Tc, Ec, w)

    ti = mtuem.utils.get_ti_additive_quadratic(p["times_elasticities"], p["times_satiations"], ess["Sth"], ess["B"], tau, tw, Tc)
    xi = mtuem.utils.get_xi_additive_quadratic(
        p["goods_elasticities"], p["goods_satiations"], p["goods_cost"], ess["Sph"], ess["H"], tw, Ec, w
    )

    assert_allclose(ti.sum(axis=1) + tw.ravel(), tau - Tc, rtol=1e-9)
    assert_allclose(xi @ vals(p["goods_cost"]), w * tw.ravel() - Ec, rtol=1e-9)


def test_additive_quadratic_allocations_decrease_linearly_in_the_shadow_prices():
    p = aq_params()
    tau, Tc, Ec, w = 168, 70, 150, 12
    ess = essentials(p)
    Tw = np.array([[40.0], [45.0]])

    ti = mtuem.utils.get_ti_additive_quadratic(p["times_elasticities"], p["times_satiations"], ess["Sth"], ess["B"], tau, Tw, Tc)
    xi = mtuem.utils.get_xi_additive_quadratic(
        p["goods_elasticities"], p["goods_satiations"], p["goods_cost"], ess["Sph"], ess["H"], Tw, Ec, w
    )

    theta, beta = vals(p["times_elasticities"]), vals(p["times_satiations"])
    phi, eta, P = vals(p["goods_elasticities"]), vals(p["goods_satiations"]), vals(p["goods_cost"])
    mu = (ess["Sth"] - (tau - Tw.ravel() - Tc)) / ess["B"]
    lam = (ess["Sph"] - (w * Tw.ravel() - Ec)) / ess["H"]

    assert ti.shape == (2, 2)
    assert_allclose(ti, np.outer(np.ones(2), theta / beta) - np.outer(mu, 1 / beta), rtol=1e-12)
    assert_allclose(xi, np.outer(np.ones(2), phi / eta) - np.outer(lam, P / eta), rtol=1e-12)


def test_values_of_time_equal_the_shadow_price_ratio_and_satisfy_vol_equals_w_plus_vtaw():
    p = aq_params()
    tau = 168
    Tc = np.array([75, 70.0])
    Ec = np.array([180, 150.0])
    w = np.array([14, 18.0])
    theta_w, beta_w = 0, 11
    ess = essentials(p)
    tw = mtuem.utils.get_tw_additive_quadratic(theta_w, beta_w, ess["Sth"], ess["Sph"], ess["B"], ess["H"], tau, Tc, Ec, w)

    vot = mtuem.utils.get_values_of_time_additive_quadratic(ess["Sth"], ess["Sph"], ess["B"], ess["H"], tw, tau, Tc, Ec, w)

    mu = (ess["Sth"] - (tau - tw.ravel() - Tc)) / ess["B"]
    lam = (ess["Sph"] - (w * tw.ravel() - Ec)) / ess["H"]
    assert np.all(np.isfinite(vot))
    assert_allclose(vot[:, 0], mu / lam, rtol=1e-12)
    assert_allclose(vot[:, 0], w + vot[:, 1], rtol=1e-12)
