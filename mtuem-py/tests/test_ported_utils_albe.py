"""Port of tests/testthat/test_utils_albe.R."""

import numpy as np
from numpy.testing import assert_allclose

import mtuem


def test_get_tw_albe_returns_the_root_of_the_alpha_beta_work_equation():
    tau = 168
    Tc = np.array([50, 60, 72.0])
    Ec = np.array([90, 120, 140.0])
    w = np.array([14, 18, 21.0])
    we = {"alpha": 0.3, "beta": 0.25}

    tw = mtuem.utils.get_tw_albe(we, tau, Tc, Ec, w)

    ec = Ec / (w * (tau - Tc))
    B = we["beta"] + we["alpha"] * ec
    u = tw.ravel() / (tau - Tc)
    assert_allclose(u**2 - 2 * B * u + ec * (2 * we["alpha"] + 2 * we["beta"] - 1), np.zeros(3), atol=1e-12)


def test_get_tw_albe_returns_an_n_by_1_array_and_a_feasible_interior_allocation():
    tau = 168
    rng = np.random.default_rng(9)
    Tc = rng.uniform(40, 90, 40)
    Ec = rng.uniform(50, 180, 40)
    w = rng.uniform(10, 28, 40)
    we = {"alpha": 0.3, "beta": 0.25}

    tw = mtuem.utils.get_tw_albe(we, tau, Tc, Ec, w)

    assert tw.shape == (40, 1)
    assert np.all(tw > 0)
    assert np.all(tw.ravel() < tau - Tc)
    assert np.all(w * tw.ravel() - Ec > 0)


def test_get_ti_albe_splits_residual_time_by_elasticity_over_one_minus_two_beta():
    tau = 168
    Tc = np.array([60, 70.0])
    beta = 0.25
    te = {"th1": 0.2, "th2": 0.1}
    Tw = np.array([[40.0], [38.0]])

    ti = mtuem.utils.get_ti_albe(te, beta, Tw, tau, Tc)

    assert ti.shape == (2, 2)
    assert_allclose(ti[:, 0], (0.2 / (1 - 2 * beta)) * (tau - Tw.ravel() - Tc), rtol=1e-12)
    assert_allclose(ti[:, 1] / ti[:, 0], np.full(2, 0.5), rtol=1e-12)


def test_get_xi_albe_splits_residual_income_by_elasticity_over_price_and_one_minus_alpha():
    alpha = 0.3
    ge = {"ph1": 0.2, "ph2": 0.1}
    gc = {"p1": 1, "p2": 2}
    Tw = np.array([[42.0], [39.0]])
    Ec = np.array([100, 90.0])
    w = np.array([16, 19.0])

    xj = mtuem.utils.get_xi_albe(ge, gc, alpha, Tw, Ec, w)

    income = w * Tw.ravel() - Ec
    assert_allclose(xj[:, 0], (0.2 / 1 / (1 - alpha)) * income, rtol=1e-12)
    assert_allclose(xj[:, 1], (0.1 / 2 / (1 - alpha)) * income, rtol=1e-12)


def test_get_values_of_time_albe_satisfies_vol_equals_w_plus_vtaw_and_its_closed_form():
    tau = 168
    Tc = np.array([60, 55.0])
    Ec = np.array([100, 130.0])
    w = np.array([15, 20.0])
    we = {"alpha": 0.3, "beta": 0.25}
    tw = mtuem.utils.get_tw_albe(we, tau, Tc, Ec, w)

    vot = mtuem.utils.get_values_of_time_albe(tw, we, Tc, Ec, w)

    # The implementation hardcodes 168 in place of tau, as in R.
    expected = ((w * tw.ravel() - Ec) / (168 - tw.ravel() - Tc)) * (1 - 2 * we["beta"]) / (1 - 2 * we["alpha"])
    assert_allclose(vot[:, 0], expected, rtol=1e-12)
    assert_allclose(vot[:, 0], w + vot[:, 1], rtol=1e-12)
