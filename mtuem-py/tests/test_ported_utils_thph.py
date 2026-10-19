"""Port of tests/testthat/test_utils_thph.R."""

import numpy as np
from numpy.testing import assert_allclose

import mtuem


def test_get_tw_thph_returns_the_root_of_the_work_first_order_condition():
    tau = 168
    Tc = np.array([50, 65, 70.0])
    Ec = np.array([100, 120, 90.0])
    w = np.array([12, 18, 22.0])
    we = {"Theta": 1.2, "Phi": 0.8, "thw": 0.15}

    tw = mtuem.utils.get_tw_thph(we, tau, Tc, Ec, w)

    ec = Ec / (w * (tau - Tc))
    base = we["Theta"] + we["Phi"] + we["thw"]
    A = we["Phi"] + we["thw"] + (we["Theta"] + we["thw"]) * ec
    u = tw.ravel() / (tau - Tc)
    assert_allclose(base * u**2 - A * u + we["thw"] * ec, np.zeros(3), atol=1e-12)


def test_get_tw_thph_collapses_to_the_closed_form_when_thw_is_zero():
    tau, Tc, Ec, w = 168, 60, 100, 15
    we = {"Theta": 1, "Phi": 2, "thw": 0}
    ec = Ec / (w * (tau - Tc))

    tw = mtuem.utils.get_tw_thph(we, tau, Tc, Ec, w)

    expected = (tau - Tc) * (we["Phi"] + we["Theta"] * ec) / (we["Theta"] + we["Phi"])
    assert_allclose(tw.ravel(), expected, rtol=1e-12)


def test_get_tw_thph_returns_an_n_by_1_array_regardless_of_input_length():
    we = {"Theta": 1, "Phi": 1, "thw": 0}

    tw1 = mtuem.utils.get_tw_thph(we, 168, 60, 100, 15)
    twN = mtuem.utils.get_tw_thph(we, 168, np.full(7, 60.0), np.full(7, 100.0), np.full(7, 15.0))

    assert tw1.shape == (1, 1)
    assert twN.shape == (7, 1)
    assert_allclose(twN.ravel(), np.repeat(tw1.ravel(), 7), rtol=1e-12)


def test_get_tw_thph_leaves_a_strictly_positive_free_time_residual():
    tau = 168
    rng = np.random.default_rng(4)
    Tc = rng.uniform(40, 90, 50)
    Ec = rng.uniform(50, 200, 50)
    w = rng.uniform(8, 30, 50)
    we = {"Theta": 1, "Phi": 1, "thw": 0}

    tw = mtuem.utils.get_tw_thph(we, tau, Tc, Ec, w).ravel()

    assert np.all(tw > 0)
    assert np.all(tw < tau - Tc)
    assert np.all(w * tw - Ec > 0)


def test_get_ti_thph_splits_residual_time_in_proportion_to_the_time_elasticities():
    tau = 168
    Tc = np.array([60, 70.0])
    Theta = 1.5
    te = {"th1": 0.3, "th2": 0.45}
    Tw = np.array([[40.0], [35.0]])

    ti = mtuem.utils.get_ti_thph(te, Theta, Tw, tau, Tc)

    assert ti.shape == (2, 2)
    # R returns a matrix with colnames th1, th2; the Python array is positional.
    assert_allclose(ti[:, 0], (0.3 / Theta) * (tau - Tw.ravel() - Tc), rtol=1e-12)
    assert_allclose(ti[:, 1] / ti[:, 0], np.full(2, 0.45 / 0.3), rtol=1e-12)


def test_get_ti_thph_exhausts_the_time_budget_when_the_elasticities_sum_to_theta():
    tau, Tc, Theta = 168, 60, 1
    te = {"th1": 0.4, "th2": 0.6}
    Tw = np.array([[40.0]])

    ti = mtuem.utils.get_ti_thph(te, Theta, Tw, tau, Tc)

    assert_allclose(ti.sum(), tau - Tw.item() - Tc, rtol=1e-12)


def test_get_xi_thph_splits_residual_income_by_elasticity_over_price():
    Phi = 2
    ge = {"ph1": 0.5, "ph2": 1.0}
    gc = {"p1": 1, "p2": 2}
    Tw = np.array([[40.0], [45.0]])
    Ec = np.array([100, 120.0])
    w = np.array([15, 20.0])

    xj = mtuem.utils.get_xi_thph(ge, gc, Phi, Tw, Ec, w)

    income = w * Tw.ravel() - Ec
    assert xj.shape == (2, 2)
    assert_allclose(xj[:, 0], (0.5 / 1 / Phi) * income, rtol=1e-12)
    assert_allclose(xj[:, 1], (1.0 / 2 / Phi) * income, rtol=1e-12)


def test_get_values_of_time_thph_satisfies_vol_equals_w_plus_vtaw():
    tau = 168
    Tc = np.array([60, 70, 55.0])
    Ec = np.array([100, 150, 80.0])
    w = np.array([15, 12, 25.0])
    we = {"Theta": 1.4, "Phi": 0.9, "thw": 0.1}
    tw = mtuem.utils.get_tw_thph(we, tau, Tc, Ec, w)

    vot = mtuem.utils.get_values_of_time_thph(tw, we, Tc, Ec, w)

    # Columns are positional (VoL, VTAW).
    assert vot.shape[1] == 2
    assert_allclose(vot[:, 0], w + vot[:, 1], rtol=1e-12)


def test_get_values_of_time_thph_equals_the_marginal_rate_of_substitution_formula():
    tau, Tc, Ec, w = 168, 60, 100, 15
    we = {"Theta": 1.4, "Phi": 0.9, "thw": 0.1}
    tw = mtuem.utils.get_tw_thph(we, tau, Tc, Ec, w)

    vot = mtuem.utils.get_values_of_time_thph(tw, we, Tc, Ec, w)

    # The implementation hardcodes 168 in place of tau, as in R.
    expected = (we["Theta"] / we["Phi"]) * (w * tw.ravel() - Ec) / (168 - tw.ravel() - Tc)
    assert_allclose(vot[:, 0], expected, rtol=1e-12)


def test_vol_is_increasing_in_theta_over_phi_holding_the_allocation_fixed():
    Tc, Ec, w = 60, 100, 15
    tw = np.array([[45.0]])

    low = mtuem.utils.get_values_of_time_thph(tw, {"Theta": 1, "Phi": 2, "thw": 0}, Tc, Ec, w)
    high = mtuem.utils.get_values_of_time_thph(tw, {"Theta": 2, "Phi": 1, "thw": 0}, Tc, Ec, w)

    assert high[0, 0] > low[0, 0]
