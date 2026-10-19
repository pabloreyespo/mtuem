"""Port of tests/testthat/test_cholesky_deltaMethod.R."""

import re

import numpy as np
import pytest
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import cov2cor, mtuem_fake_model


def test_a_cholesky_parameter_count_that_is_not_m_times_m_plus_1_over_2_is_rejected():
    model = mtuem_fake_model({"a": 1, "b": 2}, np.eye(2) * 0.01)

    with pytest.raises(ValueError, match=re.escape("M * (M + 1) / 2")):
        mtuem.mtuem_cholesky_deltaMethod(model, ["a", "b"])


def test_delta_method_recovers_the_sigmas_and_correlation_implied_by_the_cholesky_factor():
    est = {"chk11": 2, "chk21": 0.5, "chk22": 1.5}
    model = mtuem_fake_model(est, np.diag([0.01, 0.02, 0.03]))

    res = mtuem.mtuem_cholesky_deltaMethod(model, list(est))

    # Sigma = L L', so sig1 = L11, sig2 = sqrt(L21^2 + L22^2), rho = L21 * L11 / (sig1 * sig2).
    L = np.array([[2, 0], [0.5, 1.5]])
    Sigma = L @ L.T
    sig = np.sqrt(np.diag(Sigma))
    rho = Sigma[1, 0] / (sig[0] * sig[1])

    assert list(res["Expression"]) == ["sig_1", "sig_2", "rho_1_2"]
    assert_allclose(res["Value"].to_numpy(), [*sig, rho], atol=1e-4)
    assert np.all(res["s.e."] > 0)


def test_equation_names_label_the_delta_method_output_when_supplied():
    est = {"chk11": 2, "chk21": 0.5, "chk22": 1.5, "chk31": 0.3, "chk32": 0.2, "chk33": 1.1}
    model = mtuem_fake_model(est, np.diag(np.full(6, 0.01)))

    res = mtuem.mtuem_cholesky_deltaMethod(model, list(est), eq_names=["Tw", "Tf1", "Ef1"])

    assert list(res["Expression"]) == [
        "sig_Tw", "sig_Tf1", "sig_Ef1", "rho_Tw_Tf1", "rho_Tw_Ef1", "rho_Tf1_Ef1",
    ]


def test_delta_method_values_match_the_covariance_implied_by_a_3_equation_factor():
    est = {"chk11": 2, "chk21": 0.5, "chk22": 1.5, "chk31": 0.3, "chk32": 0.2, "chk33": 1.1}
    model = mtuem_fake_model(est, np.diag(np.full(6, 0.01)))
    L = np.array([[2, 0, 0], [0.5, 1.5, 0], [0.3, 0.2, 1.1]])
    Sigma = L @ L.T
    sig = np.sqrt(np.diag(Sigma))
    R = cov2cor(Sigma)

    res = mtuem.mtuem_cholesky_deltaMethod(model, list(est))

    # Sigmas first, then correlations in (1,2), (1,3), (2,3) order.
    assert_allclose(res["Value"].to_numpy(), [*sig, R[1, 0], R[2, 0], R[2, 1]], atol=1e-4)
