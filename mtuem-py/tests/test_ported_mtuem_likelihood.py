"""Port of tests/testthat/test_mtuem_likelihood.R."""

import numpy as np
import pytest
from scipy import stats
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import cov2cor, mtuem_chol_vec, mtuem_dmvnorm_ref, mtuem_fixture_thph

SIGMA3 = np.array([[4.0, 1.2, -0.6], [1.2, 2.25, 0.4], [-0.6, 0.4, 6.25]])


def test_preprocess_returns_the_component_name_falling_back_to_component_name2():
    fx = mtuem_fixture_thph(n=10)
    named = dict(fx.mtuem_settings, componentName="work_component")
    alt = dict(fx.mtuem_settings, componentName2="fallback_component")

    assert mtuem.likelihood.theta_phi(named, "preprocess", fx.inputs)["componentName"] == "work_component"
    assert mtuem.likelihood.theta_phi(alt, "preprocess", fx.inputs)["componentName"] == "fallback_component"
    default = mtuem.likelihood.theta_phi(fx.mtuem_settings, "preprocess", fx.inputs)
    assert default["componentName"] == ".."
    assert default["gradient"] is False


def test_validate_and_zero_ll_return_the_shapes_the_engine_expects():
    fx = mtuem_fixture_thph(n=17)

    v = mtuem.likelihood.theta_phi(fx.mtuem_settings, "validate", fx.inputs)
    z = mtuem.likelihood.theta_phi(fx.mtuem_settings, "zero_LL", fx.inputs)

    assert_allclose(v, np.ones(17))
    assert z.shape == (17,)
    assert np.all(np.isnan(z))


def test_an_unknown_functionality_is_rejected():
    fx = mtuem_fixture_thph(n=5)

    with pytest.raises(ValueError, match="Invalid value of argument"):
        mtuem.likelihood.theta_phi(fx.mtuem_settings, "not_a_functionality", fx.inputs)


def test_single_equation_likelihood_equals_the_univariate_normal_density():
    fx = mtuem_fixture_thph(n=40, times_elasticities={}, seed=31)
    settings = dict(fx.mtuem_settings, sig={"sig1": 2.3})

    L = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)

    err = fx.obs[:, 0] - fx.opt[:, 0]
    assert_allclose(L, stats.norm.pdf(err, scale=2.3), rtol=1e-12)


def test_with_no_sig_supplied_the_single_equation_sigma_is_inferred_from_the_residuals():
    fx = mtuem_fixture_thph(n=60, times_elasticities={}, seed=32)

    L = mtuem.likelihood.theta_phi(fx.mtuem_settings, "estimate", fx.inputs)

    err = fx.obs[:, 0] - fx.opt[:, 0]
    assert_allclose(L, stats.norm.pdf(err, scale=np.std(err, ddof=1)), rtol=1e-12)


def test_sig_without_rho_gives_the_product_of_independent_normal_densities():
    fx = mtuem_fixture_thph(n=50, seed=33)
    sig = [2.1, 1.9, 2.3]
    settings = dict(fx.mtuem_settings, sig=sig)

    L = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)

    assert_allclose(L, mtuem_dmvnorm_ref(fx.obs - fx.opt, np.diag(np.square(sig))), rtol=1e-10)


def test_sig_and_rho_reproduce_the_full_multivariate_normal_density():
    fx = mtuem_fixture_thph(n=45, seed=34)
    sig = np.array([2.0, 1.5, 2.5])
    # rho is the upper triangle read row-wise: r12, r13, r23.
    r12, r13, r23 = 0.35, -0.2, 0.15
    R = np.array([[1, r12, r13], [r12, 1, r23], [r13, r23, 1]])
    settings = dict(fx.mtuem_settings, sig=list(sig), rho=[r12, r13, r23])

    L = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)

    Sigma = np.diag(sig) @ R @ np.diag(sig)
    assert_allclose(L, mtuem_dmvnorm_ref(fx.obs - fx.opt, Sigma), rtol=1e-10)


def test_a_cholesky_specification_reproduces_the_multivariate_normal_density():
    fx = mtuem_fixture_thph(n=45, seed=35)
    settings = dict(fx.mtuem_settings, cholesky=mtuem_chol_vec(SIGMA3))

    L = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)

    assert_allclose(L, mtuem_dmvnorm_ref(fx.obs - fx.opt, SIGMA3), rtol=1e-10)


def test_a_cholesky_specification_agrees_with_the_equivalent_sig_rho_specification():
    fx = mtuem_fixture_thph(n=45, seed=36)
    R = cov2cor(SIGMA3)
    chol_settings = dict(fx.mtuem_settings, cholesky=mtuem_chol_vec(SIGMA3))
    sigrho_settings = dict(
        fx.mtuem_settings, sig=list(np.sqrt(np.diag(SIGMA3))), rho=[R[0, 1], R[0, 2], R[1, 2]]
    )

    L_chol = mtuem.likelihood.theta_phi(chol_settings, "estimate", fx.inputs)
    L_sig = mtuem.likelihood.theta_phi(sigrho_settings, "estimate", fx.inputs)

    assert_allclose(L_chol, L_sig, rtol=1e-10)


def test_a_cholesky_vector_of_the_wrong_length_is_rejected():
    fx = mtuem_fixture_thph(n=10, seed=37)
    settings = dict(fx.mtuem_settings, cholesky=[1, 0.2, 1])  # 3 entries for a 3-equation model

    with pytest.raises(ValueError, match="must match M"):
        mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)


def test_with_neither_sig_nor_rho_the_covariance_is_inferred_from_the_residuals():
    fx = mtuem_fixture_thph(n=300, seed=38)

    L = mtuem.likelihood.theta_phi(fx.mtuem_settings, "estimate", fx.inputs)

    # sigma from the residual covariance, correlation from the residual correlation
    # (already positive definite here, so the nearest-PD step is the identity).
    err = fx.obs - fx.opt
    sig = np.sqrt(np.diag(np.cov(err, rowvar=False)))
    R = np.corrcoef(err, rowvar=False)
    Sigma = np.diag(sig) @ R @ np.diag(sig)
    assert_allclose(L, mtuem_dmvnorm_ref(err, Sigma), rtol=1e-8)


def test_a_degenerate_sigma_yields_a_zero_likelihood_instead_of_an_error():
    fx1 = mtuem_fixture_thph(n=12, times_elasticities={}, seed=39)
    one_eq = dict(fx1.mtuem_settings, sig=[0])

    L1 = mtuem.likelihood.theta_phi(one_eq, "estimate", fx1.inputs)

    assert_allclose(L1, np.zeros(12))


def test_a_non_positive_definite_correlation_matrix_yields_a_zero_likelihood():
    # r12 = r13 = 0.99 with r23 = -0.99 is not a valid correlation matrix.
    fx = mtuem_fixture_thph(n=20, seed=40)
    settings = dict(fx.mtuem_settings, sig=[2, 2, 2], rho=[0.99, 0.99, -0.99])

    L = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)

    assert_allclose(L, np.zeros(20))


def test_conditionals_raw_and_output_return_the_same_likelihood_as_estimate():
    fx = mtuem_fixture_thph(n=25, seed=41)
    settings = dict(fx.mtuem_settings, sig=[2, 2, 2])

    L_est = mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs)
    for f in ("conditionals", "raw", "output"):
        assert_allclose(mtuem.likelihood.theta_phi(settings, f, fx.inputs), L_est)


def test_report_prints_the_average_values_of_time_and_still_returns_the_likelihood(capsys):
    fx = mtuem_fixture_thph(n=15, seed=42)
    settings = dict(fx.mtuem_settings, sig=[2, 2, 2])

    L = mtuem.likelihood.theta_phi(settings, "report", fx.inputs)

    assert "VoL" in capsys.readouterr().out
    assert_allclose(L, mtuem.likelihood.theta_phi(settings, "estimate", fx.inputs))


def test_prediction_returns_the_optimal_allocations_and_closes_both_budgets():
    fx = mtuem_fixture_thph(
        n=30, times_elasticities={"th1": 0.3, "th2": 0.2}, goods_elasticities={"ph1": 0.4}, seed=43
    )

    pred = mtuem.likelihood.theta_phi(fx.mtuem_settings, "prediction", fx.inputs)

    assert list(pred.columns) == ["Tw", "Tf1", "Tf2", "Tfi", "Ef1", "Xfj", "VoL", "VTAW"]
    assert_allclose(pred["Tw"] + pred["Tf1"] + pred["Tf2"] + pred["Tfi"], fx.tau - fx.Tc, rtol=1e-9)
    assert_allclose(pred["Ef1"] + pred["Xfj"], fx.w * pred["Tw"] - fx.Ec, rtol=1e-9)
    assert_allclose(pred["VoL"], fx.w + pred["VTAW"], rtol=1e-9)


def test_prediction_reports_aggregate_t_and_x_when_no_free_categories_are_given():
    fx = mtuem_fixture_thph(n=20, times_elasticities={}, seed=44)

    pred = mtuem.likelihood.theta_phi(fx.mtuem_settings, "prediction", fx.inputs)

    assert list(pred.columns) == ["Tw", "T", "X", "VoL", "VTAW"]
    assert_allclose(pred["Tw"] + pred["T"], fx.tau - fx.Tc, rtol=1e-9)
    assert_allclose(pred["X"], fx.w * pred["Tw"] - fx.Ec, rtol=1e-9)


def test_prediction_agrees_with_the_residual_computation_used_for_estimation():
    fx = mtuem_fixture_thph(n=20, seed=45)

    pred = mtuem.likelihood.theta_phi(fx.mtuem_settings, "prediction", fx.inputs)
    err = mtuem.likelihood.residuals(fx.mtuem_settings, fx.inputs)

    eqs = ["Tw", "Tf1", "Tf2"]
    assert_allclose(pred[eqs].to_numpy(), fx.database[eqs].to_numpy() - err, atol=1e-10)


def test_the_gradient_functionality_is_a_documented_stub():
    fx = mtuem_fixture_thph(n=10, seed=46)

    g = mtuem.likelihood.theta_phi(fx.mtuem_settings, "gradient", fx.inputs)

    # Gradients are not implemented, so the estimator must fall back to numeric derivatives.
    assert g == {"like": None, "grad": None}
