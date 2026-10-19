"""Port of tests/testthat/test_get_cond_err.R."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import cov2cor


def make_rho(M, seed=1):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(M, M))
    return cov2cor(A @ A.T + 2 * np.eye(M))


@pytest.mark.parametrize("M", range(2, 7))
def test_get_cond_err_reproduces_the_multivariate_quadratic_form_and_log_determinant(M):
    rho = make_rho(M, seed=M)
    mu = np.random.default_rng(100 + M).normal(size=(10, M))

    cond = mtuem.get_cond_err(mu, rho)

    assert_allclose(np.sum(cond["cond_mu"] ** 2, axis=1), np.sum((mu @ np.linalg.inv(rho)) * mu, axis=1), rtol=1e-8)
    assert_allclose(np.sum(np.log(cond["cond_sd"])), np.linalg.slogdet(rho)[1], rtol=1e-8, atol=1e-10)


def test_get_cond_err_leaves_standardised_errors_untouched_under_independence():
    M = 4
    mu = np.random.default_rng(1).normal(size=(7, M))

    cond = mtuem.get_cond_err(mu, np.eye(M))

    assert_allclose(cond["cond_mu"], mu, rtol=1e-12)
    assert_allclose(cond["cond_sd"], np.ones(M), rtol=1e-12)


def test_get_cond_err_matches_the_textbook_bivariate_conditioning():
    r = 0.6
    rho = np.array([[1, r], [r, 1]])
    mu = np.array([[1.2, 0.9], [-0.4, 2.1]])

    cond = mtuem.get_cond_err(mu, rho)

    assert_allclose(cond["cond_sd"], [1, 1 - r**2], rtol=1e-12)
    assert_allclose(cond["cond_mu"][:, 0], mu[:, 0], rtol=1e-12)
    assert_allclose(cond["cond_mu"][:, 1], (mu[:, 1] - r * mu[:, 0]) / np.sqrt(1 - r**2), rtol=1e-12)


def test_get_cond_err_preserves_the_shape_of_the_standardised_errors():
    # R also preserves colnames; numpy arrays carry no names, so only the shape is checked.
    rho = make_rho(3, seed=21)
    mu = np.random.default_rng(2).normal(size=(5, 3))

    cond = mtuem.get_cond_err(mu, rho)

    assert cond["cond_mu"].shape == mu.shape


def test_get_cond_err_is_invariant_to_the_equation_ordering_in_the_density_it_implies():
    M = 3
    rho = make_rho(M, seed=5)
    mu = np.random.default_rng(3).normal(size=(9, M))
    perm = [2, 0, 1]

    a = mtuem.get_cond_err(mu, rho)
    b = mtuem.get_cond_err(mu[:, perm], rho[np.ix_(perm, perm)])

    ll_a = -0.5 * np.sum(a["cond_mu"] ** 2, axis=1) - 0.5 * np.sum(np.log(a["cond_sd"]))
    ll_b = -0.5 * np.sum(b["cond_mu"] ** 2, axis=1) - 0.5 * np.sum(np.log(b["cond_sd"]))
    assert_allclose(ll_a, ll_b, rtol=1e-8)
