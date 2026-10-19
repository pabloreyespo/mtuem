"""Port of tests/testthat/test_lc_covariance.R."""

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import cov2cor, inputs_of, mtuem_dmvnorm_ref
from mtuem import _engine as engine
from mtuem.config import Control


def test_mtuem_residuals_refactor_preserves_mtuem_likelihood_output_closed_form_check():
    # Residuals and likelihood are derived independently from the closed-form Theta-Phi
    # equations, without calling mtuem_residuals itself.
    n, tau = 8, 168
    Tc, Ec, w = np.full(n, 60.0), np.full(n, 100.0), np.full(n, 15.0)
    we = {"Theta": 1, "Phi": 1, "thw": 0}
    te = {"th1": 0.3, "th2": 0.2}

    tw_opt = mtuem.get_tw_thph(we, tau, Tc, Ec, w)
    ti_opt = mtuem.get_ti_thph(te, we["Theta"], tw_opt, tau, Tc)

    eps = np.random.default_rng(2).normal(0, 2, (n, 3))
    database = pd.DataFrame({
        "PeID": np.arange(1, n + 1), "Tw": tw_opt.ravel() + eps[:, 0], "Tf1": ti_opt[:, 0] + eps[:, 1],
        "Tf2": ti_opt[:, 1] + eps[:, 2], "Tc": Tc, "EcI": Ec, "w": w,
    })
    sig = {"sig1": 2.1, "sig2": 1.9, "sig3": 2.3}
    settings = {
        "work_times": "Tw", "free_times": ["Tf1", "Tf2"], "work_elasticities": we, "times_elasticities": te,
        "sig": sig, "Tc": Tc, "Ec": Ec, "w": w, "tau": tau,
    }

    L_actual = mtuem.likelihood.theta_phi(settings, "estimate", inputs_of(database))

    # With explicit sig and no rho the correlation is the identity, so the likelihood is
    # a product of independent normal densities.
    opt = np.column_stack([tw_opt, ti_opt])
    err = database[["Tw", "Tf1", "Tf2"]].to_numpy() - opt
    sig_vec = np.array(list(sig.values()))
    mu = err / sig_vec
    ll = -0.5 * np.sum(mu**2, axis=1) - np.sum(np.log(sig_vec)) - 0.5 * 3 * np.log(2 * np.pi)
    assert_allclose(L_actual, np.exp(ll), rtol=1e-10)


def test_scaled_cholesky_lambda_l_preserves_correlation_and_scales_sigma():
    L = np.array([[2.0, 0, 0], [0.5, 1.5, 0], [0.3, 0.7, 1.2]])
    lam = 1.7

    Sigma = L @ L.T
    Sigma_lam = (lam * L) @ (lam * L).T

    assert_allclose(np.sqrt(np.diag(Sigma_lam)), lam * np.sqrt(np.diag(Sigma)), rtol=1e-12)
    assert_allclose(cov2cor(Sigma_lam), cov2cor(Sigma), rtol=1e-12)


def test_row_scaled_cholesky_diag_d_l_preserves_correlation_structure():
    L = np.array([[2.0, 0, 0], [0.5, 1.5, 0], [0.3, 0.7, 1.2]])
    D = np.diag([1.3, 0.8, 2.1])

    Sigma = L @ L.T
    Sigma_d = D @ Sigma @ D

    assert_allclose(cov2cor(Sigma_d), cov2cor(Sigma), rtol=1e-12)


def test_row_major_cholesky_extraction_and_reconstruction_round_trip_recovers_sigma():
    M = 3
    a = np.array([1, -0.4, 0.6])
    sigw = 2.5
    psi = np.array([1.1, 0.9, 1.3])
    Sigma = sigw**2 * np.outer(a, a) + np.diag(psi**2)
    L = np.linalg.cholesky(Sigma)

    # Row-major lower triangle: L11, L21, L22, L31, L32, L33.
    cholesky = [L[r, c] for r in range(M) for c in range(r + 1)]

    # Rebuilt through the package helper that the likelihoods use.
    Sigma_rebuilt = mtuem.utils.cholesky_to_sigma(cholesky, M)

    assert_allclose(Sigma_rebuilt, Sigma, rtol=1e-10)


def _two_class_settings(N, tau, Tc, Ec, w, we1, we2, te1, te2, p1):
    class_settings = [
        {"work_times": "Tw", "free_times": "Tf1", "work_elasticities": we1, "times_elasticities": te1,
         "Tc": Tc, "Ec": Ec, "w": w, "tau": tau},
        {"work_times": "Tw", "free_times": "Tf1", "work_elasticities": we2, "times_elasticities": te2,
         "Tc": Tc, "Ec": Ec, "w": w, "tau": tau},
    ]
    return class_settings, {"class_settings": class_settings, "classProb": [np.full(N, p1), np.full(N, 1 - p1)]}


def test_profiled_likelihood_recovers_class_specific_sigma_and_matches_the_ordinary_mixture_likelihood():
    rng = np.random.default_rng(123)
    N, tau = 5000, 168
    Tc = rng.uniform(40, 90, N)
    w = rng.uniform(8, 25, N)
    Ec = rng.uniform(50, 200, N)

    we1, we2 = {"Theta": 1, "Phi": 1, "thw": 0}, {"Theta": 1, "Phi": 2.5, "thw": 0.5}
    te1, te2 = {"th1": 0.4}, {"th1": 0.15}
    tw1 = mtuem.get_tw_thph(we1, tau, Tc, Ec, w)
    tw2 = mtuem.get_tw_thph(we2, tau, Tc, Ec, w)
    opt1 = np.column_stack([tw1, mtuem.get_ti_thph(te1, we1["Theta"], tw1, tau, Tc)])
    opt2 = np.column_stack([tw2, mtuem.get_ti_thph(te2, we2["Theta"], tw2, tau, Tc)])

    S1_true = np.array([[4.0, 1.0], [1.0, 2.25]])
    S2_true = np.array([[9.0, -1.5], [-1.5, 3.24]])
    pi1 = 0.6
    z = rng.binomial(1, 1 - pi1, N) + 1
    e1 = rng.normal(size=(N, 2)) @ np.linalg.cholesky(S1_true).T
    e2 = rng.normal(size=(N, 2)) @ np.linalg.cholesky(S2_true).T
    obs = np.where((z == 1)[:, None], opt1 + e1, opt2 + e2)

    database = pd.DataFrame({"PeID": np.arange(1, N + 1), "Tw": obs[:, 0], "Tf1": obs[:, 1], "Tc": Tc, "EcI": Ec, "w": w})
    inputs = inputs_of(database)
    class_settings, lc_settings = _two_class_settings(N, tau, Tc, Ec, w, we1, we2, te1, te2, pi1)

    P_profiled = mtuem.mtuem_lc_profiled_likelihood(lc_settings, "estimate", inputs)

    assert np.all(np.isfinite(P_profiled))
    assert np.all(P_profiled > 0)

    prof = mtuem.mtuem_get_profiled_sigma()
    S1_hat, S2_hat = prof["Sigma"]
    assert np.abs(S1_hat - S1_true).max() / np.abs(S1_true).max() < 0.2
    assert np.abs(S2_hat - S2_true).max() / np.abs(S2_true).max() < 0.2

    # The profiled likelihood equals the ordinary mixture at the recovered Sigma,
    # computed here through inv()/slogdet() rather than the Cholesky-based density.
    E1 = mtuem.likelihood.residuals(class_settings[0], inputs)
    E2 = mtuem.likelihood.residuals(class_settings[1], inputs)
    P_manual = pi1 * mtuem_dmvnorm_ref(E1, S1_hat) + (1 - pi1) * mtuem_dmvnorm_ref(E2, S2_hat)
    assert_allclose(P_profiled, P_manual, rtol=1e-8)


def test_profiled_likelihood_returns_zeros_instead_of_erroring_on_non_finite_residuals():
    # w = 0 forces a division by zero inside get_tw_thph, producing non-finite residuals.
    N, tau = 10, 168
    Tc, Ec, w = np.full(N, 50.0), np.full(N, 100.0), np.zeros(N)
    we = {"Theta": 1, "Phi": 1, "thw": 0}
    database = pd.DataFrame({"PeID": np.arange(1, N + 1), "Tw": np.full(N, 50.0), "Tc": Tc, "EcI": Ec, "w": w})
    cs = {"work_times": "Tw", "work_elasticities": we, "Tc": Tc, "Ec": Ec, "w": w, "tau": tau}
    lc_settings = {"class_settings": [cs, cs], "classProb": [np.full(N, 0.5), np.full(N, 0.5)]}

    P = mtuem.mtuem_lc_profiled_likelihood(lc_settings, "estimate", inputs_of(database))

    assert_allclose(P, np.zeros(N))


def test_profiled_likelihood_returns_zeros_instead_of_erroring_on_a_non_pd_starting_sigma():
    # Identical noise added to both equations makes the residuals perfectly collinear,
    # so the unweighted starting covariance is singular.
    rng = np.random.default_rng(7)
    N, tau = 30, 168
    Tc, Ec, w = rng.uniform(40, 80, N), rng.uniform(50, 150, N), rng.uniform(10, 20, N)
    we, te = {"Theta": 1, "Phi": 1, "thw": 0}, {"th1": 0.3}
    tw_opt = mtuem.get_tw_thph(we, tau, Tc, Ec, w)
    ti_opt = mtuem.get_ti_thph(te, we["Theta"], tw_opt, tau, Tc)
    noise = rng.normal(0, 3, N)
    database = pd.DataFrame({
        "PeID": np.arange(1, N + 1), "Tw": tw_opt.ravel() + noise, "Tf1": ti_opt.ravel() + noise,
        "Tc": Tc, "EcI": Ec, "w": w,
    })
    cs = {"work_times": "Tw", "free_times": "Tf1", "work_elasticities": we, "times_elasticities": te,
          "Tc": Tc, "Ec": Ec, "w": w, "tau": tau}
    lc_settings = {"class_settings": [cs, cs], "classProb": [np.full(N, 0.5), np.full(N, 0.5)]}

    P = mtuem.mtuem_lc_profiled_likelihood(lc_settings, "estimate", inputs_of(database))

    assert_allclose(P, np.zeros(N))


def test_profiled_likelihood_returns_zeros_instead_of_erroring_on_class_collapse():
    # Class 2 has a vanishingly small prior everywhere, so its posterior weight cannot
    # reach the minimum weight guard.
    rng = np.random.default_rng(11)
    N, tau = 40, 168
    Tc, Ec, w = rng.uniform(40, 80, N), rng.uniform(50, 150, N), rng.uniform(10, 20, N)
    we = {"Theta": 1, "Phi": 1, "thw": 0}
    tw_opt = mtuem.get_tw_thph(we, tau, Tc, Ec, w)
    database = pd.DataFrame({
        "PeID": np.arange(1, N + 1), "Tw": tw_opt.ravel() + rng.normal(0, 2, N), "Tc": Tc, "EcI": Ec, "w": w,
    })
    cs = {"work_times": "Tw", "work_elasticities": we, "Tc": Tc, "Ec": Ec, "w": w, "tau": tau}
    lc_settings = {
        "class_settings": [cs, cs], "classProb": [np.full(N, 1 - 1e-8), np.full(N, 1e-8)],
        "em_settings": {"minWeight": 5},
    }

    P = mtuem.mtuem_lc_profiled_likelihood(lc_settings, "estimate", inputs_of(database))

    assert_allclose(P, np.zeros(N))


def test_custom_multistart_normalizes_shared_unsuffixed_covariance_parameters_without_erroring():
    # A 2-class, 1-equation latent class model with one covariance parameter (chk11)
    # SHARED across both classes, i.e. not suffixed with _1/_2 in the start values.
    db = mtuem.maed()
    db = db[db["EcI"] > 0].iloc[:100].reset_index(drop=True)
    beta = {
        "Theta_1": 1, "Phi_1": 1, "thw_1": 0, "asc_1": 0, "x_female_1": 0, "x_older45_1": 0,
        "Theta_2": 1, "Phi_2": 1, "thw_2": 0, "asc_2": 0, "x_female_2": 0, "x_older45_2": 0,
        "chk11": 8,
    }
    fixed = ["Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1"]

    def lc_pars(b, inputs):
        V = {
            "class_1": b.asc_1 + b.female * b.x_female_1 + b.older45 * b.x_older45_1,
            "class_2": b.asc_2 + b.female * b.x_female_2 + b.older45 * b.x_older45_2,
        }
        return {
            "Theta": [b.Theta_1, b.Theta_2], "Phi": [b.Phi_1, b.Phi_2], "thw": [b.thw_1, b.thw_2],
            "pi_values": engine.class_alloc({"classes": {"class_1": 1, "class_2": 2}, "utilities": V}),
        }

    control = Control(model_name="test-shared-cov-normalize", indiv_id="PeID")
    inputs = engine.make_inputs(db, beta, fixed, control, lc_pars)

    def probabilities(params, inputs, functionality="estimate"):
        b = inputs.attach(params)
        P = {}
        for s in range(len(b.pi_values)):
            settings = {
                "work_times": ["Tw"],
                "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
                "cholesky": {"chk11": b.chk11},
                "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168, "componentName": f"class_{s + 1}",
            }
            P[f"class_{s + 1}"] = mtuem.likelihood.theta_phi(settings, functionality, inputs)
        P["model"] = engine.lc({"inClassProb": P, "classProb": b.pi_values}, inputs, functionality)
        return P

    result = mtuem.customMultiStart(
        beta, fixed, probabilities, inputs,
        customMultistart_settings={"nCandidates": 2},
        estimate_settings={"maxIterations": 1, "hessianRoutine": "none", "estimationRoutine": "bfgs"},
        first_em=False,
        verbose=False,
        normalization={
            "normalization": "theta_phi",
            "work_elasticities": ["Theta", "Phi", "thw"],
            "cholesky": ["chk11"],
            "covs": ["asc", "x_female", "x_older45"],
        },
        nClass=2,
        seed=1,
    )

    assert set(result) == {"best_model", "best_ll", "models"}
