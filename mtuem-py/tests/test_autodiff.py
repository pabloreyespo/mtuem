"""The JAX autodiff backend: same likelihoods as R and the numpy code, exact derivatives."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("jax")

import mtuem  # noqa: E402
from cases import build  # noqa: E402
from mtuem import apollo, numderiv  # noqa: E402
from mtuem import autodiff as ad  # noqa: E402

REF = Path(__file__).parent / "reference"
AD_CASES = ["eq1", "eq3_profiled", "eq3_sig", "eq3_full", "eq3_cholesky", "lc_cholesky", "lc_shared"]
SINGLE = ["eq1", "eq3_profiled", "eq3_sig", "eq3_full", "eq3_cholesky"]


def load_ref(case):
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


def num(x):
    return pd.to_numeric(pd.Series(x), errors="coerce")


@pytest.mark.parametrize("case", AD_CASES)
@pytest.mark.parametrize("point", [("beta_start", "ll_start"), ("beta_test", "ll_test"), ("estimate", "ll_at_estimate")])
def test_autodiff_log_likelihood_matches_r(case, point):
    ref = load_ref(case)
    _, _, pf, inputs = build(case, autodiff=True)
    ll = float(np.sum(np.asarray(pf(ref[point[0]], inputs, "estimate")["model"])))
    assert ll == pytest.approx(ref[point[1]], rel=1e-10, abs=1e-8)


@pytest.mark.parametrize("case", SINGLE)
def test_autodiff_estimation_matches_r(case):
    ref = load_ref(case)
    beta, fixed, pf, inputs = build(case, autodiff=True)
    model = apollo.estimate(beta, fixed, pf, inputs, {"silent": True})
    assert model.code == 4 and model.eigValue < 0
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    se = num(ref["se"]).dropna()
    np.testing.assert_allclose(model.estimate[se.index], num(ref["estimate"])[se.index], rtol=2e-3, atol=2e-3)
    np.testing.assert_allclose(model.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(model.robse[se.index], num(ref["robse"])[se.index], rtol=2e-3)


def test_autodiff_lc_standard_errors_match_r_at_the_r_optimum():
    ref = load_ref("lc_cholesky")
    _, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    model = apollo.estimate(ref["estimate"], fixed, pf, inputs, {"silent": True})
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    se = num(ref["se"]).dropna()
    np.testing.assert_allclose(model.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(model.robse[se.index], num(ref["robse"])[se.index], rtol=2e-3)


def test_exact_gradient_and_hessian_agree_with_numderiv():
    ref = load_ref("lc_cholesky")
    beta, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    ll = apollo._LogLike(pd.Series(ref["beta_test"], dtype=float), fixed, pf, inputs)
    x = ll.base[ll.free].to_numpy()
    g_exact = ll.scores(x).sum(axis=0)
    g_num = numderiv.gradient(ll.total, x)
    np.testing.assert_allclose(g_exact, g_num, rtol=1e-5, atol=1e-4)
    H_exact = ll.hessian(x)
    H_num = numderiv.hessian(ll.total, x)
    np.testing.assert_allclose(H_exact, H_num, rtol=1e-4, atol=1e-2)


def test_scores_stay_finite_when_a_class_probability_underflows():
    # At R's lc_shared estimate one class has probability exactly 0 for some people.
    ref = load_ref("lc_shared")
    _, fixed, pf, inputs = build("lc_shared", autodiff=True)
    assert float(np.min(inputs.attach(ref["estimate"]).pi_values[0])) == 0.0
    ll = apollo._LogLike(pd.Series(ref["estimate"], dtype=float), fixed, pf, inputs)
    assert np.all(np.isfinite(ll.scores(ll.base[ll.free].to_numpy())))


def test_posteriors_prediction_and_values_of_time_match_numpy_backend():
    ref = load_ref("lc_cholesky")
    est = pd.Series(ref["estimate"], dtype=float)
    out = {}
    for flag in (False, True):
        _, fixed, pf, inputs = build("lc_cholesky", autodiff=flag)
        model = apollo.ApolloModel(estimate=est, apollo_beta=est, apollo_fixed=fixed,
                                   apollo_control=inputs.apollo_control, maximum=0.0, LLStart=0.0)
        out[flag] = (
            apollo.lc_conditionals(model, pf, inputs)[["X1", "X2"]].to_numpy(),
            apollo.prediction(model, pf, inputs, {"silent": True}).drop(columns=["ID", "Observation"]).to_numpy(),
            mtuem.mtuem_values_of_time(model, inputs, pf, nClass=2, Ec="EcI", standard_errors=False, silent=True)["summary"]["mean"].to_numpy(),
        )
    for a, b in zip(out[False], out[True]):
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12)


def test_em_and_multistart_run_with_the_autodiff_backend():
    beta, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    em = apollo.lc_em(beta, fixed, pf, inputs, {"EMmaxIterations": 2, "postEM": 0, "silent": True})
    assert em.maximum >= em.LLStart
    beta, fixed, pf, inputs = build("eq1", autodiff=True)
    res = mtuem.customMultiStart(beta, fixed, pf, inputs, {"nCandidates": 2}, {"maxIterations": 100}, seed=1,
                                 normalization={"normalization": "theta_phi", "work_elasticities": ["Theta", "Phi", "thw"]})
    assert res["best_model"] is not None
    assert res["best_ll"] == pytest.approx(load_ref("eq1")["maximum"], abs=1e-3)


def test_degenerate_error_structures_give_minus_infinity():
    _, _, _, inputs = build("eq3_sig", autodiff=True)
    b = inputs.attach(load_ref("eq3_sig")["estimate"])
    s = {"work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"],
         "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
         "times_elasticities": [b.theta_1], "goods_elasticities": [b.phi_1],
         "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168}
    assert np.all(np.isneginf(np.asarray(ad.mtuem_loglik(s | {"sig": [-5.0, 7.0, 40.0]}, "estimate", inputs))))
    assert np.all(np.isneginf(np.asarray(ad.mtuem_loglik(s | {"sig": [5.0, 7.0, 40.0], "rho": [0.9, 0.9, -0.9]}, "estimate", inputs))))
    with pytest.raises(ValueError):
        ad.mtuem_loglik(s | {"cholesky": [1.0, 0.0, 1.0]}, "estimate", inputs)
    assert np.all(np.asarray(ad.mtuem_loglik(s, "validate", inputs)) == 0)
    assert np.all(np.isnan(np.asarray(ad.mtuem_loglik(s, "zero_LL", inputs))))
