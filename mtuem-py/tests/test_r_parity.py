"""Parity with the R package: results stored by tests/reference/make_reference.R.

The reference files hold, for each README model, the R log-likelihood at the starting values,
at a perturbed vector and at the R estimates, plus the R estimates, standard errors, values of
time and covariance decodings.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import CASES, CHK, build
from mtuem import apollo

REF = Path(__file__).parent / "reference"
SINGLE = ["eq1", "eq3_profiled", "eq3_sig", "eq3_full", "eq3_cholesky"]


def load_ref(case: str) -> dict:
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


def num(series) -> pd.Series:
    return pd.to_numeric(pd.Series(series), errors="coerce")


def total_ll(case: str, beta: dict) -> float:
    _, _, pf, inputs = build(case)
    return float(np.sum(apollo._loglik_per_indiv(pf(beta, inputs, "estimate")["model"], inputs)))


@pytest.fixture(scope="module")
def single_models():
    out = {}
    for case in SINGLE:
        beta, fixed, pf, inputs = build(case)
        out[case] = (apollo.estimate(beta, fixed, pf, inputs, {"silent": True}), pf, inputs)
    return out


@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("point", ["start", "test", "estimate"])
def test_log_likelihood_matches_r(case, point):
    ref = load_ref(case)
    beta_key = {"start": "beta_start", "test": "beta_test", "estimate": "estimate"}[point]
    ll_key = {"start": "ll_start", "test": "ll_test", "estimate": "ll_at_estimate"}[point]
    assert total_ll(case, ref[beta_key]) == pytest.approx(ref[ll_key], rel=1e-10, abs=1e-8)


@pytest.mark.parametrize("case", SINGLE)
def test_estimation_matches_r(case, single_models):
    ref = load_ref(case)
    model = single_models[case][0]
    assert model.code == 4
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    assert model.eigValue < 0
    est, se, robse = num(ref["estimate"]), num(ref["se"]), num(ref["robse"])
    np.testing.assert_allclose(model.estimate[est.index], est, rtol=1e-4, atol=1e-4)
    free = se.dropna().index
    np.testing.assert_allclose(model.se[free], se[free], rtol=2e-3)
    np.testing.assert_allclose(model.robse[free], robse[free], rtol=2e-3)


def test_lc_standard_errors_match_r_at_the_r_optimum():
    ref = load_ref("lc_cholesky")
    _, fixed, pf, inputs = build("lc_cholesky")
    model = apollo.estimate(ref["estimate"], fixed, pf, inputs, {"silent": True})
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    se, robse = num(ref["se"]).dropna(), num(ref["robse"]).dropna()
    np.testing.assert_allclose(model.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(model.robse[robse.index], robse, rtol=2e-3)


@pytest.mark.parametrize("case", ["lc_profiled", "lc_shared"])
def test_lc_estimation_reaches_at_least_the_r_likelihood(case):
    # R stops early from the README starting values (false/singular convergence); the Python
    # optimiser must reach a likelihood at least as high.
    ref = load_ref(case)
    beta, fixed, pf, inputs = build(case)
    model = apollo.estimate(beta, fixed, pf, inputs, {"silent": True, "hessianRoutine": "none"})
    assert model.maximum >= ref["maximum"] - 1e-3


@pytest.mark.parametrize("case", SINGLE)
def test_values_of_time_match_r(case, single_models):
    ref = load_ref(case)
    model, _, inputs = single_models[case]
    vot = mtuem.mtuem_values_of_time(model, inputs, Ec="EcI", silent=True)["summary"]
    r = pd.DataFrame(ref["vot"])
    np.testing.assert_allclose(vot["mean"], r["mean"], rtol=1e-3)
    np.testing.assert_allclose(vot["s.e."], r["s.e."], rtol=5e-3, atol=2e-4)


@pytest.mark.parametrize("case", ["lc_cholesky", "lc_profiled", "lc_shared"])
def test_lc_values_of_time_means_match_r(case):
    # Posterior weights and means at the R estimates (no standard errors needed).
    ref = load_ref(case)
    _, fixed, pf, inputs = build(case)
    est = pd.Series(ref["estimate"], dtype=float)
    model = apollo.ApolloModel(estimate=est, apollo_beta=est, apollo_fixed=fixed,
                               apollo_control=inputs.apollo_control, maximum=ref["maximum"], LLStart=0.0)
    vot = mtuem.mtuem_values_of_time(model, inputs, pf, nClass=2, Ec="EcI", standard_errors=False, silent=True)
    r = pd.DataFrame(ref["vot"])
    np.testing.assert_allclose(vot["summary"]["mean"], r["mean"], rtol=1e-6)


@pytest.mark.parametrize("case", ["eq3_profiled", "eq3_sig", "eq3_full", "eq3_cholesky"])
def test_residual_correlation_matches_r(case, single_models):
    ref = load_ref(case)
    model, pf, inputs = single_models[case]
    cr = mtuem.mtuem_get_corr(model, pf, inputs)
    np.testing.assert_allclose(cr["corr"].to_numpy(), np.array(ref["corr"]["corr"]), atol=2e-4)
    np.testing.assert_allclose(cr["sigma"].to_numpy(), num(ref["corr"]["sigma"]).to_numpy(), rtol=1e-4)


def test_prediction_matches_r(single_models):
    ref = load_ref("eq3_cholesky")
    model, pf, inputs = single_models["eq3_cholesky"]
    pred = apollo.prediction(model, pf, inputs, {"silent": True}).head(10)
    r = pd.DataFrame(ref["pred_head"])
    for col in ["Tw", "Tf1", "Tfi", "Ef1", "Xfj", "VoL", "VTAW"]:
        np.testing.assert_allclose(pred[col], r[col], rtol=1e-4)


def test_cholesky_delta_method_matches_r(single_models):
    ref = load_ref("eq3_cholesky")
    model = single_models["eq3_cholesky"][0]
    dm = mtuem.mtuem_cholesky_deltaMethod(model, CHK, ["Tw", "Tf1", "Ef1"])
    r = pd.DataFrame(ref["chol_dm"])
    assert list(dm["Expression"]) == list(r["Expression"])
    np.testing.assert_allclose(dm["Value"], r["Value"], rtol=1e-4, atol=2e-4)
    np.testing.assert_allclose(dm["s.e."], r["s.e."], rtol=5e-3, atol=2e-4)


def test_profiled_covariance_matches_r():
    ref = load_ref("lc_profiled")
    _, _, pf, inputs = build("lc_profiled")
    pf(ref["estimate"], inputs, "estimate")
    prof = mtuem.mtuem_get_profiled_sigma()
    for s in range(2):
        np.testing.assert_allclose(prof["Sigma"][s], np.array(ref["profiled"]["Sigma"][s]), rtol=1e-8)
    np.testing.assert_allclose(prof["shares"], ref["profiled"]["shares"], rtol=1e-10)


def _r_block(case: str, header: str) -> list[str]:
    lines = (REF / f"{case}_R_output.txt").read_text(encoding="utf-8").splitlines()
    start = lines.index(header) + 1
    end = lines.index("", start)
    return lines[start:end]


@pytest.mark.parametrize("case", CASES)
def test_estimates_table_prints_like_apollo(case):
    ref = load_ref(case)
    est = num(ref["estimate"])
    model = apollo.ApolloModel(estimate=est, apollo_beta=est, apollo_fixed=ref["fixed"],
                               apollo_control=apollo.ApolloControl(), maximum=0.0, LLStart=0.0,
                               se=num(ref["se"]), robse=num(ref["robse"]))
    r_table = _r_block(case, "Estimates:")
    if "Rob.s.e." not in r_table[0]:
        pytest.skip("R could not compute a Hessian here and printed BHHH standard errors instead")
    table = apollo.apollo_print_matrix(apollo._estimates_table(model)).splitlines()
    assert table == r_table


@pytest.mark.parametrize("case", CASES)
def test_density_mode_reproduces_r_too(case):
    # workInLogs = FALSE: likelihood functions return densities, as in R and Apollo.
    ref = load_ref(case)
    _, _, pf, inputs = build(case)
    inputs.work_in_logs = False
    L = np.asarray(pf(ref["estimate"], inputs, "estimate")["model"], dtype=float)
    assert np.all(L >= 0)
    assert float(np.sum(np.log(L))) == pytest.approx(ref["ll_at_estimate"], rel=1e-10, abs=1e-8)
