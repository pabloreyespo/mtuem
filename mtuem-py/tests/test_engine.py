"""Tests of the estimation engine: EM, output files, latent class helpers and multi-start."""

from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import CHK, build
from mtuem import _engine as engine
from mtuem import numderiv
from mtuem._engine.estimate import _escape_saddle
from mtuem._engine.latent import _loglik_per_indiv
from mtuem._engine.loglike import _hessian
from mtuem._engine.optimize import _maximise
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import LatentClassMTUEM


def small_lc(n=250, shared=False):
    """Two-class Cholesky model on the first n maed rows, classes started apart."""
    db = mtuem.maed()
    db = db[db["EcI"] > 0].iloc[:n].reset_index(drop=True)
    beta = {}
    for s, (phi, chk) in enumerate(((0.3, [6, -5, 4, 10, -20, 40]), (0.35, [6, -4, 5, 6, -12, 24])), start=1):
        beta |= {f"Theta_{s}": 1, f"Phi_{s}": phi, f"thw_{s}": 0.1 * s, f"th1_{s}": 0.7, f"ph1_{s}": 0.15}
        if not shared:
            beta |= {f"{k}_{s}": v for k, v in zip(CHK, chk)}
        beta |= {f"asc_{s}": 0, f"x_female_{s}": 0}
    if shared:
        beta |= dict(zip(CHK, [6, -5, 4, 10, -20, 40]))
    fixed = ["Theta_1", "Theta_2", "asc_1", "x_female_1"]

    def lc_pars(b, inputs):
        names = ["Theta", "Phi", "thw", "th1", "ph1"] + ([] if shared else CHK)
        out = {p: [getattr(b, f"{p}_1"), getattr(b, f"{p}_2")] for p in names}
        V = {"class_1": b.asc_1 + b.female * b.x_female_1, "class_2": b.asc_2 + b.female * b.x_female_2}
        out["pi_values"] = engine.class_alloc({"classes": {"class_1": 1, "class_2": 2}, "utilities": V})
        return out

    def pf(params, inputs, functionality="estimate"):
        b = inputs.attach(params)
        P = {}
        for s in range(len(b.pi_values)):
            st = {
                "work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"],
                "work_elasticities": {"Theta": b.Theta[s], "Phi": b.Phi[s], "thw": b.thw[s]},
                "times_elasticities": [b.th1[s]], "goods_elasticities": [b.ph1[s]],
                "cholesky": [getattr(b, k) for k in CHK] if shared else [getattr(b, k)[s] for k in CHK],
                "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
            }
            P[f"class_{s + 1}"] = mtuem.likelihood.theta_phi(st, functionality, inputs)
        P["model"] = engine.lc({"inClassProb": P, "classProb": b.pi_values}, inputs, functionality)
        return P

    control = Control(model_name="small-lc", indiv_id="PeID", output_dir=Path("output"))
    inputs = engine.make_inputs(db, beta, fixed, control, lc_pars)
    return beta, fixed, pf, inputs


@pytest.fixture(scope="module")
def eq1_model():
    beta, fixed, pf, inputs = build("eq1")
    return engine.estimate(beta, fixed, pf, inputs, silent=True), pf, inputs


@pytest.mark.slow
def test_em_never_decreases_the_likelihood_and_hands_over_to_ml():
    beta, fixed, pf, inputs = small_lc()
    em_only = engine.lc_em(beta, fixed, pf, inputs, max_iterations=3, post_em=0, estimate_kwargs={"silent": True})
    assert em_only.maximum >= em_only.LLStart
    assert em_only.hessian is None and em_only.se.isna().all()

    full = engine.lc_em(beta, fixed, pf, inputs, max_iterations=3, estimate_kwargs={"silent": True})
    assert full.maximum >= em_only.maximum - 1e-6
    assert full.routine.startswith("EM algorithm (bfgs) -> Maximum likelihood")
    assert "(EM)" in str(full.nIter)
    assert full.se[[n for n in full.estimate.index if n not in fixed]].notna().all()


def test_em_accepts_parameters_shared_across_classes():
    beta, fixed, pf, inputs = small_lc(shared=True)
    model = engine.lc_em(beta, fixed, pf, inputs, max_iterations=2, post_em=1, estimate_kwargs={"silent": True})
    assert model.maximum >= model.LLStart
    assert model.se[CHK].notna().all()


def test_lc_conditionals_are_posterior_probabilities():
    beta, fixed, pf, inputs = small_lc()
    model = engine.Estimate(estimate=pd.Series(beta, dtype=float), beta_start=pd.Series(beta, dtype=float),
                            fixed=fixed, control=inputs.control, maximum=0.0, LLStart=0.0)
    post = engine.lc_conditionals(model, pf, inputs)
    assert list(post.columns) == ["ID", "X1", "X2"]
    np.testing.assert_allclose(post[["X1", "X2"]].sum(axis=1), 1.0)


def test_lc_prediction_returns_per_class_allocations_and_shares():
    beta, fixed, pf, inputs = small_lc()
    model = engine.Estimate(estimate=pd.Series(beta, dtype=float), beta_start=pd.Series(beta, dtype=float),
                            fixed=fixed, control=inputs.control, maximum=0.0, LLStart=0.0)
    pred = engine.prediction(model, pf, inputs, silent=True)
    for col in ["ID", "Observation", "Tw_class1", "Tw_class2", "VoL_class1", "prob_class1", "prob_class2"]:
        assert col in pred.columns
    np.testing.assert_allclose(pred["prob_class1"] + pred["prob_class2"], 1.0)


def test_lc_bookkeeping_functionalities():
    _, _, _, inputs = small_lc()
    s = {"inClassProb": {"class_1": np.ones(inputs.nObs), "class_2": np.ones(inputs.nObs)},
         "classProb": [np.full(inputs.nObs, 0.5)] * 2}
    assert np.all(engine.lc(s, inputs, "validate") == 0)  # log 1: the engine always works in logs
    assert np.all(np.isnan(engine.lc(s, inputs, "zero_LL")))
    assert engine.lc(s, inputs, "preprocess")["componentName"] == "model"
    for removed in ("gradient", "report"):
        with pytest.raises(ValueError, match="Invalid value of argument"):
            engine.lc(s, inputs, removed)
    with pytest.raises(ValueError):
        engine.lc({"inClassProb": [np.ones(3)], "classProb": [np.ones(3)] * 2}, inputs)


def test_panel_data_multiplies_rows_within_individuals():
    db = pd.DataFrame({"id": [1, 1, 2], "x": [0.0, 1.0, 2.0]})
    logs = engine.make_inputs(db, {"a": 0.0}, [], Control(indiv_id="id"))
    assert logs.panel and logs.nIndivs == 2
    ll = _loglik_per_indiv(np.log([0.5, 0.5, 0.25]), logs)  # rows summed in logs
    np.testing.assert_allclose(ll, np.log([0.25, 0.25]))


def test_logs_survive_densities_that_underflow():
    # A latent class likelihood whose class densities are far below 1e-308: finite in logs.
    db = pd.DataFrame({"id": [1, 2], "x": [0.0, 0.0]})
    log_l1, log_l2 = np.array([-900.0, -2.0]), np.array([-950.0, -3.0])
    pis = engine.class_alloc({"utilities": {"c1": np.zeros(2), "c2": np.zeros(2)}})
    logs = engine.make_inputs(db, {"a": 0.0}, [], Control(indiv_id="id"))
    mixed = engine.lc({"inClassProb": [log_l1, log_l2], "classProb": pis}, logs)
    expected = np.logaddexp(np.log(0.5) + log_l1, np.log(0.5) + log_l2)
    np.testing.assert_allclose(mixed, expected)


def test_output_text_layout(eq1_model):
    model = eq1_model[0]
    out = engine.output_text(model, print_covar=True, print_corr=True)
    assert "apollo" not in out.lower()
    assert f"mtuem {version('mtuem')}" in out.splitlines()[0]
    for line in [
        "Model name                                  : eq1",
        "Estimation diagnosis                        : successful convergence",
        "Optimisation diagnosis                      : Maximum found",
        "LL(final)                                   : -2333.53",
        "AIC                                         :  4671.06 ",
        "Estimated parameters                        : 2",
        "Classical covariance matrix:",
        "Robust correlation matrix:",
    ]:
        assert line in out


def test_save_files_writes_the_output_files(eq1_model, tmp_path):
    model = eq1_model[0]
    files = engine.save_files(model, tmp_path, model.control.model_name)
    for key in ("output", "estimates", "covar", "robcovar", "corr", "robcorr"):
        assert files[key].exists()
    est = pd.read_csv(files["estimates"], index_col=0)
    assert list(est.columns) == ["Estimate", "Std.err.", "t-ratio(0)", "Rob.std.err.", "Rob.t-ratio(0)"]
    assert est.loc["Phi", "Estimate"] == pytest.approx(model.estimate["Phi"])


def test_delta_method_classical_and_robust(eq1_model):
    model = eq1_model[0]
    expr = {"expression": {"ratio": lambda p: p["thw"] / p["Phi"]}}
    rob = engine.delta_method(model, expr["expression"], silent=True)
    cls = engine.delta_method(model, expr["expression"], varcov="classical", silent=True)
    assert rob["Value"].iloc[0] == cls["Value"].iloc[0]
    assert rob["s.e."].iloc[0] != cls["s.e."].iloc[0]


def test_converged_run_reports_code_zero(eq1_model):
    model = eq1_model[0]
    assert model.code == 0 and model.message == "successful convergence"
    assert model.successfulEstimation and model.routine == "bfgs"


def test_iteration_limit_reports_a_non_converged_code():
    beta, fixed, pf, inputs = build("eq1")
    model = engine.estimate(beta, fixed, pf, inputs, max_iterations=1, hessian="none", silent=True)
    assert model.code == 10 and not model.successfulEstimation
    assert model.message == "iteration limit reached without convergence"


def test_invalid_inputs_are_rejected():
    db = mtuem.maed()
    with pytest.raises(ValueError):
        engine.make_inputs(db, {"a": 1}, ["b"], Control(indiv_id="PeID"))
    with pytest.raises(ValueError):
        engine.make_inputs(db, {"a": 1}, [], Control(indiv_id="nope"))
    with pytest.raises(ValueError):
        engine.make_inputs(db, {"a": float("nan")}, [], Control(indiv_id="PeID"))
    beta, fixed, pf, inputs = build("eq1")
    with pytest.raises(ValueError):
        engine.estimate(beta, fixed, pf, inputs, hessian="analytic", silent=True)  # needs the autodiff backend

    with pytest.raises(ValueError, match="hessian"):
        engine.estimate(beta, fixed, pf, inputs, hessian="foo", silent=True)
    with pytest.raises(ValueError, match="thw"):
        engine.estimate(beta, fixed, pf, inputs, bounds={"thw": (0.0, 2.0)}, silent=True)  # start on the bound

    def zero_pf(params, inputs, functionality="estimate"):
        return {"model": np.full(inputs.nObs, -np.inf)}  # zero likelihood, in logs

    with pytest.raises(FloatingPointError):
        engine.estimate(beta, fixed, zero_pf, inputs, silent=True)

    def no_model_pf(params, inputs, functionality="estimate"):
        return {"class_1": np.zeros(inputs.nObs)}

    with pytest.raises(ValueError, match="model"):
        engine.estimate(beta, fixed, no_model_pf, inputs, silent=True)


def test_saddle_escape_steps_off_a_stationary_saddle():
    # f(x, y) = -x^2 + y^2 - y^4 has a saddle at the origin and maxima at y = +-1/sqrt(2).
    class Toy:
        @staticmethod
        def total(z):
            return -z[0] ** 2 + z[1] ** 2 - z[1] ** 4

    origin = np.zeros(2)
    res = _maximise(lambda z: np.atleast_1d(Toy.total(z)), origin)
    assert res.f == 0.0  # BFGS alone cannot leave a point with zero gradient
    x = _escape_saddle(Toy, origin, 0.0, np.diag([-2.0, 2.0]))
    assert x is not None and Toy.total(x) > 0
    final = _maximise(lambda z: np.atleast_1d(Toy.total(z)), x)
    assert final.f == pytest.approx(0.25, abs=1e-6)
    assert _escape_saddle(Toy, origin, 0.0, np.diag([-2.0, -2.0])) is None


@pytest.mark.slow
def test_multistart_with_em_warm_start():
    values, _, _, inputs = small_lc(n=150)
    model = LatentClassMTUEM(
        work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], Tc="Tc", Ec="EcI", w="w",
        covariance="cholesky", n_classes=2, allocation=["female"],
    )
    db = inputs.database
    beta = model.beta().set(**values)
    settings = EstimationSettings(
        backend="numpy", max_iterations=30, verbose=False,
        multistart=MultiStart(n=2, first_em=True, em_iter_max=2, normalization="theta_phi", non_saddle=False, seed=3),
    )
    try:
        res = model.fit(db, beta, Control(indiv_id="PeID"), settings)
        attempts = res.starts
    except RuntimeError as err:  # no candidate converged in 30 iterations
        attempts = err.attempts
    assert len(attempts) == 2
    assert any(a is not None for a in attempts)


def test_alpha_beta_normalization_caps_work_elasticities():
    rng = np.random.default_rng(0)
    bm = pd.DataFrame({"alpha": rng.uniform(0, 1, 50), "beta": rng.uniform(0, 1, 50),
                       "t1": rng.uniform(-1, 1, 50), "g1": rng.uniform(-1, 1, 50)})
    from mtuem.multistart import _normalize

    out = _normalize(bm.copy(), rng, list(bm.columns), "alpha_beta", ["alpha", "beta"], ["t1"], ["g1"], [], [], [])
    assert (out[["alpha", "beta"]] < 0.5).all().all()
    assert (out["t1"] > 0).all() and (out["g1"] > 0).all()


def test_hessian_falls_back_to_small_steps_near_a_boundary():
    # numDeriv's first step (10% of x) leaves the region where the function is finite.
    def f(z):
        return -z[0] ** 2 if abs(z[0]) < 1.05 else -np.inf

    x = np.array([1.0])
    with np.errstate(invalid="ignore"):  # the failed first pass is expected to produce nan
        assert not np.all(np.isfinite(numderiv.hessian(f, x)))
        np.testing.assert_allclose(_hessian(f, x), [[-2.0]], rtol=1e-6)
