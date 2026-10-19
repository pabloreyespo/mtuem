"""The JAX backend (likelihood_jax): same likelihoods as R and the numpy code, exact derivatives."""

import json
import warnings
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import build, database
from mtuem import _engine as engine
from mtuem import likelihood as lk
from mtuem import likelihood_jax as lj
from mtuem import numderiv
from mtuem._engine import FallbackWarning
from mtuem._engine.loglike import _LogLike
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import MTUEM

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
    model = engine.estimate(beta, fixed, pf, inputs, silent=True)
    assert model.code == 0 and model.eigValue < 0
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    se = num(ref["se"]).dropna()
    np.testing.assert_allclose(model.estimate[se.index], num(ref["estimate"])[se.index], rtol=2e-3, atol=2e-3)
    np.testing.assert_allclose(model.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(model.robse[se.index], num(ref["robse"])[se.index], rtol=2e-3)


def test_autodiff_lc_standard_errors_match_r_at_the_r_optimum():
    ref = load_ref("lc_cholesky")
    _, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    model = engine.estimate(ref["estimate"], fixed, pf, inputs, silent=True)
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    se = num(ref["se"]).dropna()
    np.testing.assert_allclose(model.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(model.robse[se.index], num(ref["robse"])[se.index], rtol=2e-3)


def test_exact_gradient_and_hessian_agree_with_numderiv():
    ref = load_ref("lc_cholesky")
    beta, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    ll = _LogLike(pd.Series(ref["beta_test"], dtype=float), fixed, pf, inputs)
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
    ll = _LogLike(pd.Series(ref["estimate"], dtype=float), fixed, pf, inputs)
    assert np.all(np.isfinite(ll.scores(ll.base[ll.free].to_numpy())))


def test_posteriors_prediction_and_values_of_time_match_numpy_backend():
    ref = load_ref("lc_cholesky")
    est = pd.Series(ref["estimate"], dtype=float)
    out = {}
    for flag in (False, True):
        _, fixed, pf, inputs = build("lc_cholesky", autodiff=flag)
        model = engine.Estimate(estimate=est, beta_start=est, fixed=fixed,
                                  control=inputs.control, maximum=0.0, LLStart=0.0)
        out[flag] = (
            engine.lc_conditionals(model, pf, inputs)[["X1", "X2"]].to_numpy(),
            engine.prediction(model, pf, inputs, silent=True).drop(columns=["ID", "Observation"]).to_numpy(),
            mtuem.mtuem_values_of_time(model, inputs, pf, nClass=2, Ec="EcI", standard_errors=False, silent=True)["summary"]["mean"].to_numpy(),
        )
    for a, b in zip(out[False], out[True]):
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12)


def test_em_and_multistart_run_with_the_autodiff_backend():
    beta, fixed, pf, inputs = build("lc_cholesky", autodiff=True)
    em = engine.lc_em(beta, fixed, pf, inputs, max_iterations=2, post_em=0, estimate_kwargs={"silent": True})
    assert em.maximum >= em.LLStart
    model = MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w")
    settings = EstimationSettings(backend="jax", max_iterations=100, verbose=False,
                                  multistart=MultiStart(n=2, seed=1, normalization="theta_phi"))
    res = model.fit(database(), model.beta(), Control(indiv_id="PeID"), settings)
    assert res._inputs.autodiff
    assert res.loglik == pytest.approx(load_ref("eq1")["maximum"], abs=1e-3)


def test_degenerate_error_structures_give_minus_infinity():
    _, _, _, inputs = build("eq3_sig", autodiff=True)
    b = inputs.attach(load_ref("eq3_sig")["estimate"])
    s = {"work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"],
         "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
         "times_elasticities": [b.theta_1], "goods_elasticities": [b.phi_1],
         "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168}
    assert np.all(np.isneginf(np.asarray(lj.theta_phi(s | {"sig": [-5.0, 7.0, 40.0]}, "estimate", inputs))))
    assert np.all(np.isneginf(np.asarray(lj.theta_phi(s | {"sig": [5.0, 7.0, 40.0], "rho": [0.9, 0.9, -0.9]}, "estimate", inputs))))
    with pytest.raises(ValueError):
        lj.theta_phi(s | {"cholesky": [1.0, 0.0, 1.0]}, "estimate", inputs)
    assert np.all(np.asarray(lj.theta_phi(s, "validate", inputs)) == 0)
    assert np.all(np.isnan(np.asarray(lj.theta_phi(s, "zero_LL", inputs))))


# ----------------------------------------------------------------------------------------------
# Alpha-beta and additive quadratic: the JAX ports against the numpy likelihoods on maed
# ----------------------------------------------------------------------------------------------
ERRORS = ["profiled", "sig", "sig_rho", "cholesky"]

AB_POINTS = [
    {"alpha": 0.3, "beta": 0.25, "th1": 0.1, "th2": 0.05, "ph1": 0.2},
    {"alpha": 0.2, "beta": 0.3, "th1": 0.15, "th2": 0.04, "ph1": 0.25},
    {"alpha": 0.35, "beta": 0.2, "th1": 0.08, "th2": 0.03, "ph1": 0.15},
]
AQ_POINTS = [
    {"we": 0.0, "ws": 10.0, "theta_1": 800.0, "theta_2": 3000.0, "beta_1": 20.0, "beta_2": 300.0,
     "phi_1": 100.0, "phi_2": 80.0, "eta_1": 1.0, "eta_2": 0.5},
    {"we": 5.0, "ws": 8.0, "theta_1": 600.0, "theta_2": 2000.0, "beta_1": 30.0, "beta_2": 200.0,
     "phi_1": 150.0, "phi_2": 60.0, "eta_1": 2.0, "eta_2": 0.4},
    {"we": -3.0, "ws": 12.0, "theta_1": 1000.0, "theta_2": 2500.0, "beta_1": 25.0, "beta_2": 250.0,
     "phi_1": 120.0, "phi_2": 90.0, "eta_1": 1.5, "eta_2": 0.8},
]
# Standard deviations and correlations of the error terms per equation, close to the residual
# moments at the points above, so that the explicit error structures are positive definite.
SD = {"Tw": 7.0, "Tf1": 10.0, "Tf2": 4.5, "Ef1": 45.0, "Ef2": 40.0}
CORR = {("Tw", "Tf1"): -0.7, ("Tw", "Tf2"): -0.2, ("Tw", "Ef1"): 0.5, ("Tw", "Ef2"): 0.1,
        ("Tf1", "Tf2"): -0.35, ("Tf1", "Ef1"): -0.5, ("Tf1", "Ef2"): -0.1,
        ("Tf2", "Ef1"): 0.05, ("Tf2", "Ef2"): 0.0, ("Ef1", "Ef2"): 0.1}

# Model structures: the equations kept and the settings that select them.
AB_STRUCTURES = {
    "full": (["Tw", "Tf1", "Tf2", "Ef1"], {}),
    "single": (["Tw"], {}),
}
AQ_STRUCTURES = {
    "full": (["Tw", "Tf1", "Tf2", "Ef1", "Ef2"], {}),
    "omitted": (["Tw", "Tf1", "Ef1"], {"omitted_times": ["Tf2"], "omitted_goods": "Ef2"}),
    "single": (["Tw"], {"omitted_times": ["Tf1", "Tf2"], "omitted_goods": ["Ef1", "Ef2"]}),
    "observed_tw": (["Tw", "Tf1", "Tf2", "Ef1", "Ef2"], {"optimal_tw": False}),
}


@pytest.fixture(scope="module")
def maed_inputs():
    return engine.make_inputs(database(), {"b": 0.0}, [], Control(indiv_id="PeID"), autodiff=True)


def error_values(kind: str, eqs: list[str]) -> list[float]:
    """Explicit error parameters of ``kind`` for the equations ``eqs`` (empty when profiled)."""
    M = len(eqs)
    sd = np.array([SD[e] for e in eqs])
    R = np.eye(M)
    for i in range(M):
        for j in range(i + 1, M):
            R[i, j] = R[j, i] = CORR[(eqs[i], eqs[j])]
    if kind == "profiled":
        return []
    if kind == "sig":
        return list(sd)
    if kind == "sig_rho":
        return list(sd) + list(R[np.triu_indices(M, 1)])
    return list(np.linalg.cholesky(R * np.outer(sd, sd))[np.tril_indices(M)])


def error_settings(kind: str, e, M: int) -> dict:
    if kind == "profiled":
        return {}
    if kind == "sig":
        return {"sig": [e[i] for i in range(M)]}
    if kind == "sig_rho":
        return {"sig": [e[i] for i in range(M)], "rho": [e[i] for i in range(M, len(e))]}
    return {"cholesky": [e[i] for i in range(len(e))]}


def ab_settings(p: dict, structure: str, kind: str, e, db) -> dict:
    eqs, _ = AB_STRUCTURES[structure]
    s = {"work_times": ["Tw"], "work_elasticities": {"alpha": p["alpha"], "beta": p["beta"]},
         "Tc": db.Tc.to_numpy(), "Ec": db.EcI.to_numpy(), "w": db.w.to_numpy(), "tau": 168}
    if structure == "full":
        s |= {"free_times": ["Tf1", "Tf2"], "free_goods": ["Ef1"], "goods_cost": {"Ef1": 1},
              "times_elasticities": [p["th1"], p["th2"]], "goods_elasticities": [p["ph1"]]}
    return s | error_settings(kind, e, len(eqs))


def aq_settings(p: dict, structure: str, kind: str, e, db) -> dict:
    eqs, extra = AQ_STRUCTURES[structure]
    s = {"work_times": "Tw", "free_times": ["Tf1", "Tf2"], "free_goods": ["Ef1", "Ef2"], "goods_cost": 1,
         "work_elasticity": p["we"], "work_satiation": p["ws"],
         "times_elasticities": [p["theta_1"], p["theta_2"]], "times_satiations": [p["beta_1"], p["beta_2"]],
         "goods_elasticities": [p["phi_1"], p["phi_2"]], "goods_satiations": [p["eta_1"], p["eta_2"]],
         "Tc": db.Tc.to_numpy(), "Ec": db.EcI.to_numpy(), "w": db.w.to_numpy(), "tau": 168}
    return s | extra | error_settings(kind, e, len(eqs))


FORMS = {
    "alpha_beta": (lk.alpha_beta, lj.alpha_beta, ab_settings, AB_POINTS, AB_STRUCTURES),
    "additive_quadratic": (lk.additive_quadratic, lj.additive_quadratic, aq_settings, AQ_POINTS, AQ_STRUCTURES),
}


def cases_of(form: str) -> list[tuple[str, str]]:
    """(structure, error kind) pairs; a single equation has no correlations."""
    structures = FORMS[form][4]
    return [(st, kind) for st in structures for kind in ERRORS if not (len(structures[st][0]) == 1 and kind == "sig_rho")]


def assert_jax_matches_numpy(form, structure, kind, point, inputs):
    np_fn, jax_fn, make, points, structures = FORMS[form]
    e = error_values(kind, structures[structure][0])
    s = make(points[point], structure, kind, e, inputs.database)
    ll_np = np.asarray(np_fn(s, "estimate", inputs))
    ll_jax = np.asarray(jax_fn(s, "estimate", inputs))
    assert ll_jax.shape == ll_np.shape == (inputs.nObs,)
    assert np.all(np.isfinite(ll_np))
    np.testing.assert_allclose(ll_jax, ll_np, rtol=1e-10, atol=0)
    assert float(ll_jax.sum()) == pytest.approx(float(ll_np.sum()), rel=1e-10)


@pytest.mark.parametrize("point", range(3))
@pytest.mark.parametrize("structure,kind", cases_of("alpha_beta"))
def test_alpha_beta_jax_matches_numpy(structure, kind, point, maed_inputs):
    assert_jax_matches_numpy("alpha_beta", structure, kind, point, maed_inputs)


@pytest.mark.parametrize("point", range(3))
@pytest.mark.parametrize("structure,kind", cases_of("additive_quadratic"))
def test_additive_quadratic_jax_matches_numpy(structure, kind, point, maed_inputs):
    assert_jax_matches_numpy("additive_quadratic", structure, kind, point, maed_inputs)


GRADIENT_CASES = (
    [("alpha_beta", "full", kind) for kind in ERRORS]
    + [("additive_quadratic", st, kind) for st in ("full", "omitted", "observed_tw") for kind in ERRORS]
)


@pytest.mark.parametrize("form,structure,kind", GRADIENT_CASES)
def test_jax_gradient_matches_central_differences(form, structure, kind, maed_inputs):
    np_fn, jax_fn, make, points, structures = FORMS[form]
    names = list(points[0])
    x0 = np.array([points[0][n] for n in names] + error_values(kind, structures[structure][0]))
    k = len(names)

    def total(fn, x):
        s = make(dict(zip(names, x[:k])), structure, kind, x[k:], maed_inputs.database)
        return fn(s, "estimate", maed_inputs).sum()

    g_jax = np.asarray(jax.grad(lambda x: total(jax_fn, x))(jnp.asarray(x0)))
    g_num = numderiv.gradient(lambda x: float(total(np_fn, x)), x0)
    assert np.all(np.isfinite(g_jax))
    np.testing.assert_allclose(g_jax, g_num, rtol=1e-5)


@pytest.mark.parametrize("form", ["alpha_beta", "additive_quadratic"])
def test_degenerate_points_give_minus_infinity_like_numpy(form, maed_inputs):
    np_fn, jax_fn, make, points, structures = FORMS[form]
    eqs = structures["full"][0]
    M = len(eqs)
    bad_sig = [-1.0] + [SD[q] for q in eqs[1:]]
    R = np.eye(M)
    R[0, 1] = R[0, 2] = 0.9
    R[1, 2] = -0.9  # not positive definite
    bad_rho = [SD[q] for q in eqs] + list(R[np.triu_indices(M, 1)])
    for kind, e in (("sig", bad_sig), ("sig_rho", bad_rho)):
        s = make(points[0], "full", kind, e, maed_inputs.database)
        assert np.all(np.isneginf(np.asarray(np_fn(s, "estimate", maed_inputs))))
        assert np.all(np.isneginf(np.asarray(jax_fn(s, "estimate", maed_inputs))))


@pytest.mark.parametrize("form", ["alpha_beta", "additive_quadratic"])
def test_jax_ports_handle_the_other_functionalities(form, maed_inputs):
    np_fn, jax_fn, make, points, structures = FORMS[form]
    s = make(points[0], "full", "sig", error_values("sig", structures["full"][0]), maed_inputs.database)
    assert np.all(np.asarray(jax_fn(s, "validate", maed_inputs)) == 0)
    assert np.all(np.isnan(np.asarray(jax_fn(s, "zero_LL", maed_inputs))))
    pd.testing.assert_frame_equal(jax_fn(s, "prediction", maed_inputs), np_fn(s, "prediction", maed_inputs))
    with pytest.raises(ValueError):
        jax_fn(s, "nonsense", maed_inputs)
    with pytest.raises(ValueError, match="cholesky"):
        jax_fn(s | {"cholesky": [1.0, 0.0, 1.0]}, "estimate", maed_inputs)


def test_alpha_beta_get_covar_is_delegated_to_numpy(maed_inputs):
    s = ab_settings(AB_POINTS[0], "full", "profiled", [], maed_inputs.database)
    got = lj.alpha_beta(s, "get_covar", maed_inputs)
    want = lk.alpha_beta(s, "get_covar", maed_inputs)
    pd.testing.assert_frame_equal(got["covar"], want["covar"])


def test_additive_quadratic_jax_checks_omitted_categories(maed_inputs):
    s = aq_settings(AQ_POINTS[0], "full", "profiled", [], maed_inputs.database)
    with pytest.raises(ValueError, match="omitted_times"):
        lj.additive_quadratic(s | {"omitted_times": "Tf9"}, "estimate", maed_inputs)


def test_ported_forms_and_the_fallback_warning(monkeypatch):
    assert lj.PORTED == frozenset({"theta_phi", "alpha_beta", "additive_quadratic", "lc_profiled"})
    assert issubclass(FallbackWarning, UserWarning)
    monkeypatch.setattr(lj, "_warned", set())
    for form in lj.PORTED:
        assert lj.require_jax_form(form) is True
    with pytest.warns(FallbackWarning, match="numpy_only"):
        assert lj.require_jax_form("numpy_only") is False
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert lj.require_jax_form("numpy_only") is False  # warned once only


def test_cholesky_length_is_checked_before_any_functionality(maed_inputs):
    db = maed_inputs.database
    for form in ("alpha_beta", "additive_quadratic"):
        np_fn, jax_fn, make, points, _ = FORMS[form]
        s = make(points[0], "full", "cholesky", [1.0, 0.0, 1.0], db)
        for functionality in ("validate", "preprocess", "estimate"):
            with pytest.raises(ValueError, match="cholesky"):
                np_fn(s, functionality, maed_inputs)
            with pytest.raises(ValueError, match="cholesky"):
                jax_fn(s, functionality, maed_inputs)
    s = ab_settings(AB_POINTS[0], "full", "cholesky", [1.0, 0.0, 1.0], db) | {"work_elasticities": {"Theta": 1.0, "Phi": 1.0, "thw": 0.0}}
    with pytest.raises(ValueError, match="cholesky"):
        lj.theta_phi(s, "validate", maed_inputs)


# A class whose alpha-beta discriminant is negative in about half of the rows (numpy: NaN residuals there).
AB_NEGATIVE_DISCRIMINANT = {"alpha": 0.55, "beta": 0.0, "th1": 0.1, "th2": 0.05, "ph1": 0.2}


def _non_pd_rho_errors(eqs):
    M = len(eqs)
    R = np.eye(M)
    R[0, 1] = R[0, 2] = 0.9
    R[1, 2] = -0.9
    return [SD[q] for q in eqs] + list(R[np.triu_indices(M, 1)])


DEGENERATE_CLASSES = {
    "non_pd_rho": ("sig_rho", AB_POINTS[1], _non_pd_rho_errors(AB_STRUCTURES["full"][0])),
    "negative_discriminant_sig": ("sig", AB_NEGATIVE_DISCRIMINANT, error_values("sig", AB_STRUCTURES["full"][0])),
    "negative_discriminant_profiled": ("profiled", AB_NEGATIVE_DISCRIMINANT, []),
}


@pytest.mark.parametrize("case", list(DEGENERATE_CLASSES))
def test_lc_gradient_is_finite_when_one_class_is_degenerate(case, maed_inputs):
    kind2, p2, e2 = DEGENERATE_CLASSES[case]
    db, N = maed_inputs.database, maed_inputs.nObs
    names = list(AB_POINTS[0])
    e1 = error_values("sig", AB_STRUCTURES["full"][0])
    x0 = np.array([AB_POINTS[0][n] for n in names] + [p2[n] for n in names] + [0.3] + list(e2))

    def settings(x):
        s1 = ab_settings(dict(zip(names, x[0:5])), "full", "sig", e1, db)
        s2 = ab_settings(dict(zip(names, x[5:10])), "full", kind2, x[11:], db)
        return s1, s2

    def total(x):
        s1, s2 = settings(x)
        pis = lj.class_alloc({"utilities": [jnp.zeros(N), x[10] * jnp.ones(N)]})
        P = {"class_1": lj.alpha_beta(s1, "estimate", maed_inputs), "class_2": lj.alpha_beta(s2, "estimate", maed_inputs)}
        return lj.lc({"inClassProb": P, "classProb": pis}, maed_inputs).sum()

    # The value is what numpy gives, with the degenerate class contributing nothing where it
    # has zero likelihood (numpy's NaN rows are -inf here; a profiled class with an
    # incomplete row is -inf everywhere, see _loglik_from_errors).
    s1, s2 = settings(x0)
    ll1 = np.asarray(lk.alpha_beta(s1, "estimate", maed_inputs))
    ll2_np = np.asarray(lk.alpha_beta(s2, "estimate", maed_inputs))
    ll2_jax = np.asarray(lj.alpha_beta(s2, "estimate", maed_inputs))
    if case == "non_pd_rho":
        assert np.all(np.isneginf(ll2_np)) and np.all(np.isneginf(ll2_jax))
    elif case == "negative_discriminant_sig":
        assert np.any(np.isnan(ll2_np)) and np.any(np.isfinite(ll2_np))
        np.testing.assert_allclose(ll2_jax, np.where(np.isnan(ll2_np), -np.inf, ll2_np), rtol=1e-10)
    else:
        assert np.any(np.isnan(ll2_np)) and np.any(np.isfinite(ll2_np))
        assert np.all(np.isneginf(ll2_jax))
    ll2 = ll2_jax
    lp1, lp2 = -np.log1p(np.exp(0.3)), 0.3 - np.log1p(np.exp(0.3))
    expected = np.sum(np.logaddexp(lp1 + ll1, lp2 + ll2))
    assert np.isfinite(expected)
    assert float(total(jnp.asarray(x0))) == pytest.approx(expected, rel=1e-10)
    g = np.asarray(jax.grad(total)(jnp.asarray(x0)))
    assert np.all(np.isfinite(g))
    # The degenerate class has zero likelihood in the affected rows: no gradient from there.
    if case != "negative_discriminant_sig":
        assert np.all(g[5:10] == 0) and np.all(g[11:] == 0)


def _collinear_inputs(c: float):
    """maed with Tf2 replaced so that its alpha-beta residual (point 0) is a linear
    combination of the Tw and Tf1 residuals plus ``c`` times fixed noise."""
    db = database()
    p = AB_POINTS[0]
    we = {"alpha": p["alpha"], "beta": p["beta"]}
    tw = mtuem.get_tw_albe(we, 168, db.Tc.to_numpy(), db.EcI.to_numpy(), db.w.to_numpy()).ravel()
    ti = mtuem.get_ti_albe([p["th1"], p["th2"]], p["beta"], tw, 168, db.Tc.to_numpy())
    combo = 0.5 * (db.Tw.to_numpy() - tw) - 0.3 * (db.Tf1.to_numpy() - ti[:, 0])
    noise = np.random.default_rng(7).normal(size=len(db))
    db["Tf2"] = ti[:, 1] + combo + c * noise
    return engine.make_inputs(db, {"b": 0.0}, [], Control(indiv_id="PeID"), autodiff=True)


def _corr_eigenvalues(inputs):
    s = ab_settings(AB_POINTS[0], "full", "profiled", [], inputs.database)
    corr = lk.alpha_beta(s, "get_covar", inputs)["corr"].to_numpy()
    return np.linalg.eigvalsh(corr)


def test_profiled_near_collinear_residuals_agree_with_numpy_near_pd():
    """numpy projects the profiled correlation with nearPD and rejects eigenvalues <= 1e-6;
    JAX rejects the same points directly (eigenvalue <= 1e-6, or <= 1e-6 times the largest).
    """
    from scipy.optimize import brentq

    def gap(log_c):
        eig = _corr_eigenvalues(_collinear_inputs(np.exp(log_c)))
        return eig[0] - 1e-6 * (1 + eig[-1]) / 2

    # Exactly collinear; inside numpy's nearPD window 1e-6 < eig <= 1e-6 * largest; well
    # conditioned (smallest eigenvalue about 1e-3).
    c_window = np.exp(brentq(gap, np.log(1e-6), np.log(10.0), xtol=1e-6))
    for c, finite in ((0.0, False), (c_window, False), (0.5, True)):
        inputs = _collinear_inputs(c)
        eig = _corr_eigenvalues(inputs)
        if c == c_window:
            assert 1e-6 < eig[0] < 1e-6 * eig[-1]
        if finite:
            assert eig[0] > 1e-4
        s = ab_settings(AB_POINTS[0], "full", "profiled", [], inputs.database)
        ll_np = np.asarray(lk.alpha_beta(s, "estimate", inputs))
        ll_jax = np.asarray(lj.alpha_beta(s, "estimate", inputs))
        if finite:
            assert np.all(np.isfinite(ll_np))
            np.testing.assert_allclose(ll_jax, ll_np, rtol=1e-10)
        else:
            assert np.all(np.isneginf(ll_np)) and np.all(np.isneginf(ll_jax))
