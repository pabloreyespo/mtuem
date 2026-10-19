"""JAX profiled latent class likelihood (likelihood_jax.lc_profiled) against the numpy code and R."""

import json
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import build
from mtuem import _engine as engine
from mtuem import likelihood_jax as lj
from mtuem import numderiv
from mtuem._engine.latent import _loglik_per_indiv
from mtuem._engine.loglike import _LogLike

REF = Path(__file__).parent / "reference"
POINTS = [("beta_start", "ll_start"), ("beta_test", "ll_test"), ("estimate", "ll_at_estimate")]


# The maximum the optimiser reaches from the README start (log-likelihood about -7960.38,
# far above R's -8966.61, where R stopped with false convergence and no Hessian).
OPTIMUM = {
    "Theta_1": 1.0, "Phi_1": 0.195088, "thw_1": 0.4737, "th1_1": 0.767191, "ph1_1": 0.091912,
    "asc_1": 0.0, "x_female_1": 0.0, "x_older45_1": 0.0,
    "Theta_2": 1.0, "Phi_2": 0.411142, "thw_2": -0.5839, "th1_2": 0.710404, "ph1_2": 0.227713,
    "asc_2": 0.347933, "x_female_2": 1.347919, "x_older45_2": -0.850222,
}


def load_ref():
    return json.loads((REF / "lc_profiled_ref.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def models():
    _, fixed, pf_np, inputs_np = build("lc_profiled")
    _, _, pf_jax, inputs_jax = build("lc_profiled", autodiff=True)
    return fixed, pf_np, inputs_np, pf_jax, inputs_jax


def _settings(models, beta):
    """lc_settings of the profiled model at ``beta`` (JAX class probabilities)."""
    _, _, _, pf_jax, inputs = models
    captured = {}
    original = lj.lc_profiled

    def spy(lc_settings, functionality="estimate", inputs=None):
        captured["s"] = lc_settings
        return original(lc_settings, functionality, inputs)

    lj.lc_profiled = spy
    try:
        pf_jax(beta, inputs, "estimate")
    finally:
        lj.lc_profiled = original
    return captured["s"]


@pytest.mark.parametrize("point", POINTS, ids=[p[0] for p in POINTS])
def test_lc_profiled_jax_value_matches_numpy(point, models):
    ref = load_ref()
    _, pf_np, inputs_np, pf_jax, inputs_jax = models
    ll_np = np.asarray(pf_np(ref[point[0]], inputs_np, "estimate")["model"], dtype=float)
    ll_jax = np.asarray(pf_jax(ref[point[0]], inputs_jax, "estimate")["model"])
    assert ll_jax.shape == ll_np.shape
    np.testing.assert_allclose(ll_jax, ll_np, rtol=1e-10, atol=0)
    assert float(ll_jax.sum()) == pytest.approx(float(ll_np.sum()), rel=1e-10)
    assert float(ll_jax.sum()) == pytest.approx(ref[point[1]], rel=1e-10, abs=1e-8)


def test_lc_profiled_jax_value_matches_numpy_under_jit(models):
    ref = load_ref()
    fixed, pf_np, inputs_np, pf_jax, inputs_jax = models
    ll = _LogLike(pd.Series(ref["beta_test"], dtype=float), fixed, pf_jax, inputs_jax)
    x = ll.base[ll.free].to_numpy()
    want = _loglik_per_indiv(pf_np(ref["beta_test"], inputs_np, "estimate")["model"], inputs_np)
    np.testing.assert_allclose(ll.per_indiv(x), want, rtol=1e-10, atol=0)


def test_the_identical_class_start_is_an_isolated_point_in_both_backends(models):
    # At the README start both classes are identical, a fixed point of the covariance updates
    # that any parameter step breaks: the classes separate and the log-likelihood jumps by
    # about 39. Central differences there measure the separated branch. The R estimate is
    # next to such a point (the classes differ by about 1e-5 in all but the constant), and
    # steps of 1e-8 in th1_1 or ph1_1 already jump branch. The derivative tests use the test
    # point and a point halfway to the optimum; the values agree everywhere.
    ref = load_ref()
    _, pf_np, inputs_np, pf_jax, inputs_jax = models
    for h in (0.0, 1e-6, -1e-6, 1e-3):
        beta = dict(ref["beta_start"], Phi_1=ref["beta_start"]["Phi_1"] + h)
        ll_np = np.asarray(pf_np(beta, inputs_np, "estimate")["model"], dtype=float)
        ll_jax = np.asarray(pf_jax(beta, inputs_jax, "estimate")["model"])
        np.testing.assert_allclose(ll_jax, ll_np, rtol=1e-9, atol=0)
        if h:
            assert ll_np.sum() > ref["ll_start"] + 30


@pytest.mark.parametrize("point", ["beta_test", "halfway"])
def test_lc_profiled_jax_gradient_matches_central_differences(point, models):
    ref = load_ref()
    fixed, pf_np, inputs_np, pf_jax, inputs_jax = models
    if point == "halfway":
        beta = (pd.Series(ref["beta_test"], dtype=float) + pd.Series(OPTIMUM, dtype=float)) / 2
    else:
        beta = pd.Series(ref[point], dtype=float)
    ll_jax = _LogLike(beta, fixed, pf_jax, inputs_jax)
    ll_np = _LogLike(beta, fixed, pf_np, inputs_np)
    x = ll_jax.base[ll_jax.free].to_numpy()
    g_exact = ll_jax.scores(x).sum(axis=0)
    g_num = numderiv.gradient(ll_np.total, x)
    assert np.all(np.isfinite(g_exact))
    np.testing.assert_allclose(g_exact, g_num, rtol=1e-5, atol=1e-6)


def test_lc_profiled_jax_scores_match_numerical_scores(models):
    # The covariance depends on every row, so each individual's score includes the response
    # of the profiled covariance to the parameters, not only the direct effect.
    ref = load_ref()
    fixed, pf_np, inputs_np, pf_jax, inputs_jax = models
    beta = pd.Series(ref["beta_test"], dtype=float)
    ll_jax = _LogLike(beta, fixed, pf_jax, inputs_jax)
    ll_np = _LogLike(beta, fixed, pf_np, inputs_np)
    x = ll_jax.base[ll_jax.free].to_numpy()
    S_exact = ll_jax.scores(x)
    S_num = ll_np.scores(x)
    scale = np.abs(S_num).max(axis=0)
    np.testing.assert_allclose(S_exact / scale, S_num / scale, rtol=0, atol=1e-6)


def _se(H):
    return np.sqrt(np.diag(np.linalg.inv(-H)))


@pytest.mark.slow
def test_lc_profiled_jax_hessian_gives_the_numpy_standard_errors_at_the_optimum(models):
    fixed, pf_np, inputs_np, pf_jax, inputs_jax = models
    beta = pd.Series(OPTIMUM, dtype=float)
    ll_jax = _LogLike(beta, fixed, pf_jax, inputs_jax)
    ll_np = _LogLike(beta, fixed, pf_np, inputs_np)
    x = ll_jax.base[ll_jax.free].to_numpy()
    assert ll_jax.total(x) == pytest.approx(-7960.3818, abs=1e-3)
    H_exact = ll_jax.hessian(x)
    H_num = ll_np.hessian(x, mode="numeric")
    assert np.all(np.isfinite(H_exact))
    assert np.all(np.linalg.eigvalsh((H_exact + H_exact.T) / 2) < 0)
    np.testing.assert_allclose(H_exact, H_exact.T, rtol=1e-6, atol=1e-12 * np.abs(H_exact).max())
    np.testing.assert_allclose(_se(H_exact), _se(H_num), rtol=2e-3)


@pytest.mark.slow
def test_lc_profiled_jax_hessian_matches_numerical_hessian_at_the_test_point(models):
    ref = load_ref()
    fixed, pf_np, inputs_np, pf_jax, inputs_jax = models
    beta = pd.Series(ref["beta_test"], dtype=float)
    ll_jax = _LogLike(beta, fixed, pf_jax, inputs_jax)
    x = ll_jax.base[ll_jax.free].to_numpy()
    H_exact = ll_jax.hessian(x)
    # Numerical Jacobian of the exact gradient: cheap and independent of the Hessian code path.
    H_num = numderiv.jacobian(lambda z: ll_jax.scores(z).sum(axis=0), x)
    np.testing.assert_allclose(H_exact, H_num, rtol=1e-6, atol=1e-6 * np.abs(H_num).max())


@pytest.mark.slow
def test_lc_profiled_jax_estimation_reaches_the_r_likelihood(models):
    ref = load_ref()
    fixed, _, _, pf_jax, inputs_jax = models
    model = engine.estimate(ref["beta_start"], fixed, pf_jax, inputs_jax, silent=True, hessian="none")
    assert model.maximum >= ref["maximum"] - 1e-3


def test_lc_profiled_jax_follows_the_em_settings_like_numpy(models):
    ref = load_ref()
    _, pf_np, inputs_np, _, inputs_jax = models
    s_jax = _settings(models, ref["beta_test"])
    s_np = {"class_settings": [lj._to_numpy_settings(c) for c in s_jax["class_settings"]],
            "classProb": [np.asarray(p) for p in s_jax["classProb"]]}
    for em in ({"maxIterations": 1, "tol": 0}, {"maxIterations": 7}, {"tol": 1e-3}):
        got = np.asarray(lj.lc_profiled(dict(s_jax, em_settings=em), "estimate", inputs_jax))
        want = mtuem.lc_profiled.lc_profiled(dict(s_np, em_settings=em), "estimate", inputs_np)
        np.testing.assert_allclose(got, want, rtol=1e-10, atol=0)


@pytest.mark.slow
def test_lc_profiled_jax_degenerate_points_give_minus_infinity_with_finite_gradients(models):
    ref = load_ref()
    _, _, inputs_np, _, inputs_jax = models
    s_jax = _settings(models, ref["beta_test"])
    N = inputs_jax.nObs
    # A minimum class weight above the sample size: numpy rejects the point.
    too_heavy = dict(s_jax, em_settings={"minWeight": N + 1})
    s_np = {"class_settings": [lj._to_numpy_settings(c) for c in s_jax["class_settings"]],
            "classProb": [np.asarray(p) for p in s_jax["classProb"]], "em_settings": {"minWeight": N + 1}}
    assert np.all(np.isneginf(mtuem.lc_profiled.lc_profiled(s_np, "estimate", inputs_np)))
    assert np.all(np.isneginf(np.asarray(lj.lc_profiled(too_heavy, "estimate", inputs_jax))))

    # A negative discriminant in the optimal work time of class 2 for a few rows: numpy has
    # NaN residuals there and rejects the point.
    def negative_discriminant(phi, base=s_jax):
        c1, c2 = base["class_settings"]
        c2 = dict(c2, work_elasticities=dict(c2["work_elasticities"], Phi=phi, thw=3.0))
        return dict(base, class_settings=[c1, c2])

    def heavy(phi):
        return dict(too_heavy, class_settings=[
            dict(c, work_elasticities=dict(c["work_elasticities"], Phi=phi)) for c in s_jax["class_settings"]])

    s_bad = negative_discriminant(-0.1)
    s_bad_np = dict(s_np, em_settings={},
                    class_settings=[lj._to_numpy_settings(c) for c in s_bad["class_settings"]])
    assert np.all(np.isneginf(mtuem.lc_profiled.lc_profiled(s_bad_np, "estimate", inputs_np)))
    assert np.all(np.isneginf(np.asarray(lj.lc_profiled(s_bad, "estimate", inputs_jax))))

    for settings_of, phi in ((negative_discriminant, -0.1), (heavy, 1.1)):
        def total(p, settings_of=settings_of):
            return jnp.sum(lj.lc_profiled(settings_of(p), "estimate", inputs_jax))

        assert np.isneginf(float(total(phi)))
        assert np.isfinite(float(jax.grad(total)(phi)))
        assert np.all(np.isfinite(np.asarray(jax.hessian(total)(phi))))


def test_lc_profiled_jax_answers_the_other_functionalities(models, capsys):
    ref = load_ref()
    _, pf_np, inputs_np, _, inputs_jax = models
    s = _settings(models, ref["beta_test"])
    N = inputs_jax.nObs
    assert np.all(np.asarray(lj.lc_profiled(s, "validate", inputs_jax)) == 0)
    assert np.all(np.isnan(np.asarray(lj.lc_profiled(s, "zero_LL", inputs_jax))))
    assert lj.lc_profiled(dict(s, componentName="lc"), "preprocess", inputs_jax)["componentName"] == "lc"
    with pytest.raises(ValueError, match="Invalid value of argument"):
        lj.lc_profiled(s, "nonsense", inputs_jax)
    want = pf_np(ref["beta_test"], inputs_np, "prediction")["model"]
    pd.testing.assert_frame_equal(lj.lc_profiled(s, "prediction", inputs_jax), want)
    for functionality in ("output", "conditionals"):
        got = np.asarray(lj.lc_profiled(s, functionality, inputs_jax))
        np.testing.assert_allclose(got, pf_np(ref["beta_test"], inputs_np, "estimate")["model"], rtol=1e-10)
    capsys.readouterr()
    lj.lc_profiled(s, "report", inputs_jax)
    assert "Class 2 VoL" in capsys.readouterr().out
    assert len(np.asarray(lj.lc_profiled(s, "estimate", inputs_jax))) == N


def test_profiled_sigma_after_a_jax_evaluation_matches_numpy(models):
    ref = load_ref()
    _, pf_np, inputs_np, pf_jax, inputs_jax = models
    pf_np(ref["estimate"], inputs_np, "estimate")
    want = mtuem.lc_profiled.get_profiled_sigma()
    pf_jax(ref["estimate"], inputs_jax, "estimate")
    got = mtuem.lc_profiled.get_profiled_sigma()
    for s in range(2):
        np.testing.assert_allclose(got["Sigma"][s], want["Sigma"][s], rtol=1e-10)
        np.testing.assert_allclose(got["Sigma"][s], np.array(ref["profiled"]["Sigma"][s]), rtol=1e-8)
    np.testing.assert_allclose(got["posterior"], want["posterior"], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(got["shares"], want["shares"], rtol=1e-10)


def test_lc_profiled_is_ported():
    assert "lc_profiled" in lj.PORTED
    assert lj.require_jax_form("lc_profiled") is True
