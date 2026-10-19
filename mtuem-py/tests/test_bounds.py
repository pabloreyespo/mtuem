"""Parameter bounds: the adaptive log-barrier of R's constrOptim in BFGS and in the EM M step."""

import json
from pathlib import Path

import numpy as np
import pytest

from cases import build
from mtuem import _engine as engine
from mtuem._engine.optimize import _maximise, maximise_bounded
from test_engine import small_lc

REF = Path(__file__).parent / "reference"


def load_ref(case: str) -> dict:
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def eq1_free():
    beta, fixed, pf, inputs = build("eq1")
    return engine.estimate(beta, fixed, pf, inputs, silent=True)


def test_unbinding_bounds_match_unbounded(eq1_free):
    beta, fixed, pf, inputs = build("eq1")
    # Theta is fixed: its bound is ignored even though the start lies outside it.
    bounds = {"Phi": (0.0, 10.0), "thw": (-5.0, None), "Theta": (5.0, 6.0)}
    model = engine.estimate(beta, fixed, pf, inputs, bounds=bounds, silent=True)
    assert model.code == 0
    np.testing.assert_allclose(model.estimate, eq1_free.estimate, rtol=0, atol=1e-8)
    assert model.maximum == pytest.approx(eq1_free.maximum, abs=1e-8)


@pytest.mark.parametrize("autodiff", [False, True])
def test_binding_bound_matches_r(autodiff):
    ref = load_ref("eq1_bounded")
    lower = ref["bounds"]["thw"]["lower"]
    _, fixed, pf, inputs = build("eq1", autodiff=autodiff)
    model = engine.estimate(ref["beta_start"], fixed, pf, inputs, bounds={"thw": (lower, None)}, silent=True)
    assert model.code == 0
    assert model.estimate["thw"] >= lower
    # The reference is the constrained maximum (thw fixed at the binding bound in R).
    names = list(ref["estimate"])
    np.testing.assert_allclose(model.estimate[names], [ref["estimate"][n] for n in names], rtol=1e-4)
    assert model.maximum == pytest.approx(ref["maximum"], abs=1e-3)
    # Standard errors come from the unpenalised log-likelihood Hessian, as in maxLik and Apollo.
    free = [n for n in names if n not in fixed]
    rtol = 5e-3 if autodiff else 2e-3
    np.testing.assert_allclose(model.se[free], [ref["se"][n] for n in free], rtol=rtol)
    np.testing.assert_allclose(model.robse[free], [ref["robse"][n] for n in free], rtol=rtol)
    # maxLik's own barrier run (Apollo bfgs with the constraint) stalls short of that maximum:
    # it reaches the bound but not the best Phi.
    barrier = ref["barrier"]
    assert model.estimate["thw"] == pytest.approx(barrier["estimate"]["thw"], rel=1e-4)
    assert model.maximum >= barrier["maximum"]


def test_start_on_bound_raises():
    beta, fixed, pf, inputs = build("eq1")
    with pytest.raises(ValueError, match="thw"):
        engine.estimate(beta, fixed, pf, inputs, bounds={"thw": (0.0, None)}, silent=True)
    with pytest.raises(ValueError, match="Phi"):
        engine.estimate(beta, fixed, pf, inputs, bounds={"Phi": (None, 0.5)}, silent=True)
    with pytest.raises(ValueError, match="nope"):
        engine.estimate(beta, fixed, pf, inputs, bounds={"nope": (0.0, 1.0)}, silent=True)
    with pytest.raises(ValueError, match="Phi"):
        engine.estimate(beta, fixed, pf, inputs, bounds={"Phi": (2.0, 0.0)}, silent=True)


@pytest.mark.slow
def test_em_respects_bounds():
    beta, fixed, pf, inputs = small_lc(n=150)
    # From Phi_1 = 0.5 the unbounded EM takes Phi_1 to about 0.27; a lower bound of 0.45 binds.
    lower = 0.45
    beta = beta | {"Phi_1": 0.5}
    bounds = {"Phi_1": (lower, None)}
    em_only = engine.lc_em(beta, fixed, pf, inputs, max_iterations=2, post_em=0, bounds=bounds,
                           estimate_kwargs={"silent": True})
    assert em_only.estimate["Phi_1"] >= lower
    assert em_only.estimate["Phi_1"] == pytest.approx(lower, abs=1e-6)
    assert em_only.maximum >= em_only.LLStart
    # The bounds are passed on to the estimation that follows the EM iterations.
    full = engine.lc_em(em_only.estimate.to_dict(), fixed, pf, inputs, max_iterations=1, bounds=bounds,
                        estimate_kwargs={"silent": True, "max_iterations": 3, "hessian": "none"})
    assert full.estimate["Phi_1"] >= lower


def test_em_rejects_a_start_outside_the_bounds():
    beta, fixed, pf, inputs = small_lc(n=50)
    with pytest.raises(ValueError, match="Phi_1"):
        engine.lc_em(beta, fixed, pf, inputs, bounds={"Phi_1": (0.6, None)}, estimate_kwargs={"silent": True})


def test_maximise_bounded_on_a_quadratic():
    # f(x) = -(x0 - 2)^2 - (x1 + 1)^2, maximum at (2, -1).
    def f(x):
        return np.array([-(x[0] - 2.0) ** 2, -(x[1] + 1.0) ** 2])

    res = maximise_bounded(f, np.array([0.0, 0.0]), np.array([-np.inf, -0.5]), np.array([1.0, np.inf]))
    assert res.code == 0
    np.testing.assert_allclose(res.x, [1.0, -0.5], atol=1e-5)
    assert res.x[0] < 1.0 and res.x[1] > -0.5
    assert res.f == pytest.approx(-1.25, abs=1e-4)

    # No finite bound: the plain optimiser.
    free = maximise_bounded(f, np.zeros(2), np.full(2, -np.inf), np.full(2, np.inf))
    plain = _maximise(f, np.zeros(2))
    np.testing.assert_array_equal(free.x, plain.x)

    with pytest.raises(ValueError):
        maximise_bounded(f, np.array([1.0, 0.0]), np.array([-np.inf, -0.5]), np.array([1.0, np.inf]))


def coupled(x):
    # f = -(x0 - 2)^2 - 10 (x0 x1 - 1)^2 - (x1 + 1)^2; with x0 < 1 and x1 > -0.5 the maximum is
    # at x0 = 1, x1 = 9/11, f = -51/11 (the x0 bound binds).
    return np.array([-(x[0] - 2.0) ** 2, -10.0 * (x[0] * x[1] - 1.0) ** 2, -(x[1] + 1.0) ** 2])


COUPLED_BOUNDS = (np.array([-np.inf, -0.5]), np.array([1.0, np.inf]))


@pytest.mark.parametrize("level", [0.0, 1000.0])
def test_barrier_converges_on_a_coupled_problem(level):
    # A constant shifts f to log-likelihood scale; the binding parameter must still land on its bound.
    def shifted(x):
        v = coupled(x)
        v[0] -= level
        return v

    res = maximise_bounded(shifted, np.zeros(2), *COUPLED_BOUNDS)
    assert res.code == 0 and res.message == "successful convergence"
    assert res.iterations <= 500
    assert 0.0 < 1.0 - res.x[0] < 1e-7
    # x1 is unbounded: its accuracy is that of BFGS's relative tolerance on f (~1e-10 |f|).
    assert res.x[1] == pytest.approx(9 / 11, abs=1e-6 if level == 0 else 1e-4)
    assert res.f == pytest.approx(-51 / 11 - level, abs=1e-7)


@pytest.mark.parametrize("sc", [1e4, 1e6])
def test_barrier_converges_on_a_badly_scaled_problem(sc):
    # f = -(x0 - 2)^2 - sc (x1 - 3 / sqrt(sc))^2 with x0 < 1: curvatures 2 and 2 sc.
    def f(x):
        with np.errstate(over="ignore"):  # trial points far along x1 overflow to -inf
            return np.array([-(x[0] - 2.0) ** 2, -sc * (x[1] - 3.0 / np.sqrt(sc)) ** 2])

    res = maximise_bounded(f, np.zeros(2), np.full(2, -np.inf), np.array([1.0, np.inf]))
    assert res.code == 0
    assert 0.0 < 1.0 - res.x[0] < 1e-7
    assert res.x[1] == pytest.approx(3.0 / np.sqrt(sc), rel=1e-6)


def test_exhausted_budget_inside_the_barrier_reports_the_iteration_limit():
    # The free run takes 16 iterations, so the barrier runs out of the 33 available.
    res = maximise_bounded(coupled, np.zeros(2), *COUPLED_BOUNDS, max_iterations=33)
    assert res.code == 10 and res.message == "iteration limit reached without convergence"
    assert res.iterations <= 33
    # A budget used up by the free run leaves none for the barrier: its point (x0 = 2.26),
    # brought inside the bounds, is returned.
    res = maximise_bounded(coupled, np.zeros(2), *COUPLED_BOUNDS, max_iterations=16)
    assert res.code == 10 and res.iterations <= 16
    assert 0.0 < 1.0 - res.x[0] < 1e-11


def test_result_is_kept_strictly_inside_the_bounds():
    res = maximise_bounded(coupled, np.array([1.0 - 1e-15, 0.0]), *COUPLED_BOUNDS)
    assert res.code == 0
    assert 1.0 - res.x[0] > 5e-13  # 1e-12 inside, up to rounding


@pytest.mark.parametrize("level, d, c", [(1e3, 1e-6, 1e6), (1e3, 1e-5, 1e4), (1e5, 1e-3, 1e2)])
def test_interior_parameter_next_to_a_bound_stays_interior(level, d, c):
    # x2 < 1 binds; at x2 = 1 the maximum in x0 is d, just inside its bound x0 > 0.
    def f(x):
        w = x[0] - d - 0.3 * (x[2] - 1.0)
        return np.array([-level - c * w * w, -(x[1] - 1.0) ** 2, -(x[2] - 2.0) ** 2])

    with np.errstate(over="ignore", invalid="ignore"):
        res = maximise_bounded(f, np.array([0.5, 0.0, 0.0]), np.array([0.0, -np.inf, -np.inf]),
                               np.array([np.inf, np.inf, 1.0]))
    assert res.code == 0
    assert abs(res.x[0] - d) <= 1e-2 * d
    assert 0.0 < 1.0 - res.x[2] < 1e-7
    assert res.f == pytest.approx(-level - 1.0, abs=1e-9 * level + 1e-9)
