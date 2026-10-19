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
