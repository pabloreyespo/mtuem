"""The activity participation example (examples/activity_participation.py) against the R models
of the article and the fixtures of tests/reference/make_article_reference.R.

The example is a script, not part of the package: it is imported by adding examples/ to
sys.path, which keeps its classes picklable and importable by name (``activity_participation``).
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem
from mtuem import _engine as engine
from mtuem._engine.latent import _loglik_per_indiv
from mtuem.models import LatentClassMTUEM

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
ARTICLE = Path(__file__).parent / "reference" / "article"
sys.path.insert(0, str(EXAMPLES))

import activity_participation as ap  # noqa: E402

RUNS = ["2015_p5", "2015-new_p5", "2022_p5", "2022-bigcities_p4"]
BACKENDS = ["numpy", "jax"]


def fixture(run: str) -> dict:
    return json.loads((ARTICLE / f"{run}.json").read_text(encoding="utf-8"))


def model_of(ref: dict, **kw) -> ap.ActivityParticipationMTUEM:
    participation = [[mix[t] for t in ref["free_times"]] for mix in ref["act_mix"]]
    args = dict(
        work_times=ref["work_times"], free_times=ref["free_times"], participation=participation,
        allocation=ref["allocation"], penalization=ref["penalization"],
    )
    return ap.article_model(**(args | kw))


def compiled(model, data: pd.DataFrame, backend: str, indiv_id: str = "id_persona"):
    prob_fn, lc_pars, autodiff = model._compile(backend)
    beta = model.beta()
    control = mtuem.Control(model_name="ap", indiv_id=indiv_id)
    inputs = engine.make_inputs(data, beta.values(), beta.fixed_names, control, lc_pars, autodiff=autodiff)
    return prob_fn, inputs


def classes_at(model, data, backend, values) -> dict[str, np.ndarray]:
    prob_fn, inputs = compiled(model, data, backend)
    P = prob_fn(values, inputs, "estimate")
    return {k: np.asarray(v) for k, v in P.items()}


def unpenalized(model) -> LatentClassMTUEM:
    """The same latent class model without the participation penalty."""
    fields = model.model_dump(exclude={"penalization", "likelihood_floor"})
    return LatentClassMTUEM(**fields)


# ----------------------------------------------------------------------------------------------
# The penalty and the floor
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("mixes", [None, [[True, False, False], [True, True, True]]], ids=["article", "two_absent"])
def test_penalty_matches_r_formula(backend, mixes):
    # R: P_s = density_s * exp(-penalization * rowSums(hours of the absent activities)); in logs
    # the class log-likelihood loses penalization times those hours, row by row. The second case
    # puts two absent activities in one class, which the article fixtures no longer exercise.
    ref = fixture("2015_p5")
    data = pd.read_parquet(ARTICLE / "2015_p5.parquet")
    if mixes is None:
        model = model_of(ref)
        mixes = [[mix[t] for t in ref["free_times"]] for mix in ref["act_mix"]]
    else:
        model = model_of(ref, participation=mixes)
    values = ref["r_estimate"]
    P = classes_at(model, data, backend, values)
    B = classes_at(unpenalized(model), data, backend, values)
    for s, mix in enumerate(mixes, start=1):
        absent = [t for t, takes_part in zip(ref["free_times"], mix) if not takes_part]
        hours = data[absent].sum(axis=1).to_numpy(float) if absent else np.zeros(len(data))
        assert model.absent_times(s) == absent
        np.testing.assert_allclose(P[f"class_{s}"], B[f"class_{s}"] - ref["penalization"] * hours, rtol=1e-12)
    assert np.any(data[model.absent_times(1)].sum(axis=1) > 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_floor_applies_in_logs(backend):
    # R: pmax(P_s, floor) on densities is max(log P_s, log floor) in logs, in every class.
    ref = fixture("2015_p5")
    data = pd.read_parquet(ARTICLE / "2015_p5.parquet")
    values = ref["r_estimate"]
    free = classes_at(model_of(ref), data, backend, values)
    floor = math.exp(float(np.median(free["class_1"])))
    floored = classes_at(model_of(ref, likelihood_floor=floor), data, backend, values)
    for s in (1, 2):
        np.testing.assert_allclose(floored[f"class_{s}"], np.maximum(free[f"class_{s}"], math.log(floor)), rtol=1e-12)
    clamped = free["class_1"] < math.log(floor)
    assert clamped.sum() >= len(data) // 3
    assert np.all(floored["class_1"][clamped] == math.log(floor))


def test_floor_and_penalization_must_not_be_negative():
    ref = fixture("2015_p5")
    with pytest.raises(ValueError, match="penalization"):
        model_of(ref, penalization=-1.0)
    with pytest.raises(ValueError, match="likelihood_floor"):
        model_of(ref, likelihood_floor=-1e-300)


# ----------------------------------------------------------------------------------------------
# Start values and the fixture models
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("run", RUNS)
def test_article_beta_matches_r_start_values(run):
    ref = fixture(run)
    beta = ap.article_beta(model_of(ref))
    assert beta.names == list(ref["beta_start"])
    assert beta.values() == pytest.approx(ref["beta_start"], abs=0)
    assert set(beta.fixed_names) == set(ref["fixed"])


def test_article_beta_keeps_the_locks():
    beta = ap.article_beta(model_of(fixture("2015_p5")))
    assert {"th2_1", "chk31_1", "chk32_1", "chk33_1", "chk43_1"} == set(beta.locked)
    with pytest.raises(ValueError, match="locked"):
        beta.free("th2_1")


@pytest.mark.parametrize("run", RUNS)
def test_from_fixture_rebuilds_the_r_model(run):
    ref = fixture(run)
    model, beta, data, control = ap.from_fixture(ARTICLE / f"{run}.json")
    assert isinstance(model, ap.ActivityParticipationMTUEM)
    assert model.penalization == ref["penalization"] and model.likelihood_floor == ref["class_likelihood_floor"]
    assert model.n_classes == ref["n_classes"] and model.allocation == ref["allocation"]
    assert beta.values() == ref["beta_start"] and set(beta.fixed_names) == set(ref["fixed"])
    assert beta.locked == model.beta().locked
    assert len(data) == ref["r_nobs"] and set(model.columns()) <= set(data.columns)
    assert control.indiv_id == ref["indiv_id"] and control.model_name == ref["model_name"]


def test_from_fixture_rejects_a_fixture_that_frees_a_locked_parameter(tmp_path):
    ref = fixture("2015_p5")
    ref["fixed"].remove("th2_1")
    path = tmp_path / "2015_p5.json"
    path.write_text(json.dumps(ref), encoding="utf-8")
    (tmp_path / "2015_p5.parquet").write_bytes((ARTICLE / "2015_p5.parquet").read_bytes())
    with pytest.raises(ValueError, match="th2_1"):
        ap.from_fixture(path)


# ----------------------------------------------------------------------------------------------
# The Python model is the R model: the log-likelihood at the R estimates
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("run", RUNS)
def test_loglik_at_r_estimate_matches_r(run, backend):
    ref = fixture(run)
    model, beta, data, control = ap.from_fixture(ARTICLE / f"{run}.json")
    prob_fn, inputs = compiled(model, data, backend, control.indiv_id)
    P = prob_fn(ref["r_estimate"], inputs, "estimate")
    ll = float(np.sum(_loglik_per_indiv(P["model"], inputs)))
    assert ll == pytest.approx(ref["r_loglik"], rel=1e-10)
