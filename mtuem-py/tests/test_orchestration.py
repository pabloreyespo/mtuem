"""Estimation orchestration of fit: EM with fallback, and multistart."""

import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import database, declarative, latent_class
from mtuem import _engine as engine
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import LatentClassMTUEM
from mtuem.multistart import _candidates

REF = Path(__file__).parent / "reference"
CONTROL = Control(model_name="orchestration", indiv_id="PeID")


def load_ref(case: str) -> dict:
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


def subset(n: int) -> pd.DataFrame:
    return database().iloc[:n].reset_index(drop=True)


def lc_start(model):
    return model.beta().set(**load_ref("lc_cholesky")["beta_start"])


# A box where every eq1 candidate has a finite starting log-likelihood.
NARROW = {"lower": {"Phi": 0.9, "thw": -0.1}, "upper": {"Phi": 1.1, "thw": 0.1}}


def multistart(**kw) -> EstimationSettings:
    """Quiet NumPy settings with a multistart."""
    return EstimationSettings(backend="numpy", verbose=False, multistart=MultiStart(**({"n": 3, "seed": 1} | kw)))


# ----------------------------------------------------------------------------------------------
# EM and fallback
# ----------------------------------------------------------------------------------------------
def test_em_then_bfgs_matches_direct_bfgs_optimum():
    model = latent_class("lc_cholesky")
    em = model.fit(database(), lc_start(model), CONTROL, EstimationSettings(method="em", verbose=False))
    assert em.converged and em.notes == []
    assert em._estimate.routine.startswith("EM algorithm")
    # EM hands over to BFGS: BFGS started at its estimates stays at the same optimum.
    direct = model.fit(database(), em.beta.model_copy(deep=True).set(**em.estimates.to_dict()), CONTROL,
                       EstimationSettings(verbose=False))
    assert em.loglik == pytest.approx(direct.loglik, abs=1e-3)
    # From the R start values EM reaches at least the R optimum.
    assert em.loglik >= load_ref("lc_cholesky")["maximum"] - 1e-3


@pytest.fixture
def small_lc():
    """A small two-class, one-equation model (quick to estimate) with the classes started apart."""
    model = LatentClassMTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", covariance="cholesky", n_classes=2,
                             allocation=["female"])
    return model, subset(150), model.beta().set(Phi_2=0.8, thw_2=0.1)


def test_fallback_runs_when_em_fails(small_lc, monkeypatch):
    model, data, beta = small_lc
    bfgs = model.fit(data, beta, CONTROL, EstimationSettings(backend="numpy", hessian="none", verbose=False))

    def failing_em(*args, **kwargs):
        raise FloatingPointError("boom")

    monkeypatch.setattr(engine, "lc_em", failing_em)
    result = model.fit(data, beta, CONTROL, EstimationSettings(method="em", backend="numpy", hessian="none",
                                                               verbose=False))
    assert result.notes == ["EM failed (boom); estimated by BFGS from the start values"]
    assert result.loglik == bfgs.loglik
    assert result.beta.values() == beta.values()


@pytest.mark.parametrize("error", [ValueError("bad"), np.linalg.LinAlgError("singular")])
def test_fallback_catches_value_and_linalg_errors(small_lc, monkeypatch, error):
    model, data, beta = small_lc

    def failing_em(*args, **kwargs):
        raise error

    monkeypatch.setattr(engine, "lc_em", failing_em)
    result = model.fit(data, beta, CONTROL, EstimationSettings(method="em", backend="numpy", hessian="none",
                                                               verbose=False))
    assert result.notes == [f"EM failed ({error}); estimated by BFGS from the start values"]


def _not_converged(real):
    def em(*args, **kwargs):
        est = real(*args, **kwargs)
        return dataclasses.replace(est, code=1, message="iteration limit", successfulEstimation=False)

    return em


def test_fallback_runs_when_em_does_not_converge(small_lc, monkeypatch):
    model, data, beta = small_lc
    monkeypatch.setattr(engine, "lc_em", _not_converged(engine.lc_em))
    settings = EstimationSettings(method="em", backend="numpy", hessian="none", em_max_iterations=1,
                                  verbose=False)
    result = model.fit(data, beta, CONTROL, settings)
    assert result.notes == [
        "EM failed (not converged: iteration limit); estimated by BFGS from the start values"
    ]
    assert result._estimate.routine == "bfgs"


def test_without_fallback_em_errors_propagate(small_lc, monkeypatch):
    model, data, beta = small_lc

    def failing_em(*args, **kwargs):
        raise FloatingPointError("boom")

    monkeypatch.setattr(engine, "lc_em", failing_em)
    with pytest.raises(FloatingPointError, match="boom"):
        model.fit(data, beta, CONTROL, EstimationSettings(method="em", fallback=None, verbose=False))


def test_without_fallback_a_non_converged_em_result_is_returned_with_a_note(small_lc, monkeypatch):
    model, data, beta = small_lc
    monkeypatch.setattr(engine, "lc_em", _not_converged(engine.lc_em))
    settings = EstimationSettings(method="em", fallback=None, backend="numpy", hessian="none",
                                  em_max_iterations=1, verbose=False)
    result = model.fit(data, beta, CONTROL, settings)
    assert not result.converged
    assert result.notes == ["EM did not converge (iteration limit); no fallback"]


def test_em_on_single_class_raises():
    model = declarative("eq1")
    with pytest.raises(ValueError, match="method='em' needs a latent class model"):
        model.fit(subset(50), model.beta(), CONTROL, EstimationSettings(method="em", verbose=False))
    settings = EstimationSettings(verbose=False, multistart=MultiStart(n=2, first_em=True))
    with pytest.raises(ValueError, match="first_em needs a latent class model"):
        model.fit(subset(50), model.beta(), CONTROL, settings)


def test_em_on_a_profiled_model_raises():
    model = latent_class("lc_profiled")
    with pytest.raises(ValueError, match="profiled"):
        model.fit(subset(50), model.beta(), CONTROL, EstimationSettings(method="em", verbose=False))


def test_notes_are_saved_loaded_and_reported(small_lc, monkeypatch, tmp_path):
    model, data, beta = small_lc

    def failing_em(*args, **kwargs):
        raise FloatingPointError("boom")

    monkeypatch.setattr(engine, "lc_em", failing_em)
    control = Control(model_name="notes", indiv_id="PeID", output_dir=tmp_path)
    result = model.fit(data, beta, control, EstimationSettings(method="em", backend="numpy", hessian="none",
                                                               verbose=False))
    assert "EM failed (boom)" in result.summary_text()
    files = result.save()
    loaded = mtuem.load(files["model"])
    assert loaded.notes == result.notes
    assert loaded.summary_text() == result.summary_text()


def test_a_plain_fit_has_no_notes():
    model = declarative("eq1")
    result = model.fit(subset(80), model.beta(), CONTROL, EstimationSettings(backend="numpy", verbose=False))
    assert result.notes == [] and result.starts is None


# ----------------------------------------------------------------------------------------------
# Multistart
# ----------------------------------------------------------------------------------------------
def test_multistart_same_seed_same_result():
    model = declarative("eq1")
    data = subset(120)
    a = model.fit(data, model.beta(), CONTROL, multistart(seed=7))
    b = model.fit(data, model.beta(), CONTROL, multistart(seed=7))
    assert a.estimates.equals(b.estimates)
    assert a.beta == b.beta
    assert [s.beta.values() for s in a.starts] == [s.beta.values() for s in b.starts]
    c = model.fit(data, model.beta(), CONTROL, multistart(seed=8))
    assert [s.beta.values() for s in a.starts] != [s.beta.values() for s in c.starts]


def _scripted(real, outcomes):
    """``engine.estimate`` whose k-th call has its maximum, code and Hessian eigenvalue replaced."""
    calls = iter(outcomes)

    def estimate(*args, **kwargs):
        est = real(*args, **kwargs)
        outcome = next(calls)
        if outcome is None:
            raise FloatingPointError("scripted failure")
        maximum, code, eig = outcome
        return dataclasses.replace(est, maximum=maximum, code=code, successfulEstimation=code == 0, eigValue=eig)

    return estimate


def test_multistart_keeps_best_non_saddle(monkeypatch):
    model = declarative("eq1")
    outcomes = [(-100.0, 0, -1.0), (-50.0, 0, 1.0), (-10.0, 1, -1.0), None, (-80.0, 0, -2.0)]
    monkeypatch.setattr(engine, "estimate", _scripted(engine.estimate, outcomes))
    result = model.fit(subset(120), model.beta(), CONTROL, multistart(n=5, **NARROW))
    assert result.loglik == -80.0
    assert len(result.starts) == 5 and result.starts[3] is None
    assert result.beta == result.starts[4].beta
    assert result.notes == ["Best of 5 multistart candidates: candidate 5"]
    assert result.starts[4].starts is None


def test_multistart_takes_a_saddle_point_when_non_saddle_is_off(monkeypatch):
    model = declarative("eq1")
    outcomes = [(-100.0, 0, -1.0), (-50.0, 0, 1.0)]
    monkeypatch.setattr(engine, "estimate", _scripted(engine.estimate, outcomes))
    result = model.fit(subset(120), model.beta(), CONTROL, multistart(n=2, non_saddle=False, **NARROW))
    assert result.loglik == -50.0


def test_multistart_all_fail_raises_with_attempts(monkeypatch):
    model = declarative("eq1")

    def failing(*args, **kwargs):
        raise FloatingPointError("scripted failure")

    monkeypatch.setattr(engine, "estimate", failing)
    with pytest.raises(RuntimeError, match="no multistart candidate") as info:
        model.fit(subset(120), model.beta(), CONTROL, multistart(n=3, **NARROW))
    assert info.value.attempts == [None, None, None]


def test_multistart_with_non_saddle_needs_a_hessian():
    model = declarative("eq1")
    settings = EstimationSettings(backend="numpy", hessian="none", verbose=False, multistart=MultiStart(n=2))
    with pytest.raises(ValueError, match="non_saddle"):
        model.fit(subset(50), model.beta(), CONTROL, settings)


def test_multistart_prints_one_line_per_candidate_when_verbose(capsys):
    model = declarative("eq1")
    settings = EstimationSettings(backend="numpy", multistart=MultiStart(n=2, seed=1))
    model.fit(subset(120), model.beta(), CONTROL, settings)
    out = capsys.readouterr().out
    assert "2 candidates available" in out
    assert "Candidate 1:" in out and "Candidate 2:" in out
    assert "Initial LL" not in out  # the estimation of each candidate is quiet


def test_multistart_is_quiet_when_not_verbose(capsys):
    model = declarative("eq1")
    model.fit(subset(120), model.beta(), CONTROL, multistart(n=2))
    assert capsys.readouterr().out == ""


# Candidate draws ---------------------------------------------------------------------------------
def test_draws_only_perturb_free_parameters():
    model = LatentClassMTUEM(
        work_times=["Tw"], free_times=["Tf1", "Tf2"], Tc="Tc", Ec="EcI", w="w", covariance="cholesky",
        n_classes=2, participation=[[True, False], [True, True]], allocation=["female"],
    )
    beta = model.beta()
    draws = _candidates(model, beta, MultiStart(n=20, seed=1, normalization="theta_phi"))
    for name in beta.fixed_names:
        assert (draws[name] == beta[name].value).all(), name
    assert "th2_1" in beta.locked and (draws["th2_1"] == 0).all()
    for name in ("Phi_1", "th1_2", "chk11_1"):
        assert draws[name].nunique() > 10


def test_draws_follow_lower_and_upper_and_reject_unknown_or_fixed_names():
    model = declarative("eq1")
    draws = _candidates(model, model.beta(), MultiStart(n=50, seed=1, lower={"Phi": 2.0}, upper={"Phi": 2.5}))
    assert draws["Phi"].between(2.0, 2.5).all()
    assert draws["thw"].between(-1.0, 1.0).all()
    with pytest.raises(ValueError, match="unknown parameters: gamma"):
        _candidates(model, model.beta(), MultiStart(n=2, lower={"gamma": 0.0}))
    with pytest.raises(ValueError, match="fixed parameters: Theta"):
        _candidates(model, model.beta(), MultiStart(n=2, upper={"Theta": 2.0}))
    with pytest.raises(ValueError, match="Phi"):
        _candidates(model, model.beta(), MultiStart(n=2, lower={"Phi": 3.0}, upper={"Phi": 2.0}))


def test_draws_stay_inside_beta_bounds():
    model = declarative("eq1")
    beta = model.beta().lower("Phi", 0.5).upper("thw", 0.2)
    draws = _candidates(model, beta, MultiStart(n=200, seed=1))
    assert (draws["Phi"] >= 0.5).all() and (draws["thw"] <= 0.2).all()


def test_latent_class_normalization_caps_each_class_separately():
    model = LatentClassMTUEM(
        work_times=["Tw"], free_times=["Tf1", "Tf2"], free_goods=["Ef1"], Tc="Tc", Ec="EcI", w="w",
        covariance="cholesky", n_classes=2, allocation=["female"],
    )
    beta = model.beta().free(["Theta_1", "Theta_2"])
    lower = {f"{p}_{s}": 0.1 for p in ("Theta", "th1", "th2") for s in (1, 2)}
    upper = {f"{p}_{s}": 1.5 for p in ("Theta", "th1", "th2") for s in (1, 2)}
    draws = _candidates(model, beta, MultiStart(n=100, seed=2, lower=lower, upper=upper, normalization="theta_phi"))
    for s in (1, 2):
        assert (draws[f"th1_{s}"] + draws[f"th2_{s}"] <= draws[f"Theta_{s}"] + 1e-12).all()
        # Off-diagonal Cholesky terms start at zero.
        assert (draws[f"chk21_{s}"] == 0).all()
        assert (draws[f"chk11_{s}"] > 0).all()
