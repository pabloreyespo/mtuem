"""MTUEM.fit, the Result object, saving and loading."""

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import database, declarative
from mtuem.beta import Beta
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import MTUEM
from mtuem.results import Result

NUMPY = EstimationSettings(backend="numpy", verbose=False)
JSON_KEYS = {
    "mtuem_version", "model_class", "model", "beta", "control", "settings", "estimates", "se", "robse",
    "varcov", "robvarcov", "loglik", "loglik_start", "converged", "n_obs", "n_indivs",
}


@pytest.fixture(scope="module")
def fitted(tmp_path_factory):
    folder = tmp_path_factory.mktemp("eq1")
    model = declarative("eq1")
    control = Control(model_name="eq1", description="one equation", indiv_id="PeID", output_dir=folder)
    return model.fit(database(), model.beta(), control, NUMPY)


@pytest.fixture(scope="module")
def saved(fitted):
    return fitted.save()


# ----------------------------------------------------------------------------------------------
# fit
# ----------------------------------------------------------------------------------------------
def test_fit_returns_a_result(fitted):
    assert isinstance(fitted, Result)
    assert fitted.converged
    assert fitted.starts is None
    assert list(fitted.estimates.index) == ["Theta", "Phi", "thw"]
    assert np.isnan(fitted.se["Theta"]) and fitted.se["Phi"] > 0
    assert list(fitted.varcov.index) == ["Phi", "thw"] and list(fitted.robvarcov.columns) == ["Phi", "thw"]
    assert fitted.n_free == 2
    assert fitted.n_obs == fitted.n_indivs == len(database())
    assert fitted.loglik > fitted.loglik_start
    assert fitted.aic == pytest.approx(-2 * fitted.loglik + 4)
    assert fitted.bic == pytest.approx(-2 * fitted.loglik + 2 * np.log(fitted.n_obs))


def test_fit_keeps_the_start_beta_and_configuration(fitted):
    assert isinstance(fitted.beta, Beta)
    assert fitted.beta.values() == {"Theta": 1, "Phi": 1, "thw": 0}
    assert fitted.control.model_name == "eq1"
    assert fitted.settings == NUMPY
    assert fitted.model == declarative("eq1")


def test_fit_does_not_change_the_start_beta():
    model = declarative("eq1")
    beta = model.beta()
    before = beta.model_dump()
    result = model.fit(database(), beta, Control(indiv_id="PeID"), NUMPY)
    assert beta.model_dump() == before
    result.beta.set(Phi=3)
    assert beta["Phi"].value == 1


def test_fit_rejects_beta_with_other_names():
    model = declarative("eq1")
    beta = Beta.from_dict({"Theta": 1, "Phi": 1, "gamma": 0}, fixed=["Theta"])
    with pytest.raises(ValueError, match=r"missing: thw.*extra: gamma"):
        model.fit(database(), beta, Control(indiv_id="PeID"), NUMPY)


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        (EstimationSettings(method="em", verbose=False), "method='em' is not available"),
        (EstimationSettings(multistart=MultiStart(n=5), verbose=False), "multistart is not available"),
    ],
)
def test_fit_rejects_settings_not_available_yet(settings, message):
    model = declarative("eq1")
    with pytest.raises(NotImplementedError, match=message):
        model.fit(database(), model.beta(), Control(indiv_id="PeID"), settings)


def test_fit_applies_beta_bounds():
    model = declarative("eq1")
    beta = model.beta().set(thw=0.5).lower("thw", 0.4)
    result = model.fit(database(), beta, Control(indiv_id="PeID"), NUMPY)
    assert result.estimates["thw"] > 0.4
    assert result.estimates["thw"] == pytest.approx(0.4, abs=1e-3)


def test_fit_default_settings_use_jax_and_print(capsys):
    model = declarative("eq1")
    result = model.fit(database(), model.beta(), Control(indiv_id="PeID"))
    assert result.settings == EstimationSettings()
    assert result.converged
    assert "Initial LL" in capsys.readouterr().out


def test_quiet_fit_prints_nothing(capsys):
    model = declarative("eq1")
    model.fit(database(), model.beta(), Control(indiv_id="PeID"), NUMPY)
    assert capsys.readouterr().out == ""


# ----------------------------------------------------------------------------------------------
# Report and tables
# ----------------------------------------------------------------------------------------------
def test_summary_header_has_no_apollo(fitted, capsys):
    fitted.summary()
    out = capsys.readouterr().out
    assert out.startswith("mtuem ")
    assert "apollo" not in out.lower()
    assert "Model name                                  : eq1" in out
    assert out == fitted.summary_text()


def test_summary_text_prints_covariance_and_correlation_on_request(fitted):
    assert "Classical covariance matrix" not in fitted.summary_text()
    text = fitted.summary_text(print_covar=True, print_corr=True)
    assert "Robust covariance matrix" in text and "Robust correlation matrix" in text


def test_to_frame_columns(fitted):
    frame = fitted.to_frame()
    assert list(frame.columns) == ["Estimate", "Std.err.", "t-ratio(0)", "Rob.std.err.", "Rob.t-ratio(0)"]
    assert list(frame.index) == ["Theta", "Phi", "thw"]
    assert frame.loc["Phi", "Estimate"] == fitted.estimates["Phi"]
    assert frame.loc["Phi", "Rob.t-ratio(0)"] == pytest.approx(fitted.estimates["Phi"] / fitted.robse["Phi"])


def test_delta_method(fitted):
    dm = fitted.delta_method({"ratio": lambda b: b["thw"] / b["Phi"]})
    assert list(dm["Expression"]) == ["ratio"]
    assert dm["Value"].iloc[0] == pytest.approx(round(fitted.estimates["thw"] / fitted.estimates["Phi"], 4))
    classical = fitted.delta_method({"ratio": lambda b: b["thw"] / b["Phi"]}, varcov="classical")
    assert classical["s.e."].iloc[0] != dm["s.e."].iloc[0]


def test_delta_method_rejects_an_unknown_covariance(fitted):
    with pytest.raises(ValueError, match="robust"):
        fitted.delta_method({"Phi": lambda b: b["Phi"]}, varcov="sandwich")


def test_repr_names_the_model(fitted):
    assert repr(fitted).startswith("<Result MTUEM 'eq1': loglik=")


def test_values_of_time_and_predict(fitted):
    vot = fitted.values_of_time()
    assert set(vot) == {"summary", "individual", "delta"}
    assert list(vot["summary"]["quantity"]) == ["VoL", "VTAW", "wage"]
    pred = fitted.predict()
    assert len(pred) == fitted.n_obs and "Tw" in pred.columns


def test_values_of_time_needs_the_theta_phi_form():
    model = MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", utility="alpha_beta")
    result = model.fit(database(), model.beta(), Control(indiv_id="PeID"), NUMPY)
    with pytest.raises(ValueError, match="theta_phi"):
        result.values_of_time()


def test_values_of_time_reject_elasticity_covariates():
    model = MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", elasticity_covariates={"Phi": ["female"]})
    result = model.fit(database(), model.beta(), Control(indiv_id="PeID"), NUMPY)
    with pytest.raises(ValueError, match="elasticity_covariates"):
        result.values_of_time()


def test_correlations_need_a_cholesky_model(fitted):
    with pytest.raises(ValueError, match="cholesky"):
        fitted.correlations()


# ----------------------------------------------------------------------------------------------
# Saving and loading
# ----------------------------------------------------------------------------------------------
def test_save_writes_engine_files_and_json(fitted, saved):
    folder = fitted.control.output_dir
    assert set(saved) == {"output", "estimates", "covar", "robcovar", "corr", "robcorr", "model"}
    assert saved["model"] == folder / "eq1_model.json"
    assert saved["output"] == folder / "eq1_output.txt"
    for path in saved.values():
        assert path.exists()
    data = json.loads(saved["model"].read_text(encoding="utf-8"))
    assert set(data) >= JSON_KEYS
    assert data["model_class"] == "mtuem.models:MTUEM"
    assert data["model"] == declarative("eq1").model_dump(mode="json")
    assert data["beta"] == fitted.beta.model_dump(mode="json")
    assert data["loglik"] == fitted.loglik
    assert data["se"]["Theta"] is None


def test_save_to_another_folder(fitted, tmp_path):
    files = fitted.save(tmp_path / "elsewhere")
    assert files["model"] == tmp_path / "elsewhere" / "eq1_model.json"
    assert files["estimates"].exists()


def test_save_keeps_beta_locks_and_bounds(tmp_path):
    model = declarative("eq1")
    beta = model.beta().set(thw=0.5).lower("thw", 0.1)
    beta.locked["Theta"] = "normalisation"
    result = model.fit(database(), beta, Control(indiv_id="PeID", output_dir=tmp_path), NUMPY)
    loaded = mtuem.load(result.save()["model"])
    assert loaded.beta == beta
    assert loaded.beta.locked == {"Theta": "normalisation"}
    assert loaded.beta["thw"].lower == 0.1


def test_load_round_trip(fitted, saved):
    loaded = mtuem.load(saved["model"], data=database())
    assert loaded.model == fitted.model
    assert loaded.beta == fitted.beta
    # The output folder is where the json was read from, not a path of the saving machine.
    assert loaded.control == fitted.control.model_copy(update={"output_dir": saved["model"].parent})
    assert loaded.settings == fitted.settings
    pd.testing.assert_series_equal(loaded.estimates, fitted.estimates, check_exact=True)
    pd.testing.assert_series_equal(loaded.se, fitted.se, check_exact=True)
    pd.testing.assert_series_equal(loaded.robse, fitted.robse, check_exact=True)
    pd.testing.assert_frame_equal(loaded.varcov, fitted.varcov, check_exact=True)
    pd.testing.assert_frame_equal(loaded.robvarcov, fitted.robvarcov, check_exact=True)
    assert loaded.loglik == fitted.loglik and loaded.loglik_start == fitted.loglik_start
    assert loaded.converged == fitted.converged
    assert (loaded.n_obs, loaded.n_indivs, loaded.n_free) == (fitted.n_obs, fitted.n_indivs, fitted.n_free)
    assert loaded.aic == pytest.approx(fitted.aic) and loaded.bic == pytest.approx(fitted.bic)
    pd.testing.assert_frame_equal(loaded.to_frame(), fitted.to_frame(), check_exact=True)
    vot, again = fitted.values_of_time(), loaded.values_of_time()
    pd.testing.assert_frame_equal(again["summary"], vot["summary"], check_exact=True)
    pd.testing.assert_frame_equal(loaded.predict(), fitted.predict(), check_exact=True)


def _numbers_and_lists(obj, numbers, lists):
    if isinstance(obj, dict):
        for v in obj.values():
            _numbers_and_lists(v, numbers, lists)
    elif isinstance(obj, list):
        lists.append(len(obj))
        for v in obj:
            _numbers_and_lists(v, numbers, lists)
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        numbers.add(float(obj))


def test_saved_json_holds_no_machine_path_and_no_data(fitted, saved):
    text = saved["model"].read_text(encoding="utf-8")
    folder = fitted.control.output_dir
    assert str(folder) not in text and json.dumps(str(folder))[1:-1] not in text
    assert folder.name not in text
    saved_json = json.loads(text)
    assert saved_json["control"]["output_dir"] == "."
    numbers, lists = set(), []
    _numbers_and_lists(saved_json, numbers, lists)
    assert max(lists) < fitted.n_obs
    data = database()[declarative("eq1").columns()].to_numpy(dtype=float).ravel()
    data = {float(v) for v in data if v != round(v)}
    assert not numbers & data


def test_loaded_result_saves_next_to_its_json(saved, tmp_path):
    folder = tmp_path / "moved"
    folder.mkdir()
    path = folder / "eq1_model.json"
    path.write_text(saved["model"].read_text(encoding="utf-8"), encoding="utf-8")
    loaded = mtuem.load(path)
    assert loaded.control.output_dir == folder
    assert loaded.save()["output"] == folder / "eq1_output.txt"


def test_result_keeps_its_own_copy_of_the_model():
    model = declarative("eq1")
    result = model.fit(database(), model.beta(), Control(indiv_id="PeID"), NUMPY)
    model.tau = 100
    assert result.model.tau == 168


def test_save_writes_infinite_bounds_as_no_bound(tmp_path):
    model = declarative("eq1")
    beta = model.beta().lower("Phi", -np.inf).upper("Phi", np.inf)
    result = model.fit(database(), beta, Control(indiv_id="PeID", output_dir=tmp_path), NUMPY)
    loaded = mtuem.load(result.save()["model"])
    assert loaded.beta["Phi"].lower is None and loaded.beta["Phi"].upper is None


@pytest.mark.parametrize("kwargs", [{"tau": np.inf}, {"extra": {"scale": np.inf}}, {"extra": {"s": [1.0, np.nan]}}])
def test_save_rejects_other_infinite_values(tmp_path, kwargs):
    model = MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", **kwargs)
    result = Result(model, model.beta(), Control(output_dir=tmp_path), NUMPY, _fake_estimate(model))
    with pytest.raises(ValueError, match="model.*not finite"):
        result.save()


def _fake_estimate(model):
    from mtuem import _engine as engine
    est = pd.Series(model.beta().values(), dtype=float)
    free = model.beta().free_names
    vc = pd.DataFrame(np.eye(len(free)), index=free, columns=free)
    return engine.Estimate(estimate=est, beta_start=est, fixed=model.beta().fixed_names, control=Control(),
                           maximum=-1.0, LLStart=-2.0, se=est * 0, robse=est * 0, varcov=vc, robvarcov=vc,
                           corrmat=vc, robcorrmat=vc, nObs=1, nIndivs=1, nFreeParams=len(free))


def test_load_rejects_a_class_that_is_not_a_model(saved):
    with pytest.raises(TypeError, match="MTUEM"):
        mtuem.load(saved["model"], model_class=dict)


def test_loaded_summary_matches_the_fitted_one(fitted, saved):
    loaded = mtuem.load(saved["model"])
    assert loaded.summary_text(True, True) == fitted.summary_text(True, True)


def test_load_accepts_the_folder(saved):
    loaded = mtuem.load(saved["model"].parent)
    assert loaded.control.model_name == "eq1"


def test_load_rejects_a_folder_without_one_model_file(tmp_path):
    with pytest.raises(ValueError, match="_model.json"):
        mtuem.load(tmp_path)


def test_load_without_data_methods_raise(saved):
    loaded = mtuem.load(saved["model"])
    for method in (loaded.values_of_time, loaded.predict):
        with pytest.raises(ValueError, match="pass data= to mtuem.load"):
            method()
    # Methods that only need the estimates work without data.
    assert len(loaded.delta_method({"Phi": lambda b: b["Phi"]})) == 1


def test_load_rejects_data_of_another_size(saved):
    with pytest.raises(ValueError, match="rows"):
        mtuem.load(saved["model"], data=database().head(10))


def test_load_unimportable_class_raises_with_hint(saved, tmp_path):
    data = json.loads(saved["model"].read_text(encoding="utf-8"))
    data["model_class"] = "my_models.missing:CustomMTUEM"
    path = tmp_path / "eq1_model.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ImportError, match=r"my_models\.missing:CustomMTUEM.*model_class="):
        mtuem.load(path)
    data["model_class"] = "mtuem.models:NoSuchClass"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ImportError, match="model_class="):
        mtuem.load(path)
    loaded = mtuem.load(path, model_class=MTUEM)
    assert loaded.model == declarative("eq1")


def test_load_rejects_a_model_that_does_not_match_the_estimates(saved, tmp_path):
    data = json.loads(saved["model"].read_text(encoding="utf-8"))
    data["model"]["free_times"] = ["Tf1"]
    path = tmp_path / "eq1_model.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match=r"th1.*estimates"):
        mtuem.load(path)
    data = json.loads(saved["model"].read_text(encoding="utf-8"))
    data["beta"]["params"]["gamma"] = {"value": 0.0}
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match=r"saved beta.*gamma"):
        mtuem.load(path)


def test_values_of_time_with_constant_budgets(fitted, saved, tmp_path):
    # A model whose Tc is a constant: values of time use that constant for every row.
    data = json.loads(saved["model"].read_text(encoding="utf-8"))
    data["model"]["Tc"] = 60.0
    path = tmp_path / "eq1_model.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    loaded = mtuem.load(path, data=database())
    vot = loaded.values_of_time()
    expected = mtuem.values_of_time(
        fitted._estimate, SimpleNamespace(database=database().assign(Tc=60.0)), Ec="EcI", silent=True
    )
    pd.testing.assert_frame_equal(vot["summary"], expected["summary"])


def test_save_rejects_an_infinite_beta_value(tmp_path):
    model = declarative("eq1")
    beta = model.beta().set(Phi=np.inf)
    result = Result(model, beta, Control(output_dir=tmp_path), NUMPY, _fake_estimate(model))
    with pytest.raises(ValueError, match="beta.*not finite"):
        result.save()
