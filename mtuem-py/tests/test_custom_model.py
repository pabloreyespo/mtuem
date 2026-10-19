"""mtuem.CustomModel: hand-written likelihoods through the same fit, Result, save and load."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import build, database, declarative, latent_class, renamed
from mtuem import CustomModel
from mtuem.beta import Beta
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import MTUEM

REF = json.loads((Path(__file__).parent / "reference" / "lc_cholesky_ref.json").read_text("utf-8"))
NUMPY = EstimationSettings(backend="numpy", verbose=False)


def eq3_model():
    beta, fixed, probabilities, inputs = build("eq3_cholesky")
    model = CustomModel(probabilities, columns=["Tw", "Tf1", "Ef1", "Tc", "EcI", "w"])
    return model, Beta.from_dict(beta, fixed=fixed)


def lc_model():
    beta, fixed, probabilities, inputs = build("lc_cholesky")
    return CustomModel(probabilities, lc_pars=inputs.lc_pars), Beta.from_dict(beta, fixed=fixed)


def test_custom_model_reproduces_declarative_eq3(tmp_path):
    control = Control(model_name="eq3", indiv_id="PeID", output_dir=tmp_path)
    model, beta = eq3_model()
    custom = model.fit(database(), beta, control, NUMPY)
    decl = declarative("eq3_cholesky")
    dbeta = decl.beta().set(**renamed({k: beta[k].value for k in beta.names}))
    ref = decl.fit(database(), dbeta, control, NUMPY)
    assert custom.converged and ref.converged
    np.testing.assert_allclose(custom.loglik, ref.loglik, rtol=1e-8)
    renamed_est = custom.estimates.rename({"theta_1": "th1", "phi_1": "ph1"})
    np.testing.assert_allclose(renamed_est[ref.estimates.index], ref.estimates, rtol=1e-8)


@pytest.mark.slow
def test_custom_latent_class_with_lc_pars():
    model, beta = lc_model()
    beta = beta.set(**REF["estimate"])
    control = Control(model_name="lc", indiv_id="PeID")
    result = model.fit(database(), beta, control, NUMPY)
    assert result.loglik == pytest.approx(REF["maximum"], abs=1e-3)
    assert result.class_shares().sum() == pytest.approx(1)
    assert result.posterior().shape[1] == 3
    decl = latent_class("lc_cholesky")
    ref = decl.fit(database(), decl.beta().set(**REF["estimate"]), control, NUMPY)
    np.testing.assert_allclose(result.loglik, ref.loglik, rtol=1e-6)
    np.testing.assert_allclose(result.estimates[ref.estimates.index], ref.estimates, rtol=1e-6)


def test_custom_backend_mismatch_raises():
    model, beta = eq3_model()
    control = Control(model_name="eq3", indiv_id="PeID")
    with pytest.raises(ValueError, match="backend"):
        model.fit(database(), beta, control, EstimationSettings(backend="jax", verbose=False))
    jax_model = CustomModel(model._probabilities, backend="jax")
    with pytest.raises(ValueError, match="backend"):
        jax_model.fit(database(), beta, control, NUMPY)


def test_custom_jax_backend_fits():
    beta, fixed, probabilities, _ = build("eq3_cholesky", autodiff=True)
    model = CustomModel(probabilities, backend="jax")
    control = Control(model_name="eq3", indiv_id="PeID")
    jax = model.fit(database(), Beta.from_dict(beta, fixed=fixed), control, EstimationSettings(verbose=False))
    numpy = eq3_model()[0].fit(database(), Beta.from_dict(beta, fixed=fixed), control, NUMPY)
    np.testing.assert_allclose(jax.loglik, numpy.loglik, rtol=1e-6)


def test_custom_em_needs_lc_pars():
    model, beta = eq3_model()
    control = Control(model_name="eq3", indiv_id="PeID")
    with pytest.raises(ValueError, match="latent class"):
        model.fit(database(), beta, control, EstimationSettings(backend="numpy", method="em", verbose=False))


def test_custom_validates_data_columns_and_indiv_id():
    model, beta = eq3_model()
    control = Control(model_name="eq3", indiv_id="PeID")
    with pytest.raises(ValueError, match="lack the columns: Tw"):
        model.fit(database().drop(columns="Tw"), beta, control, NUMPY)
    with pytest.raises(ValueError, match="lack the columns: PeID"):
        model.fit(database().drop(columns="PeID"), beta, control, NUMPY)
    bad = database()
    bad.loc[0, "Tf1"] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        model.fit(bad, beta, control, NUMPY)


def test_custom_without_columns_validates_only_the_id():
    model = CustomModel(eq3_model()[0]._probabilities)
    model._validate_data(database(), Control(indiv_id="PeID"))
    model._validate_data(pd.DataFrame({"a": [1.0]}), Control())


def test_custom_multistart_uses_normalization_groups():
    model, beta = eq3_model()
    model = CustomModel(
        model._probabilities, columns=model.columns(),
        normalization_groups={"cholesky": ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]},
    )
    assert model.normalization_groups()["cholesky"][0] == "chk11"
    ms = MultiStart(n=2, seed=1, normalization="theta_phi", non_saddle=False)
    settings = EstimationSettings(backend="numpy", verbose=False, multistart=ms, hessian="none")
    result = model.fit(database(), beta, Control(model_name="eq3", indiv_id="PeID"), settings)
    assert result.starts is not None and len(result.starts) == 2


def test_result_keeps_a_custom_model():
    model, beta = eq3_model()
    result = model.fit(database(), beta, Control(model_name="eq3", indiv_id="PeID"), NUMPY)
    assert isinstance(result.model, CustomModel)
    assert "CustomModel" in repr(result)
    assert result.predict().shape[0] == result.n_obs


def test_custom_values_of_time_need_explicit_columns():
    model, beta = eq3_model()
    result = model.fit(database(), beta, Control(model_name="eq3", indiv_id="PeID"), NUMPY)
    with pytest.raises(ValueError, match="Tc, Ec and w"):
        result.values_of_time()
    vot = result.values_of_time(Tc="Tc", Ec="EcI", w="w")
    assert list(vot["summary"]["quantity"]) == ["VoL", "VTAW", "wage"]


def test_custom_save_and_load_with_and_without_the_model(tmp_path):
    model, beta = eq3_model()
    control = Control(model_name="eq3", indiv_id="PeID", output_dir=tmp_path)
    result = model.fit(database(), beta, control, NUMPY)
    files = result.save()
    saved = json.loads(files["model"].read_text("utf-8"))
    assert saved["model_class"] == "mtuem.models:CustomModel" and saved["model"] is None
    assert {"model", "beta", "estimates"} <= set(saved)
    for kind, path in files.items():
        assert path.exists(), kind

    bare = mtuem.load(files["model"])
    pd.testing.assert_series_equal(bare.estimates, result.estimates)
    assert bare.summary_text() == result.summary_text()
    with pytest.raises(ValueError, match="model="):
        bare.predict()
    with pytest.raises(ValueError, match="model="):
        mtuem.load(files["model"], data=database())
    bare.save(tmp_path / "again")
    again = json.loads((tmp_path / "again" / "eq3_model.json").read_text("utf-8"))
    assert again["model_class"] == "mtuem.models:CustomModel" and again["model"] is None

    full = mtuem.load(files["model"], data=database(), model=model)
    assert full.predict().shape[0] == result.n_obs
    assert full.loglik == pytest.approx(result.loglik)
    assert isinstance(full.model, CustomModel)


def test_load_rejects_a_model_that_does_not_match(tmp_path):
    model, beta = eq3_model()
    result = model.fit(database(), beta, Control(model_name="eq3", output_dir=tmp_path), NUMPY)
    path = result.save()["model"]
    with pytest.raises(TypeError, match="model"):
        mtuem.load(path, model=object())
    with pytest.raises(ValueError, match="rows"):
        mtuem.load(path, data=database().head(5), model=model)


def test_load_accepts_a_declarative_model_instance(tmp_path):
    decl = declarative("eq1")
    control = Control(model_name="eq1", indiv_id="PeID", output_dir=tmp_path)
    result = decl.fit(database(), decl.beta(), control, NUMPY)
    path = result.save()["model"]
    plain = mtuem.load(path, data=database())
    attached = mtuem.load(path, data=database(), model=decl)
    assert isinstance(plain.model, MTUEM) and isinstance(attached.model, MTUEM)
    assert attached.values_of_time()["summary"].equals(plain.values_of_time()["summary"])
    other = declarative("eq3_cholesky")
    with pytest.raises(ValueError, match="parameters"):
        mtuem.load(path, model=other)


def test_declarative_model_json_unchanged(tmp_path):
    result = declarative("eq1").fit(
        database(), declarative("eq1").beta(), Control(model_name="eq1", output_dir=tmp_path), NUMPY
    )
    saved = json.loads(result.save()["model"].read_text("utf-8"))
    assert saved["model_class"] == "mtuem.models:MTUEM" and saved["model"]["Tc"] == "Tc"


def test_custom_beta_checks():
    model, beta = eq3_model()
    control = Control(model_name="eq3", indiv_id="PeID")
    with pytest.raises(ValueError, match="no parameters"):
        model.fit(database(), Beta(params={}), control, NUMPY)
    bad = CustomModel(model._probabilities, normalization_groups={"cholesky": ["chk99"]})
    with pytest.raises(ValueError, match="chk99"):
        bad.fit(database(), beta, control, NUMPY)


def test_constructor_checks_callables():
    with pytest.raises(TypeError, match="probabilities"):
        CustomModel("not a function")
    with pytest.raises(TypeError, match="lc_pars"):
        CustomModel(lambda p, i, f: {}, lc_pars=3)


@pytest.mark.slow
def test_custom_values_of_time_checks_columns_and_classes(tmp_path):
    model, beta = eq3_model()
    result = model.fit(database(), beta, Control(model_name="eq3", indiv_id="PeID"), NUMPY)
    with pytest.raises(ValueError, match="nope"):
        result.values_of_time(Tc="Tc", Ec="nope", w="w")
    lc, lbeta = lc_model()
    lresult = lc.fit(database(), lbeta.set(**REF["estimate"]), Control(model_name="lc", indiv_id="PeID"), NUMPY)
    with pytest.raises(ValueError, match="latent class"):
        lresult.values_of_time(Tc="Tc", Ec="EcI", w="w")


def test_load_model_mismatches(tmp_path):
    model, beta = eq3_model()
    control = Control(model_name="eq3", output_dir=tmp_path / "c")
    custom_path = model.fit(database(), beta, control, NUMPY).save()["model"]
    with pytest.raises(ValueError, match="model_class"):
        mtuem.load(custom_path, model_class=MTUEM)
    decl = declarative("eq1")
    dcontrol = Control(model_name="eq1", output_dir=tmp_path / "d")
    decl_path = decl.fit(database(), decl.beta(), dcontrol, NUMPY).save()["model"]
    with pytest.raises(ValueError, match="does not match"):
        mtuem.load(decl_path, model=model)
    with pytest.raises(ValueError, match="does not match"):
        mtuem.load(custom_path, model=decl)
