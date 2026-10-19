"""The declarative single-class model mtuem.MTUEM."""

import json
import warnings
from pathlib import Path

import jax
import numpy as np
import pandas as pd
import pytest

from cases import SINGLE, database, declarative, renamed
from mtuem import _engine as engine
from mtuem import likelihood, likelihood_jax
from mtuem._engine import FallbackWarning
from mtuem._engine.latent import _loglik_per_indiv
from mtuem.beta import Beta
from mtuem.config import Control
from mtuem.models import MTUEM

REF = Path(__file__).parent / "reference"
CONTROL = Control(model_name="m", indiv_id="PeID")


def load_ref(case: str) -> dict:
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


def eq3(**kw) -> MTUEM:
    args = dict(work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"], Tc="Tc", Ec="EcI", w="w")
    return MTUEM(**(args | kw))


def compiled(model: MTUEM, backend: str, data: pd.DataFrame | None = None):
    """(prob_fn, inputs, autodiff) of ``model`` on ``data`` (the README database by default)."""
    data = database() if data is None else data
    prob_fn, lc_pars, autodiff = model._compile(backend)
    beta = model.beta()
    inputs = engine.make_inputs(data, beta.values(), beta.fixed_names, CONTROL, lc_pars, autodiff=autodiff)
    return prob_fn, inputs, autodiff


def rows_ll(prob_fn, inputs, values: dict) -> np.ndarray:
    return np.asarray(prob_fn(values, inputs, "estimate")["model"], dtype=float)


def total_ll(prob_fn, inputs, values: dict) -> float:
    return float(np.sum(_loglik_per_indiv(prob_fn(values, inputs, "estimate")["model"], inputs)))


# ----------------------------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------------------------
def test_beta_names_match_r_convention():
    beta = eq3(covariance="cholesky").beta()
    assert beta.names == ["Theta", "Phi", "thw", "th1", "ph1", "chk11", "chk21", "chk22", "chk31", "chk32", "chk33"]


def test_beta_defaults_theta_phi():
    beta = eq3(covariance="cholesky").beta()
    assert isinstance(beta, Beta)
    assert beta.fixed_names == ["Theta"]
    assert beta.values() == {
        "Theta": 1, "Phi": 1, "thw": 0, "th1": 0.3, "ph1": 0.3,
        "chk11": 10, "chk21": 0, "chk22": 10, "chk31": 0, "chk32": 0, "chk33": 10,
    }


def test_beta_is_fresh_on_every_call():
    model = eq3()
    first = model.beta()
    first.set(Phi=5)
    assert model.beta()["Phi"].value == 1


def test_beta_names_diagonal_full_and_profiled():
    two = dict(free_times=["Tf1", "Tf2"], free_goods=["Ef1"])
    assert eq3(**two).beta().names == ["Theta", "Phi", "thw", "th1", "th2", "ph1"]
    assert eq3(covariance="diagonal").beta().names[-3:] == ["sig1", "sig2", "sig3"]
    full = eq3(covariance="full").beta()
    assert full.names[5:] == ["sig1", "sig2", "sig3", "rho12", "rho13", "rho23"]
    assert full.values()["sig2"] == 10 and full.values()["rho13"] == 0


def test_beta_names_with_elasticity_covariates_follow_their_base():
    model = eq3(elasticity_covariates={"Phi": ["female"], "th1": ["female", "older45"]})
    assert model.beta().names == ["Theta", "Phi", "Phi_female", "thw", "th1", "th1_female", "th1_older45", "ph1"]
    assert model.beta().values()["th1_older45"] == 0


def test_beta_alpha_beta_names_and_normalisation():
    beta = eq3(utility="alpha_beta", covariance="cholesky").beta()
    assert beta.names[:4] == ["alpha", "beta", "th1", "ph1"]
    assert beta.fixed_names == []
    assert beta.values()["alpha"] == 0.25 and beta.values()["beta"] == 0.25


def test_beta_additive_quadratic_names_and_normalisation():
    model = MTUEM(
        work_times=["Tw"], free_times=["T1", "Ti"], free_goods=["E1", "Ei"], Tc="Tc", Ec="EcI", w="w",
        utility="additive_quadratic", covariance="cholesky",
        extra={"omitted_times": ["Ti"], "omitted_goods": ["Ei"]},
    )
    beta = model.beta()
    assert beta.names == [
        "theta_w", "beta_w", "th1", "beta1", "th2", "beta2", "ph1", "eta1", "ph2", "eta2",
        "chk11", "chk21", "chk22", "chk31", "chk32", "chk33",
    ]
    assert beta.fixed_names == ["theta_w", "eta1"]
    assert beta.values()["theta_w"] == 0 and beta.values()["eta1"] == 1


def test_columns_lists_every_data_column():
    model = eq3(
        free_goods=["Ef1", "Ef2"], goods_cost={"Ef1": 1, "Ef2": "price2"}, Tc=60.0,
        elasticity_covariates={"Phi": ["female"], "ph2": ["female", "older45"]},
    )
    assert model.columns() == ["Tw", "Tf1", "Ef1", "Ef2", "EcI", "w", "price2", "female", "older45"]


def test_normalization_groups():
    groups = eq3(covariance="cholesky", elasticity_covariates={"Phi": ["female"]}).normalization_groups()
    assert groups == {
        "work_elasticities": ["Theta", "Phi", "thw"],
        "times_elasticities": ["th1"],
        "goods_elasticities": ["ph1"],
        "sig": [],
        "rho": [],
        "cholesky": ["chk11", "chk21", "chk22", "chk31", "chk32", "chk33"],
    }
    full = eq3(covariance="full").normalization_groups()
    assert full["sig"] == ["sig1", "sig2", "sig3"] and full["rho"] == ["rho12", "rho13", "rho23"]


def test_model_round_trips_through_json():
    model = eq3(
        covariance="full", goods_cost={"Ef1": "price"}, tau=150, elasticity_covariates={"Phi": ["female"]},
        extra={"componentName": "x"},
    )
    again = MTUEM.model_validate(json.loads(json.dumps(model.model_dump())))
    assert again == model
    assert again.Tc == "Tc" and again.tau == 150


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"work_times": []}, "work_times"),
        ({"work_times": ["Tw", "Tw2"]}, "work_times"),
        ({"free_times": ["Tw"]}, "more than once"),
        ({"elasticity_covariates": {"Gamma": ["female"]}}, "Gamma"),
        ({"elasticity_covariates": {"th2": ["female"]}}, "th2"),
        ({"goods_cost": {"Ef9": 1}}, "Ef9"),
        ({"utility": "alpha_beta", "optimal_tw": False}, "optimal_tw"),
        ({"utility": "additive_quadratic", "elasticity_covariates": {"th1": ["female"]}}, "additive_quadratic"),
        ({"utility": "additive_quadratic", "goods_cost": {"Ef1": "price"}}, "additive_quadratic"),
        ({"extra": {"sig": [1]}}, "sig"),
        ({"bogus": 1}, "bogus"),
        ({"extra": {"omitted_times": ["Tf1"]}}, "additive_quadratic"),
        ({"elasticity_covariates": {"Phi": ["female", "female"]}}, "more than once"),
    ],
)
def test_invalid_specifications_are_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        eq3(**kwargs)


# ----------------------------------------------------------------------------------------------
# Compiled likelihood
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("backend", ["numpy", "jax"])
@pytest.mark.parametrize("case", SINGLE)
@pytest.mark.parametrize("point", ["start", "test", "estimate"])
def test_compiled_loglik_matches_reference(case, point, backend):
    ref = load_ref(case)
    beta_key = {"start": "beta_start", "test": "beta_test", "estimate": "estimate"}[point]
    ll_key = {"start": "ll_start", "test": "ll_test", "estimate": "ll_at_estimate"}[point]
    model = declarative(case)
    prob_fn, inputs, autodiff = compiled(model, backend)
    values = renamed(ref[beta_key])
    assert sorted(values) == sorted(model.beta().names)
    assert autodiff == (backend == "jax")
    assert total_ll(prob_fn, inputs, values) == pytest.approx(ref[ll_key], rel=1e-10, abs=1e-8)


def test_compile_returns_no_lc_pars():
    _, lc_pars, autodiff = eq3()._compile("numpy")
    assert lc_pars is None and autodiff is False


CHOL = {"chk11": 7, "chk21": -5, "chk22": 5, "chk31": 14, "chk32": -13, "chk33": 37}
COVARIATE_CASES = {
    "theta_phi": ({"Theta": 1, "Phi": 0.6, "thw": 0.02, "th1": 0.5, "ph1": 0.2},
                  {"Phi": 0.1, "thw": -0.03, "th1": 0.07, "ph1": -0.05}),
    "alpha_beta": ({"alpha": 0.3, "beta": 0.2, "th1": 0.4, "ph1": 0.2},
                   {"alpha": 0.02, "beta": -0.03, "th1": 0.05, "ph1": -0.04}),
}


@pytest.mark.parametrize("backend", ["numpy", "jax"])
@pytest.mark.parametrize("utility", ["theta_phi", "alpha_beta"])
def test_elasticity_covariates_enter_linearly(backend, utility):
    # With a Cholesky covariance each row's log-likelihood depends on that row alone, so the
    # model with covariates equals the model without them evaluated on each group with the
    # shifted elasticities.
    data = database()
    base, shift = COVARIATE_CASES[utility]
    base = base | CHOL
    with_cov = eq3(utility=utility, covariance="cholesky", elasticity_covariates={k: ["female"] for k in shift})
    plain = eq3(utility=utility, covariance="cholesky")
    values = base | {f"{k}_female": v for k, v in shift.items()}

    prob_fn, inputs, _ = compiled(with_cov, backend, data)
    ll = rows_ll(prob_fn, inputs, values)

    female = data["female"].to_numpy()
    expected = np.empty(len(data))
    for flag, point in ((False, base), (True, base | {k: base[k] + v for k, v in shift.items()})):
        sub = data[female == flag].reset_index(drop=True)
        pf, sub_inputs, _ = compiled(plain, backend, sub)
        expected[female == flag] = rows_ll(pf, sub_inputs, point)
    np.testing.assert_allclose(ll, expected, rtol=1e-12)

    zero = rows_ll(prob_fn, inputs, base | {f"{k}_female": 0.0 for k in shift})
    pf, all_inputs, _ = compiled(plain, backend, data)
    np.testing.assert_allclose(zero, rows_ll(pf, all_inputs, base), rtol=1e-14)


def test_jax_backend_is_differentiable_with_covariates():
    model = eq3(covariance="cholesky", elasticity_covariates={"Phi": ["female"], "th1": ["older45"]})
    prob_fn, inputs, autodiff = compiled(model, "jax")
    assert autodiff
    values = model.beta().values() | {"Phi_female": 0.1, "th1_older45": 0.05}
    names = list(values)

    def total(x):
        return prob_fn(dict(zip(names, x)), inputs, "estimate")["model"].sum()

    x = np.array([values[n] for n in names], dtype=float)
    grad = jax.jit(jax.grad(total))(x)
    assert np.all(np.isfinite(np.asarray(grad)))
    np.testing.assert_allclose(float(jax.jit(total)(x)), float(total(x)), rtol=1e-12)


@pytest.mark.parametrize("covariance", ["profiled", "diagonal", "full", "cholesky"])
def test_alpha_beta_compiles_to_the_alpha_beta_likelihood(covariance):
    data = database()
    model = eq3(utility="alpha_beta", covariance=covariance)
    values = model.beta().values() | {"alpha": 0.3, "beta": 0.2, "th1": 0.4, "ph1": 0.2}
    settings = {
        "work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"], "goods_cost": {"Ef1": 1},
        "work_elasticities": {"alpha": 0.3, "beta": 0.2},
        "times_elasticities": {"th1": 0.4}, "goods_elasticities": {"ph1": 0.2},
        "Tc": data["Tc"].to_numpy(), "Ec": data["EcI"].to_numpy(), "w": data["w"].to_numpy(), "tau": 168,
    }
    if covariance in ("diagonal", "full"):
        settings["sig"] = [10.0, 10.0, 10.0]
    if covariance == "full":
        settings["rho"] = [0.0, 0.0, 0.0]
    if covariance == "cholesky":
        settings["cholesky"] = [10.0, 0, 10.0, 0, 0, 10.0]
    for backend in ("numpy", "jax"):
        prob_fn, inputs, _ = compiled(model, backend, data)
        direct = likelihood.alpha_beta(settings, "estimate", inputs)
        np.testing.assert_allclose(rows_ll(prob_fn, inputs, values), direct, rtol=1e-10)


@pytest.mark.parametrize("backend", ["numpy", "jax"])
def test_additive_quadratic_compiles_to_the_additive_quadratic_likelihood(backend):
    data = database().assign(T1=lambda d: d["Tf1"], Ti=lambda d: d["Tf2"], E1=lambda d: d["Ef1"],
                             Ei=lambda d: d["Ef2"] + d["Ef3"])
    model = MTUEM(
        work_times=["Tw"], free_times=["T1", "Ti"], free_goods=["E1", "Ei"], Tc="Tc", Ec="EcI", w="w",
        utility="additive_quadratic", covariance="cholesky",
        extra={"omitted_times": ["Ti"], "omitted_goods": ["Ei"]},
    )
    values = model.beta().values() | {"beta_w": 1.1, "th2": 1.0, "beta2": 0.9, "ph1": 1.0, "ph2": 1.2, "eta2": 0.8}
    settings = {
        "work_times": ["Tw"], "free_times": ["T1", "Ti"], "free_goods": ["E1", "Ei"],
        "omitted_times": ["Ti"], "omitted_goods": ["Ei"], "goods_cost": {"E1": 1, "Ei": 1},
        "work_elasticity": {"theta_w": 0}, "work_satiation": {"beta_w": 1.1},
        "times_elasticities": {"th1": 0.3, "th2": 1.0}, "times_satiations": {"beta1": 1, "beta2": 0.9},
        "goods_elasticities": {"ph1": 1.0, "ph2": 1.2}, "goods_satiations": {"eta1": 1, "eta2": 0.8},
        "cholesky": [10.0, 0, 10.0, 0, 0, 10.0],
        "Tc": data["Tc"].to_numpy(), "Ec": data["EcI"].to_numpy(), "w": data["w"].to_numpy(), "tau": 168,
    }
    prob_fn, inputs, _ = compiled(model, backend, data)
    direct = likelihood.additive_quadratic(settings, "estimate", inputs)
    np.testing.assert_allclose(rows_ll(prob_fn, inputs, values), direct, rtol=1e-10)


def test_constant_context_values_are_broadcast():
    data = database()
    model = eq3(covariance="cholesky", tau=150, Tc=50.0)
    prob_fn, inputs, _ = compiled(model, "numpy", data)
    values = model.beta().values()
    settings = {
        "work_times": ["Tw"], "free_times": ["Tf1"], "free_goods": ["Ef1"],
        "work_elasticities": {"Theta": 1, "Phi": 1, "thw": 0},
        "times_elasticities": {"th1": 0.3}, "goods_elasticities": {"ph1": 0.3},
        "cholesky": [10.0, 0, 10.0, 0, 0, 10.0],
        "Tc": np.full(len(data), 50.0), "Ec": data["EcI"].to_numpy(), "w": data["w"].to_numpy(), "tau": 150,
    }
    np.testing.assert_allclose(rows_ll(prob_fn, inputs, values), likelihood.theta_phi(settings, "estimate", inputs),
                               rtol=1e-12)


def test_jax_backend_falls_back_to_numpy_for_forms_without_a_jax_version(monkeypatch):
    monkeypatch.setattr(likelihood_jax, "PORTED", frozenset({"theta_phi"}))
    monkeypatch.setattr(likelihood_jax, "_warned", set())
    with pytest.warns(FallbackWarning, match="alpha_beta"):
        prob_fn, lc_pars, autodiff = eq3(utility="alpha_beta")._compile("jax")
    assert autodiff is False and lc_pars is None
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        eq3(utility="alpha_beta")._compile("jax")  # warns once only
    data = database()
    inputs = engine.make_inputs(data, eq3(utility="alpha_beta").beta().values(), [], CONTROL)
    assert isinstance(prob_fn(eq3(utility="alpha_beta").beta().values(), inputs, "estimate")["model"], np.ndarray)


def test_compile_rejects_an_unknown_backend():
    with pytest.raises(ValueError, match="backend"):
        eq3()._compile("torch")


# ----------------------------------------------------------------------------------------------
# Data validation and fit
# ----------------------------------------------------------------------------------------------
def test_fit_rejects_nan_column():
    data = database()
    data.loc[3, "Tf1"] = np.nan
    data.loc[5, "w"] = np.inf
    with pytest.raises(ValueError, match="Tf1.*w|w.*Tf1"):
        eq3()._validate_data(data)
    with pytest.raises(ValueError, match="Tf1"):
        eq3().fit(data, eq3().beta(), CONTROL)


def test_fit_rejects_missing_column():
    data = database().drop(columns=["Ef1", "EcI"])
    with pytest.raises(ValueError, match="Ef1.*EcI|EcI.*Ef1"):
        eq3()._validate_data(data)
    with pytest.raises(ValueError, match="Ef1"):
        eq3().fit(data, eq3().beta(), CONTROL)


def test_fit_rejects_non_numeric_column():
    data = database().assign(Tf1=lambda d: d["Tf1"].astype(str))
    with pytest.raises(ValueError, match="numeric.*Tf1"):
        eq3()._validate_data(data)


def test_fit_rejects_missing_indiv_id():
    data = database().drop(columns=["PeID"])
    with pytest.raises(ValueError, match="PeID"):
        eq3()._validate_data(data, CONTROL)


def test_validate_accepts_boolean_covariates():
    eq3(elasticity_covariates={"Phi": ["female"]})._validate_data(database(), CONTROL)


def test_fit_is_not_implemented_yet():
    with pytest.raises(NotImplementedError, match="MTUEM.fit is not implemented yet"):
        eq3().fit(database(), eq3().beta(), CONTROL)
