"""The declarative latent class model mtuem.LatentClassMTUEM and its class_loglik hook."""

import json
from pathlib import Path
from types import SimpleNamespace

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import pytest

import mtuem
from cases import CHK, LATENT, build, database, latent_class
from mtuem import _engine as engine
from mtuem._engine.latent import _loglik_per_indiv
from mtuem.config import Control, EstimationSettings
from mtuem.models import Data, LatentClassMTUEM
from mtuem.results import Result

REF = Path(__file__).parent / "reference"
CONTROL = Control(model_name="lc", indiv_id="PeID")
NUMPY = EstimationSettings(backend="numpy", verbose=False)
JAX = EstimationSettings(backend="jax", verbose=False)


def load_ref(case: str) -> dict:
    return json.loads((REF / f"{case}_ref.json").read_text(encoding="utf-8"))


def num(series) -> pd.Series:
    return pd.to_numeric(pd.Series(series), errors="coerce")


def compiled(model, backend: str, data: pd.DataFrame | None = None, control: Control = CONTROL):
    data = database() if data is None else data
    prob_fn, lc_pars, autodiff = model._compile(backend)
    beta = model.beta()
    inputs = engine.make_inputs(data, beta.values(), beta.fixed_names, control, lc_pars, autodiff=autodiff)
    return prob_fn, inputs


def total_ll(prob_fn, inputs, values: dict) -> float:
    return float(np.sum(_loglik_per_indiv(prob_fn(values, inputs, "estimate")["model"], inputs)))


def hobbies(**kw) -> LatentClassMTUEM:
    """Four equations as in the article's hobbies model, class 1 without the second free time."""
    args = dict(
        work_times=["Tw"], free_times=["Tf1", "Tf2", "Td"], Tc="Tc", Ec="EcI", w="w",
        covariance="cholesky", n_classes=2, participation=[[True, False, True], [True, True, True]],
        allocation=["female", "older45"],
    )
    return LatentClassMTUEM(**(args | kw))


def result_at(model, values: dict, data: pd.DataFrame | None = None, backend: str = "numpy",
              control: Control = CONTROL) -> Result:
    """A ``Result`` at fixed parameter values (no estimation), with a made-up covariance."""
    prob_fn, inputs = compiled(model, backend, data, control)
    beta = model.beta().set(**{k: v for k, v in values.items() if k not in model.beta().locked})
    est = pd.Series(values, dtype=float)
    free = [n for n in est.index if n not in beta.fixed_names]
    varcov = pd.DataFrame(np.eye(len(free)) * 0.01, index=free, columns=free)
    se = pd.Series(np.nan, index=est.index)
    se[free] = 0.1
    ll = total_ll(prob_fn, inputs, values)
    estimate = engine.Estimate(
        estimate=est, beta_start=est, fixed=beta.fixed_names, control=control, maximum=ll, LLStart=ll,
        LLout={"model": ll}, se=se, robse=se, varcov=varcov, robvarcov=varcov, corrmat=varcov * 100,
        robcorrmat=varcov * 100, nObs=inputs.nObs, nIndivs=inputs.nIndivs, nFreeParams=len(free),
    )
    return Result(model, beta, control, EstimationSettings(backend=backend, verbose=False), estimate, inputs, prob_fn)


# ----------------------------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("case", LATENT)
def test_beta_names_and_fixed_match_r(case):
    ref = load_ref(case)
    beta = latent_class(case).beta()
    assert sorted(beta.names) == sorted(ref["beta_start"])
    assert sorted(beta.fixed_names) == sorted(ref["fixed"])
    assert beta.locked == {}


def test_beta_order_and_defaults_per_class():
    beta = latent_class("lc_cholesky").beta()
    assert beta.names[:12] == [
        "Theta_1", "Phi_1", "thw_1", "th1_1", "ph1_1", *[f"{c}_1" for c in CHK], "asc_1",
    ]
    assert beta.names[12:14] == ["x_female_1", "x_older45_1"]
    assert beta["chk11_2"].value == 10 and beta["chk21_2"].value == 0 and beta["asc_2"].value == 0
    assert beta["x_older45_2"].value == 0 and not beta["x_older45_2"].fixed


def test_shared_names_have_no_suffix():
    beta = latent_class("lc_shared").beta()
    assert beta.names[-6:] == CHK
    assert "chk11_1" not in beta


def test_reference_class_fixes_its_allocation_without_locking():
    beta = hobbies(reference_class=2).beta()
    assert {"asc_2", "x_female_2", "x_older45_2"} <= set(beta.fixed_names)
    assert "asc_1" in beta.free_names
    beta.free("asc_2")
    assert "asc_2" in beta.free_names


def test_absent_activity_parameters_are_locked():
    model = hobbies()
    beta = model.beta()
    assert model.absent_times(1) == ["Tf2"] and model.absent_times(2) == []
    # Equation numbers are those of the full model (1 work, 2 Tf1, 3 Tf2, 4 Td), as in R.
    locked = ["th2_1", "chk31_1", "chk32_1", "chk33_1", "chk43_1"]
    assert sorted(beta.locked) == sorted(locked)
    assert set(beta.locked.values()) == {"class 1 does not take part in Tf2"}
    for name in locked:
        assert beta[name].fixed and beta[name].value == 0
    assert beta["chk44_1"].value == 10 and beta["chk41_1"].value == 0 and not beta["chk41_1"].fixed
    assert beta["chk33_2"].value == 10 and not beta["th2_2"].fixed
    # The fixed set of R/models/hobbies_cholesky.R (get_apollo_fixed) without its covariates.
    assert sorted(beta.fixed_names) == sorted(
        ["Theta_1", "Theta_2", *locked, "asc_1", "x_female_1", "x_older45_1"]
    )
    with pytest.raises(ValueError, match="class 1 does not take part in Tf2"):
        beta.free("th2_1")
    with pytest.raises(ValueError, match="class 1 does not take part in Tf2"):
        beta.set(th2_1=0.3)
    with pytest.raises(ValueError, match="class 1 does not take part in Tf2"):
        beta.free("th2")


def test_absent_activity_locks_its_covariate_terms():
    beta = hobbies(elasticity_covariates={"th2": ["female"]}).beta()
    assert beta.locked["th2_female_1"] == "class 1 does not take part in Tf2"
    assert "th2_female_2" in beta.free_names


def test_absent_activity_diagonal_and_full_covariances():
    diag = hobbies(covariance="diagonal").beta()
    assert sorted(diag.locked) == ["sig3_1", "th2_1"]
    assert diag["sig3_1"].value == 0 and diag["sig4_1"].value == 10
    full = hobbies(covariance="full").beta()
    assert sorted(full.locked) == ["rho13_1", "rho23_1", "rho34_1", "sig3_1", "th2_1"]


def test_columns_include_allocation_covariates():
    assert hobbies().columns() == ["Tw", "Tf1", "Tf2", "Td", "Tc", "EcI", "w", "female", "older45"]


def test_model_round_trips_through_json():
    model = hobbies(shared=["chk11"], reference_class=2)
    assert LatentClassMTUEM.model_validate_json(model.model_dump_json()) == model


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"n_classes": 1}, "n_classes"),
        ({"participation": [[True, True, True]]}, "one row per class"),
        ({"participation": [[True, True], [True, True]]}, "one entry per free time"),
        ({"reference_class": 3}, "reference_class"),
        ({"allocation": ["female", "female"]}, "more than once"),
        ({"allocation": ["Female", "female"]}, "x_female"),
        ({"shared": ["gamma"]}, "shared names unknown parameters: gamma"),
        ({"shared": ["th2"]}, "class 1 does not take part in Tf2"),
        ({"shared": ["chk33"]}, "class 1 does not take part in Tf2"),
        ({"shared": ["Phi", "Phi"]}, "more than once"),
        ({"covariance": "profiled"}, "profiled"),
        ({"participation": None, "covariance": "profiled", "utility": "alpha_beta"}, "theta_phi"),
        ({"utility": "additive_quadratic"}, "additive_quadratic"),
    ],
)
def test_invalid_specifications_are_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        hobbies(**kwargs)


def test_fit_rejects_a_beta_that_frees_an_absent_parameter():
    model = hobbies()
    values = model.beta().values()
    beta = mtuem.beta.Beta.from_dict(values, fixed=[n for n in model.beta().fixed_names if n != "th2_1"])
    with pytest.raises(ValueError, match="th2_1.*class 1 does not take part in Tf2"):
        model.fit(database(), beta, CONTROL, NUMPY)


# ----------------------------------------------------------------------------------------------
# Likelihood
# ----------------------------------------------------------------------------------------------
@pytest.mark.parametrize("backend", ["numpy", "jax"])
@pytest.mark.parametrize("point", ["beta_start", "beta_test", "estimate"])
@pytest.mark.parametrize("case", LATENT)
def test_compiled_lc_loglik_matches_reference(case, point, backend):
    ref = load_ref(case)
    prob_fn, inputs = compiled(latent_class(case), backend)
    values = ref[point]
    ll = total_ll(prob_fn, inputs, values)
    _, _, pf, hand = build(case, autodiff=backend == "jax")
    assert ll == pytest.approx(total_ll(pf, hand, values), rel=1e-10)
    r = {"beta_start": "ll_start", "beta_test": "ll_test", "estimate": "ll_at_estimate"}[point]
    assert ll == pytest.approx(ref[r], rel=1e-10, abs=1e-8)


def test_compiled_lc_returns_lc_pars_and_class_components():
    prob_fn, inputs = compiled(latent_class("lc_cholesky"), "numpy")
    P = prob_fn(load_ref("lc_cholesky")["estimate"], inputs, "estimate")
    assert list(P) == ["class_1", "class_2", "model"]
    assert len(inputs.attach(load_ref("lc_cholesky")["estimate"]).pi_values) == 2


def test_absent_activity_drops_its_equation():
    # Class 1 of the hobbies model is the three equation model without Tf2.
    model = hobbies()
    beta = model.beta()
    values = beta.values() | {"th1_1": 0.4, "th3_1": 0.2, "chk21_1": 1.0, "chk41_1": 2.0, "chk42_1": 0.5}
    prob_fn, inputs = compiled(model, "numpy")
    class_1 = np.asarray(prob_fn(values, inputs, "estimate")["class_1"])
    single = mtuem.MTUEM(work_times=["Tw"], free_times=["Tf1", "Td"], Tc="Tc", Ec="EcI", w="w",
                         covariance="cholesky")
    sv = {"Theta": 1, "Phi": 1, "thw": 0, "th1": 0.4, "th2": 0.2,
          "chk11": 10, "chk21": 1.0, "chk22": 10, "chk31": 2.0, "chk32": 0.5, "chk33": 10}
    spf, _, _ = single._compile("numpy")
    np.testing.assert_allclose(class_1, np.asarray(spf(sv, inputs, "estimate")["model"]), rtol=1e-12)


@pytest.mark.parametrize("backend", ["numpy", "jax"])
def test_elasticity_covariates_are_per_class(backend):
    model = latent_class("lc_cholesky").model_copy(update={"elasticity_covariates": {"Phi": ["female"]}})
    model = LatentClassMTUEM.model_validate(model.model_dump())
    beta = model.beta()
    assert "Phi_female_1" in beta and "Phi_female_2" in beta
    values = load_ref("lc_cholesky")["estimate"] | {"Phi_female_1": 0.1, "Phi_female_2": -0.05}
    prob_fn, inputs = compiled(model, backend)
    ll = total_ll(prob_fn, inputs, values)
    # Same as a hand-built per-row Phi.
    data = database()
    _, _, pf, hand = build("lc_cholesky", autodiff=backend == "jax")
    assert ll != pytest.approx(total_ll(pf, hand, values), rel=1e-6)
    seen = {}

    class Spy(LatentClassMTUEM):
        def class_loglik(self, s, params, data):
            seen[s] = (np.asarray(params.Phi), params.Phi_female)
            return super().class_loglik(s, params, data)

    spy = Spy.model_validate(model.model_dump())
    spf, sinputs = compiled(spy, "numpy")
    assert total_ll(spf, sinputs, values) == pytest.approx(ll, rel=1e-10)
    phi = values["Phi_1"] + 0.1 * data["female"].to_numpy(float)
    np.testing.assert_allclose(seen[1][0], phi)
    assert seen[1][1] == 0.1


@pytest.mark.parametrize("backend", ["numpy", "jax"])
def test_class_loglik_override_is_used(backend):
    ref = load_ref("lc_cholesky")

    class Shifted(LatentClassMTUEM):
        shift: float = 0.5

        def class_loglik(self, s, params, data):
            return super().class_loglik(s, params, data) - self.shift

    base = latent_class("lc_cholesky")
    shifted = Shifted.model_validate(base.model_dump())
    prob_fn, inputs = compiled(base, backend)
    sprob_fn, sinputs = compiled(shifted, backend)
    values = ref["estimate"]
    # Every class loses the same constant per row, so the mixture loses it too.
    n = len(database())
    assert total_ll(sprob_fn, sinputs, values) == pytest.approx(total_ll(prob_fn, inputs, values) - 0.5 * n, rel=1e-10)


def test_class_loglik_override_sees_params_and_data_and_runs_under_jit():
    calls = []

    class Penalized(LatentClassMTUEM):
        def class_loglik(self, s, params, data):
            calls.append((s, type(data.xp)))
            assert isinstance(params, SimpleNamespace) and isinstance(data, Data)
            ll = super().class_loglik(s, params, data)
            return ll - 0.01 * s * data["Tf1"] - 0.0 * data[["Tf1", "Ef1"]].sum(axis=1)

    model = Penalized.model_validate(latent_class("lc_cholesky").model_dump())
    prob_fn, inputs = compiled(model, "jax")
    values = load_ref("lc_cholesky")["estimate"]
    free = [n for n in values if n not in model.beta().fixed_names]

    def total(x):
        b = dict(values) | dict(zip(free, x))
        return jnp.sum(prob_fn(b, inputs, "estimate")["model"])

    x0 = jnp.asarray([values[n] for n in free])
    g = jax.jit(jax.grad(total))(x0)
    assert np.all(np.isfinite(np.asarray(g)))
    assert {c[0] for c in calls} == {1, 2}
    # The penalty enters each class: class 2 loses twice what class 1 loses per row.
    P = prob_fn(values, inputs, "estimate")
    base_fn, base_inputs = compiled(latent_class("lc_cholesky"), "jax")
    B = base_fn(values, base_inputs, "estimate")
    tf1 = database()["Tf1"].to_numpy(float)
    np.testing.assert_allclose(np.asarray(P["class_2"]), np.asarray(B["class_2"]) - 0.02 * tf1, rtol=1e-12)


def test_data_wrapper():
    inputs = engine.make_inputs(database(), {"a": 1.0}, [], CONTROL)
    data = Data(inputs, np)
    assert data.xp is np and data.n_obs == len(database())
    np.testing.assert_array_equal(data["Tf1"], database()["Tf1"].to_numpy(float))
    assert data[["Tf1", "Tf2"]].shape == (len(database()), 2)
    jdata = Data(inputs, jnp)
    assert isinstance(jdata["Tf1"], jax.Array)


def test_profiled_rejects_a_class_loglik_override():
    class Custom(LatentClassMTUEM):
        def class_loglik(self, s, params, data):
            return super().class_loglik(s, params, data)

    model = Custom.model_validate(latent_class("lc_profiled").model_dump())
    with pytest.raises(ValueError, match="class_loglik"):
        model._compile("numpy")


# ----------------------------------------------------------------------------------------------
# Estimation
# ----------------------------------------------------------------------------------------------
@pytest.mark.slow
def test_lc_cholesky_matches_r_at_r_optimum():
    ref = load_ref("lc_cholesky")
    model = latent_class("lc_cholesky")
    beta = model.beta().set(**ref["estimate"])
    result = model.fit(database(), beta, Control(model_name="lc_cholesky", indiv_id="PeID"), NUMPY)
    assert result.loglik == pytest.approx(ref["maximum"], abs=1e-3)
    se, robse = num(ref["se"]).dropna(), num(ref["robse"]).dropna()
    np.testing.assert_allclose(result.se[se.index], se, rtol=2e-3)
    np.testing.assert_allclose(result.robse[robse.index], robse, rtol=2e-3)
    assert list(result._estimate.classShares.index) == ["class_1", "class_2"]
    assert result.class_shares().sum() == pytest.approx(1)


@pytest.mark.slow
@pytest.mark.parametrize("case", ["lc_profiled", "lc_shared"])
def test_lc_reaches_at_least_r_likelihood(case):
    ref = load_ref(case)
    model = latent_class(case)
    beta = model.beta().set(**ref["beta_start"])
    settings = EstimationSettings(backend="numpy", hessian="none", verbose=False)
    result = model.fit(database(), beta, Control(model_name=case, indiv_id="PeID"), settings)
    assert result.loglik >= ref["maximum"] - 1e-3


@pytest.mark.slow
def test_class_without_free_activities_fits():
    model = LatentClassMTUEM(
        work_times=["Tw"], free_times=["Tf1"], Tc="Tc", Ec="EcI", w="w", covariance="cholesky",
        n_classes=2, participation=[[False], [True]], allocation=["female"],
    )
    beta = model.beta()
    assert sorted(beta.locked) == ["chk21_1", "chk22_1", "th1_1"]
    assert "chk11_1" in beta.free_names
    prob_fn, inputs = compiled(model, "numpy")
    values = beta.values()
    class_1 = np.asarray(prob_fn(values, inputs, "estimate")["class_1"])
    one = mtuem.MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", covariance="cholesky")
    opf, _, _ = one._compile("numpy")
    one_ll = opf({"Theta": 1, "Phi": 1, "thw": 0, "chk11": 10}, inputs, "estimate")["model"]
    np.testing.assert_allclose(class_1, np.asarray(one_ll), rtol=1e-12)
    # Without a participation penalty the work-only class takes everyone (its density covers
    # one equation), so the fit uses the article's penalty and data where a third of the
    # individuals spend no time on Tf1.
    data = database()
    data.loc[data.index % 3 == 0, "Tf1"] = 0.0
    penalized = Participation.model_validate(model.model_dump() | {"allocation": [], "penalization": 1.0})
    start = penalized.beta().set(Phi_1=0.2, thw_1=0.3, chk11_1=6, Phi_2=0.18, thw_2=0.31, th1_2=0.7,
                                 chk11_2=6.4, chk21_2=-5, chk22_2=5)
    result = penalized.fit(data, start, CONTROL, JAX)
    assert result.converged and result._estimate.eigValue < 0
    assert np.isnan(result.se["th1_1"]) and result.estimates["th1_1"] == 0
    assert result.se["chk11_1"] > 0 and result.se["chk22_2"] > 0
    assert result.class_shares()["class_1"] == pytest.approx(1 / 3, abs=0.05)


class Participation(LatentClassMTUEM):
    """The participation penalty of the activity participation article."""

    penalization: float = 0.0

    def class_loglik(self, s, params, data):
        ll = super().class_loglik(s, params, data)
        absent = self.absent_times(s)
        return ll - self.penalization * data[absent].sum(axis=1) if absent else ll


# ----------------------------------------------------------------------------------------------
# Results
# ----------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def panel():
    db = database()
    return pd.concat([db, db]).sort_index(kind="stable").reset_index(drop=True)


@pytest.mark.parametrize("case", ["lc_cholesky", "lc_shared"])
def test_posterior_matches_lc_conditionals(case):
    ref = load_ref(case)
    result = result_at(latent_class(case), ref["estimate"])
    post = result.posterior()
    assert list(post.columns) == ["ID", "class_1", "class_2"]
    expected = engine.lc_conditionals(result._estimate, result._prob_fn, result._inputs)
    np.testing.assert_allclose(post[["class_1", "class_2"]].to_numpy(), expected[["X1", "X2"]].to_numpy())
    np.testing.assert_array_equal(post["ID"], database()["PeID"])


def test_posterior_of_profiled_model():
    ref = load_ref("lc_profiled")
    result = result_at(latent_class("lc_profiled"), ref["estimate"])
    post = result.posterior()
    assert list(post.columns) == ["ID", "class_1", "class_2"] and len(post) == len(database())
    _, _, pf, inputs = build("lc_profiled")
    pf(ref["estimate"], inputs, "estimate")
    expected = mtuem.lc_profiled.get_profiled_sigma()["posterior"]
    np.testing.assert_allclose(post[["class_1", "class_2"]].to_numpy(), expected, rtol=1e-10)
    np.testing.assert_array_equal(post["ID"], database()["PeID"])
    np.testing.assert_allclose(result.class_shares().to_numpy(), ref["profiled"]["shares"], rtol=1e-10)


def test_profiled_covariance_matches_r():
    ref = load_ref("lc_profiled")
    result = result_at(latent_class("lc_profiled"), ref["estimate"])
    prof = result.profiled_covariance()
    eqs = ["Tw", "Tf1", "Ef1"]
    for s in range(2):
        assert list(prof["Sigma"][s].index) == eqs and list(prof["Sigma"][s].columns) == eqs
        np.testing.assert_allclose(prof["Sigma"][s].to_numpy(), np.array(ref["profiled"]["Sigma"][s]), rtol=1e-8)
        sigma = np.array(ref["profiled"]["Sigma"][s])
        sd = np.sqrt(np.diag(sigma))
        np.testing.assert_allclose(prof["sig"][s].to_numpy(), sd, rtol=1e-8)
        np.testing.assert_allclose(prof["rho"][s].to_numpy(), sigma / np.outer(sd, sd), rtol=1e-8)
    assert list(prof["shares"].index) == ["class_1", "class_2"]
    np.testing.assert_allclose(prof["shares"].to_numpy(), ref["profiled"]["shares"], rtol=1e-10)
    # Information criteria counting the profiled-out covariance terms, 2 * 3 * 4 / 2 of them.
    assert prof["k"] == ref["ic"]["k"]
    assert prof["AIC"] == pytest.approx(ref["ic"]["AIC"], rel=1e-8)
    assert prof["BIC"] == pytest.approx(ref["ic"]["BIC"], rel=1e-8)


def test_profiled_covariance_needs_covariance_profiled():
    ref = load_ref("lc_cholesky")
    result = result_at(latent_class("lc_cholesky"), ref["estimate"])
    with pytest.raises(ValueError, match="profiled"):
        result.profiled_covariance()
    with pytest.raises(ValueError, match="single-class"):
        result.residual_covariance()


def test_posterior_of_a_degenerate_profiled_model_raises():
    result = result_at(latent_class("lc_profiled"), load_ref("lc_profiled")["estimate"])
    result.posterior()  # leaves posteriors from a good evaluation behind
    n = len(database())
    result._prob_fn = lambda beta, inputs, functionality: {"model": np.full(n, -np.inf)}
    with pytest.raises(ValueError, match="degenerate"):
        result.posterior()


def test_posterior_one_row_per_individual_on_panel_data(panel):
    ref = load_ref("lc_cholesky")
    model = latent_class("lc_cholesky")
    post = result_at(model, ref["estimate"], data=panel).posterior()
    ids = database()["PeID"].to_numpy()
    assert len(post) == len(ids)
    np.testing.assert_array_equal(post["ID"], ids)
    np.testing.assert_allclose(post[["class_1", "class_2"]].sum(axis=1), 1)
    # Two identical rows per individual: the class likelihood of an individual is squared.
    prob_fn, inputs = compiled(model, "numpy")
    P = prob_fn(ref["estimate"], inputs, "estimate")
    pis = inputs.attach(ref["estimate"]).pi_values.log
    a = np.column_stack([pis[s] + 2 * np.asarray(P[f"class_{s + 1}"]) for s in range(2)])
    expected = np.exp(a - np.logaddexp(a[:, 0], a[:, 1])[:, None])
    np.testing.assert_allclose(post[["class_1", "class_2"]].to_numpy(), expected, rtol=1e-10)


def test_profiled_model_rejects_panel_data(panel):
    model = latent_class("lc_profiled")
    with pytest.raises(ValueError, match="panel"):
        model.fit(panel, model.beta(), CONTROL, NUMPY)


def test_posterior_needs_a_latent_class_model():
    result = result_at(latent_class("lc_cholesky"), load_ref("lc_cholesky")["estimate"])
    single = mtuem.MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w")
    other = Result(single, single.beta(), CONTROL, NUMPY, result._estimate, result._inputs, result._prob_fn)
    with pytest.raises(ValueError, match="latent class"):
        other.posterior()
    with pytest.raises(ValueError, match="latent class"):
        other.class_shares()


@pytest.mark.parametrize("case", LATENT)
def test_values_of_time_use_posterior_weights(case):
    ref = load_ref(case)
    vot = result_at(latent_class(case), ref["estimate"]).values_of_time(standard_errors=False)
    np.testing.assert_allclose(vot["summary"]["mean"], pd.DataFrame(ref["vot"])["mean"], rtol=1e-6)
    assert sorted(set(vot["summary"]["class"])) == [1, 2]


def test_values_of_time_on_panel_data(panel):
    ref = load_ref("lc_cholesky")
    vot = result_at(latent_class("lc_cholesky"), ref["estimate"], data=panel).values_of_time(standard_errors=False)
    assert len(vot["individual"]) == 2 * len(panel)


def test_values_of_time_need_class_specific_work_elasticities():
    model = LatentClassMTUEM.model_validate(latent_class("lc_cholesky").model_dump() | {"shared": ["Phi"]})
    values = {k: v for k, v in load_ref("lc_cholesky")["estimate"].items() if k not in ("Phi_1", "Phi_2")}
    with pytest.raises(ValueError, match="Phi"):
        result_at(model, values | {"Phi": 0.3}).values_of_time()


def test_correlations_per_class():
    result = result_at(hobbies(), hobbies().beta().values())
    cr = result.correlations()
    assert list(cr.columns[:2]) == ["class", "Expression"]
    first = cr[cr["class"] == 1]
    assert list(first["Expression"]) == ["sig_Tw", "sig_Tf1", "sig_Td", "rho_Tw_Tf1", "rho_Tw_Td", "rho_Tf1_Td"]
    np.testing.assert_allclose(first["Value"].iloc[:3], 10)
    assert len(cr[cr["class"] == 2]) == 10


def test_correlations_need_a_cholesky_model():
    result = result_at(latent_class("lc_profiled"), load_ref("lc_profiled")["estimate"])
    with pytest.raises(ValueError, match="cholesky"):
        result.correlations()


def test_save_and_load_round_trip(tmp_path):
    ref = load_ref("lc_cholesky")
    model = latent_class("lc_cholesky")
    result = result_at(model, ref["estimate"])
    files = result.save(tmp_path)
    loaded = mtuem.load(files["model"], data=database())
    assert type(loaded.model) is LatentClassMTUEM and loaded.model == model
    pd.testing.assert_frame_equal(loaded.posterior(), result.posterior())


class Penalty(LatentClassMTUEM):
    penalization: float = 0.0

    def class_loglik(self, s, params, data):
        return super().class_loglik(s, params, data) - self.penalization * data["Tf1"]


def test_load_user_subclass_with_model_class(tmp_path):
    model = Penalty.model_validate(hobbies().model_dump() | {"penalization": 0.1})
    result = result_at(model, model.beta().values())
    path = result.save(tmp_path)["model"]
    loaded = mtuem.load(path, data=database(), model_class=Penalty)
    assert loaded.model.penalization == 0.1
    assert loaded.beta.locked == model.beta().locked
    pd.testing.assert_frame_equal(loaded.posterior(), result.posterior())


def test_load_resolves_a_subclass_whose_module_is_imported(tmp_path):
    model = Penalty.model_validate(hobbies().model_dump() | {"penalization": 0.1})
    path = result_at(model, model.beta().values()).save(tmp_path)["model"]
    assert type(mtuem.load(path).model) is Penalty


def test_normalization_groups_are_per_class():
    groups = latent_class("lc_shared").normalization_groups()
    assert groups["work_elasticities"] == ["Theta_1", "Phi_1", "thw_1", "Theta_2", "Phi_2", "thw_2"]
    assert groups["times_elasticities"] == ["th1_1", "th1_2"]
    assert groups["cholesky"] == CHK


def test_absent_times_rejects_an_unknown_class():
    with pytest.raises(ValueError, match="from 1 to 2"):
        hobbies().absent_times(3)


@pytest.mark.parametrize("backend", ["numpy", "jax"])
def test_predict_gives_class_allocations_and_probabilities(backend):
    model = hobbies()
    values = model.beta().values() | {"th1_1": 0.4, "th3_1": 0.2, "th1_2": 0.3, "th2_2": 0.1, "th3_2": 0.2}
    pred = result_at(model, values, backend=backend).predict()
    assert len(pred) == len(database())
    # Class 1 does not take part in Tf2, so it predicts no Tf2 allocation.
    assert {"Tw_class1", "Tf1_class1", "Td_class1", "prob_class1", "Tf2_class2", "prob_class2"} <= set(pred.columns)
    assert "Tf2_class1" not in pred.columns
    np.testing.assert_allclose(pred["prob_class1"] + pred["prob_class2"], 1)
