"""Port of tests/testthat/test_values_of_time.R."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

import mtuem
from helper_mtuem import inputs_of


def _db(n, seed):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"Tc": rng.uniform(40, 80, n), "EcI": rng.uniform(50, 150, n), "w": rng.uniform(10, 25, n)})


def _model(**estimate):
    return SimpleNamespace(estimate=pd.Series(estimate, dtype=float))


def _mean(res, quantity, cls=1):
    s = res["summary"]
    return s.loc[(s["quantity"] == quantity) & (s["class"] == cls), "mean"].iloc[0]


def test_single_class_values_of_time_match_the_closed_form():
    n = 50
    db = _db(n, 7)
    we = {"Theta": 1, "Phi": 1.4, "thw": 0.2}
    model = _model(**we)

    res = mtuem.mtuem_values_of_time(model, inputs_of(db), Ec="EcI", standard_errors=False, silent=True)

    Tw = mtuem.get_tw_thph(we, 168, db["Tc"], db["EcI"], db["w"]).ravel()
    expected_vol = np.mean((db["w"] * Tw - db["EcI"]) / (168 - Tw - db["Tc"])) * we["Theta"] / we["Phi"]
    expected_vtaw = np.mean((db["w"] * Tw - db["EcI"]) / Tw) * we["thw"] / we["Phi"]

    assert_allclose(_mean(res, "VoL"), expected_vol)
    assert_allclose(_mean(res, "VTAW"), expected_vtaw)
    assert_allclose(_mean(res, "wage"), db["w"].mean())
    assert len(res["individual"]) == n


def test_vol_minus_vtaw_is_the_wage_rate_per_observation():
    db = _db(40, 8)
    model = _model(Theta=1, Phi=1.1, thw=0.35)

    res = mtuem.mtuem_values_of_time(model, inputs_of(db), Ec="EcI", standard_errors=False, silent=True)

    ind = res["individual"]
    assert_allclose(ind["VoL"].to_numpy() - ind["VTAW"].to_numpy(), db["w"].to_numpy())


def test_class_weights_change_the_reported_means():
    n = 60
    db = _db(n, 9)
    model = _model(Theta_1=1, Phi_1=1.4, thw_1=0.2, Theta_2=1, Phi_2=1.4, thw_2=0.2)

    # Identical parameters in both classes, so only the weights can differ:
    # class 1 keeps the first half of the sample, class 2 the second half.
    W = np.column_stack([np.r_[np.ones(n // 2), np.zeros(n // 2)], np.r_[np.zeros(n // 2), np.ones(n // 2)]])
    res = mtuem.mtuem_values_of_time(
        model, inputs_of(db), nClass=2, weights=W, Ec="EcI", standard_errors=False, silent=True
    )

    assert_allclose(_mean(res, "wage", 1), db["w"].iloc[: n // 2].mean())
    assert_allclose(_mean(res, "wage", 2), db["w"].iloc[n // 2 :].mean())
    assert not np.isclose(_mean(res, "VoL", 1), _mean(res, "VoL", 2))


def test_equal_weights_reproduce_the_unweighted_mean():
    n = 30
    db = _db(n, 10)
    model = _model(Theta=1, Phi=1.2, thw=0.1)

    weighted = mtuem.mtuem_values_of_time(
        model, inputs_of(db), nClass=1, weights=np.full(n, 0.5), Ec="EcI", standard_errors=False, silent=True
    )
    plain = mtuem.mtuem_values_of_time(
        model, inputs_of(db), nClass=1, Ec="EcI", standard_errors=False, silent=True
    )

    assert_allclose(weighted["summary"]["mean"].to_numpy(), plain["summary"]["mean"].to_numpy())


def test_a_missing_class_parameter_is_reported_by_name():
    db = pd.DataFrame({"Tc": [60.0], "EcI": [100.0], "w": [15.0]})
    model = _model(Theta=1, Phi=1, thw=0)

    # The Python port raises a KeyError naming the parameter; R raises a generic error.
    with pytest.raises(KeyError, match="Theta_1"):
        mtuem.mtuem_values_of_time(
            model, inputs_of(db), nClass=2, weights=np.ones((1, 2)), Ec="EcI", standard_errors=False, silent=True
        )


def test_weights_of_the_wrong_shape_are_rejected():
    db = pd.DataFrame({"Tc": [60.0, 62.0], "EcI": [100.0, 110.0], "w": [15.0, 16.0]})
    model = _model(Theta_1=1, Phi_1=1, thw_1=0, Theta_2=1, Phi_2=1, thw_2=0)

    with pytest.raises(ValueError, match="must have 2 rows"):
        mtuem.mtuem_values_of_time(
            model, inputs_of(db), nClass=2, weights=np.ones((3, 2)), Ec="EcI", standard_errors=False, silent=True
        )
