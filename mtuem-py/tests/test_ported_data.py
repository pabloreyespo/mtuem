"""Port of tests/testthat/test_data.R."""

import numpy as np
import pandas as pd
import pytest

import mtuem


def test_all_documented_datasets_ship_with_the_package():
    for name in ("maed", "enut.i", "enut.i.raw", "enut.ii", "enut.ii.raw"):
        assert name in mtuem.DATASETS
        assert isinstance(mtuem.load_dataset(name), pd.DataFrame)


def test_maed_carries_the_columns_the_model_examples_rely_on():
    d = mtuem.maed()

    assert isinstance(d, pd.DataFrame)
    assert len(d) > 0
    needed = ["PeID", "Tw", "Tc", "Tf1", "Tf2", "EcI", "w", "Ef1", "Ef2", "ta",
              "female", "older45", "professional", "self_employed"]
    assert set(needed) <= set(d.columns)
    assert not d["PeID"].duplicated().any()


def test_maed_has_no_missing_values():
    assert int(mtuem.maed().isna().sum().sum()) == 0


def test_maed_satisfies_the_sign_restrictions_the_likelihood_assumes():
    d = mtuem.maed()

    assert (d["w"] > 0).all()
    assert (d["Tw"] > 0).all()
    assert (d["Tc"] > 0).all()
    assert (d["Tf1"] >= 0).all()
    assert (d["Tf2"] >= 0).all()
    assert (d["ta"] == 168).all()
    assert (d["Tw"] + d["Tc"] < d["ta"]).all()


def test_the_maed_estimation_subset_has_strictly_positive_committed_expenditure():
    d = mtuem.maed()

    keep = d[d["EcI"] > 0]

    assert len(keep) > 500
    # Disposable income is negative for a minority of the sample; pin the share.
    disposable = keep["w"] * keep["Tw"] - keep["EcI"]
    assert np.mean(disposable <= 0) < 0.1


def test_maed_time_use_exhausts_the_time_budget_so_the_3_equation_time_model_is_unidentified():
    d = mtuem.maed()

    slack = d["ta"] - (d["Tw"] + d["Tf1"] + d["Tf2"] + d["Tc"])

    assert (slack.abs() <= 1.05).all()


@pytest.mark.parametrize("v", ["female", "older45", "professional", "self_employed"])
def test_the_maed_binary_covariates_really_are_binary(v):
    # R also covers factor columns through as.numeric(); parquet columns are numeric here.
    assert set(pd.to_numeric(mtuem.maed()[v]).unique()) <= {0, 1}


@pytest.mark.parametrize("name", ["enut.bigcities", "enut.bigcities.raw"])
def test_enut_bigcities_savings_and_expenditures_come_from_the_same_donors(name):
    d = mtuem.load_dataset(name)
    con = d[d["ing_personal"] > 0]

    slack = d["ing_personal"] - (d["Ef1"] + d["Ef2"] + d["Ec"] + d["savings"])
    rate = con["savings"] / con["ing_personal"]

    assert (slack.abs() <= 0.05).all()
    assert (rate > -1).all()
    assert ((d["total_expenses"] - (d["Ef1"] + d["Ef2"] + d["Ec"])).abs() <= 0.05).all()


@pytest.mark.parametrize("name", ["enut.i", "enut.i.raw", "enut.ii", "enut.ii.raw"])
def test_enut_i_and_ii_jdi_columns_exhaust_income_and_stay_in_the_epf_support(name):
    d = mtuem.load_dataset(name)
    if "Ef1_jdi" not in d.columns:
        pytest.skip(f"{name} predates the _jdi columns: rerun its pipeline and data-raw/update_from_pipelines.R")
    jdi = ["Ef1_jdi", "Ef2_jdi", "Ec_jdi", "savings_jdi"]
    fmnl = ["Ef1_fmnl", "Ef2_fmnl", "Ec_fmnl", "savings_fmnl"]
    assert set(jdi + fmnl) <= set(d.columns)
    con = d[(d["ing_personal"] > 0) & d["savings_jdi"].notna()]

    slack_jdi = con["ing_personal"] - con[jdi].sum(axis=1)
    slack_fmnl = con["ing_personal"] - con[fmnl].sum(axis=1)
    rate = con["savings_jdi"] / con["ing_personal"]

    assert (slack_jdi.abs() <= 0.05).all()
    assert (slack_fmnl.abs() <= 0.1).all()
    assert (rate > -1).all()
