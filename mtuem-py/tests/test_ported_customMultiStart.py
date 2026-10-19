"""Port of tests/testthat/test_customMultiStart.R."""

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

import mtuem
from mtuem import apollo

FAST = {"maxIterations": 2, "hessianRoutine": "none", "estimationRoutine": "bfgs"}


def cms_setup():
    """Single-class, single-equation model on a small maed subset."""
    db = mtuem.maed()
    db = db[db["EcI"] > 0].iloc[:120].reset_index(drop=True)
    beta = {"Theta": 1, "Phi": 1, "thw": 0, "sig1": 8}
    fixed = ["Theta"]
    control = {"modelName": "test-cms", "indivID": "PeID", "noValidation": True, "noDiagnostics": True}
    inputs = apollo.validate_inputs(db, beta, fixed, control, silent=True)

    def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
        b = apollo_inputs.attach(apollo_beta)
        settings = {
            "work_times": ["Tw"],
            "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
            "sig": {"sig1": b.sig1},
            "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
        }
        P = {"model": mtuem.mtuem_likelihood(settings, functionality, apollo_inputs)}
        return apollo.prepare_prob(P, apollo_inputs, functionality)

    return beta, fixed, apollo_probabilities, inputs


def test_custom_multistart_returns_the_best_model_its_log_likelihood_and_every_candidate():
    beta, fixed, pf, inputs = cms_setup()

    res = mtuem.customMultiStart(
        beta, fixed, pf, inputs, customMultistart_settings={"nCandidates": 3}, estimate_settings=FAST,
        verbose=False, seed=1,
    )

    assert set(res) == {"best_model", "best_ll", "models"}
    # Unlike R, only candidates that pass the initial screening get an entry in models.
    assert 0 < len(res["models"]) <= 3
    succeeded = [m for m in res["models"].values() if isinstance(m, apollo.ApolloModel)]
    assert len(succeeded) > 0
    if res["best_model"] is not None:
        assert res["best_ll"] == res["best_model"].maximum
        assert res["best_ll"] >= max(m.maximum for m in succeeded) - 1e-8


def test_custom_multistart_samples_candidates_inside_the_requested_box():
    # Pinning the box to a single point makes every candidate identical, so all of
    # them must reach the same log-likelihood.
    beta, fixed, pf, inputs = cms_setup()
    pinned = pd.Series(beta, dtype=float)

    res = mtuem.customMultiStart(
        beta, fixed, pf, inputs,
        customMultistart_settings={"nCandidates": 2, "apolloBetaMin": pinned, "apolloBetaMax": pinned},
        estimate_settings=FAST, verbose=False, seed=1,
    )

    lls = [m.maximum for m in res["models"].values() if isinstance(m, apollo.ApolloModel)]
    assert len(lls) > 0
    assert_allclose(np.ptp(lls), 0, atol=1e-8)


def test_theta_phi_normalization_keeps_every_sampled_candidate_feasible(capsys):
    # A box wide enough to produce negative Phi and negative sigma. Without
    # normalization those candidates give a non-finite log-likelihood and get
    # discarded; with it, all of them must survive the initial screening.
    beta, fixed, pf, inputs = cms_setup()
    b = pd.Series(beta, dtype=float)

    mtuem.customMultiStart(
        beta, fixed, pf, inputs,
        customMultistart_settings={"nCandidates": 5, "apolloBetaMin": b - 12, "apolloBetaMax": b + 12},
        estimate_settings=FAST, verbose=False,
        normalization={"normalization": "theta_phi", "work_elasticities": ["Theta", "Phi", "thw"], "sig": ["sig1"]},
        seed=1,
    )

    assert "5 candidates available" in capsys.readouterr().out


def test_the_deprecated_estimation_settings_alias_is_still_honoured():
    beta, fixed, pf, inputs = cms_setup()

    res = mtuem.customMultiStart(
        beta, fixed, pf, inputs, customMultistart_settings={"nCandidates": 2},
        estimation_settings=FAST, verbose=False, seed=1,
    )

    succeeded = [m for m in res["models"].values() if isinstance(m, apollo.ApolloModel)]
    assert len(succeeded) > 0
    assert all(m.estimationRoutine == "bfgs" for m in succeeded)
