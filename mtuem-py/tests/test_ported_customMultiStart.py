"""Port of tests/testthat/test_customMultiStart.R, on fit with EstimationSettings.multistart."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import mtuem
from mtuem.config import Control, EstimationSettings, MultiStart
from mtuem.models import MTUEM
from mtuem.results import Result

CONTROL = Control(model_name="test-cms", indiv_id="PeID")


def cms_setup():
    """Single-class, single-equation model on a small maed subset."""
    db = mtuem.maed()
    db = db[db["EcI"] > 0].iloc[:120].reset_index(drop=True)
    model = MTUEM(work_times=["Tw"], Tc="Tc", Ec="EcI", w="w", covariance="diagonal")
    beta = model.beta().set(sig1=8)
    return model, db, beta


def settings(**kw) -> EstimationSettings:
    fit = {k: kw.pop(k) for k in ("max_iterations", "hessian", "verbose") if k in kw}
    multistart = MultiStart(**({"seed": 1} | kw))
    return EstimationSettings(backend="numpy", **({"verbose": False} | fit), multistart=multistart)


def test_multistart_returns_the_best_model_its_log_likelihood_and_every_candidate():
    model, db, beta = cms_setup()

    res = model.fit(db, beta, CONTROL, settings(n=3))

    assert len(res.starts) == 3
    succeeded = [r for r in res.starts if isinstance(r, Result)]
    assert len(succeeded) > 0
    assert res.converged
    assert res.loglik >= max(r.loglik for r in succeeded if r.converged) - 1e-8
    assert any(r.loglik == res.loglik and r.beta == res.beta for r in succeeded)


def test_multistart_samples_candidates_inside_the_requested_box():
    # Pinning the box to a single point makes every candidate identical, so all of
    # them must reach the same log-likelihood.
    model, db, beta = cms_setup()
    pinned = {n: beta[n].value for n in beta.free_names}

    res = model.fit(db, beta, CONTROL, settings(n=2, lower=pinned, upper=pinned))

    lls = [r.loglik for r in res.starts if isinstance(r, Result)]
    assert len(lls) == 2
    assert_allclose(np.ptp(lls), 0, atol=1e-8)
    assert all(r.beta.values() == beta.values() for r in res.starts)


def test_theta_phi_normalization_keeps_every_sampled_candidate_feasible(capsys):
    # A box wide enough to produce negative Phi and negative sigma. Without
    # normalization those candidates give a non-finite log-likelihood and get
    # discarded; with it, all of them must survive the initial screening.
    model, db, beta = cms_setup()
    lower = {n: beta[n].value - 12 for n in beta.free_names}
    upper = {n: beta[n].value + 12 for n in beta.free_names}

    with pytest.raises(RuntimeError) as info:
        model.fit(db, beta, CONTROL, settings(n=5, lower=lower, upper=upper, normalization="theta_phi",
                                              max_iterations=2, hessian="none", non_saddle=False, verbose=True))

    assert "5 candidates available" in capsys.readouterr().out
    assert len(info.value.attempts) == 5


def test_estimation_settings_apply_to_every_candidate():
    # max_iterations 2 and no Hessian: no candidate converges, every attempt is kept.
    model, db, beta = cms_setup()

    with pytest.raises(RuntimeError, match="no multistart candidate") as info:
        model.fit(db, beta, CONTROL, settings(n=2, max_iterations=2, hessian="none", non_saddle=False))

    succeeded = [r for r in info.value.attempts if isinstance(r, Result)]
    assert len(succeeded) > 0
    assert all(r._estimate.hessian is None and int(r._estimate.nIter) <= 2 for r in succeeded)
