from pathlib import Path

import pytest
from pydantic import ValidationError

from mtuem.config import Control, EstimationSettings, MultiStart


def test_control_defaults():
    c = Control()
    assert c.model_name == "model"
    assert c.description == ""
    assert c.indiv_id is None
    assert c.output_dir == Path(".")


def test_settings_reject_unknown_field():
    with pytest.raises(ValidationError):
        EstimationSettings(maxIterations=5)


@pytest.mark.parametrize("field", ["max_iterations", "max_function_evals", "em_max_iterations"])
def test_settings_reject_nonpositive_iterations(field):
    with pytest.raises(ValidationError):
        EstimationSettings(**{field: 0})


def test_numpy_backend_rejects_analytic_hessian():
    with pytest.raises(ValidationError):
        EstimationSettings(backend="numpy", hessian="analytic")


def test_resolved_hessian_defaults():
    assert EstimationSettings().resolved_hessian() == "analytic"
    assert EstimationSettings(backend="numpy").resolved_hessian() == "numeric"
    assert EstimationSettings(hessian="none").resolved_hessian() == "none"


def test_settings_json_round_trip():
    s = EstimationSettings(method="em", multistart=MultiStart(n=5, seed=1, lower={"a": -1.0}, normalization="theta_phi"))
    assert EstimationSettings.model_validate_json(s.model_dump_json()) == s
