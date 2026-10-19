"""Loading a saved result back from its ``_model.json``."""

from __future__ import annotations

import importlib
import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from . import _engine as engine
from ._engine.estimate import _corr
from .beta import Beta
from .config import Control, EstimationSettings
from .models import MTUEM
from .results import Result


def load(path: Path, data: pd.DataFrame | None = None, model_class: type | None = None) -> Result:
    """Restore a saved result from ``<model_name>_model.json`` (or the folder holding it).

    Without ``data`` the estimates, report, files and delta method are available; values of
    time and predictions need the estimation data. ``model_class`` replaces the model class
    stored in the file, for custom models whose module cannot be imported here.
    """
    path = _model_file(Path(path))
    saved = json.loads(path.read_text(encoding="utf-8"))
    cls = model_class if model_class is not None else _import_class(saved["model_class"], path)
    if not (isinstance(cls, type) and issubclass(cls, MTUEM)):
        raise TypeError(f"the model class must be a subclass of mtuem.MTUEM, not {cls!r}")
    model = cls.model_validate(saved["model"])
    beta = Beta.model_validate(saved["beta"])
    names = model.beta().names
    for what, saved_names in (("saved estimates", list(saved["estimates"])), ("saved beta", beta.names)):
        if sorted(names) != sorted(saved_names):
            raise ValueError(
                f"the saved model has the parameters {', '.join(names)}, "
                f"but the {what} are of {', '.join(saved_names)}"
            )
    # Files of a loaded result go next to its json, not to the folder of the machine that saved it.
    control = Control.model_validate(saved["control"] | {"output_dir": path.parent})
    settings = EstimationSettings.model_validate(saved["settings"])
    estimate = _estimate(saved, beta, control)
    inputs = prob_fn = None
    if data is not None:
        if len(data) != estimate.nObs:
            raise ValueError(f"the data have {len(data)} rows, the saved result was estimated on {estimate.nObs}")
        model._validate_data(data, control)
        prob_fn, lc_pars, autodiff = model._compile(settings.backend)
        inputs = engine.make_inputs(data, beta.values(), beta.fixed_names, control, lc_pars, autodiff=autodiff)
    return Result(model, beta, control, settings, estimate, inputs, prob_fn, notes=saved.get("notes", []))


def _model_file(path: Path) -> Path:
    if not path.is_dir():
        return path
    found = sorted(path.glob("*_model.json"))
    if len(found) != 1:
        raise ValueError(f"{path} holds {len(found)} *_model.json files; pass the path of the one to load")
    return found[0]


def _import_class(spec: str, path: Path) -> type:
    module_name, _, qualname = spec.partition(":")
    try:
        obj = importlib.import_module(module_name)
        for part in qualname.split("."):
            obj = getattr(obj, part)
    except (ImportError, AttributeError) as err:
        raise ImportError(
            f"cannot import the model class {spec} stored in {path}; "
            "pass model_class= to mtuem.load with the class to use"
        ) from err
    return obj


def _float(x) -> float:
    return math.nan if x is None else float(x)


def _series(d: dict) -> pd.Series:
    return pd.Series({k: _float(v) for k, v in d.items()}, dtype=float)


def _matrix(d: dict) -> pd.DataFrame:
    names = d["names"]
    values = np.array([[_float(v) for v in row] for row in d["values"]], dtype=float).reshape(len(names), len(names))
    return pd.DataFrame(values, index=names, columns=names)


def _estimate(saved: dict, beta: Beta, control: Control) -> engine.Estimate:
    diag = saved["diagnostics"]
    varcov, robvarcov = _matrix(saved["varcov"]), _matrix(saved["robvarcov"])
    eig = diag["hessian_eigenvalues"]
    if eig is not None:
        eig = np.array([_float(v) for v in eig["real"]]) + (
            1j * np.array([_float(v) for v in eig["imag"]]) if "imag" in eig else 0
        )
    n_obs, k = int(saved["n_obs"]), len(varcov)
    maximum = _float(saved["loglik"])
    return engine.Estimate(
        estimate=_series(saved["estimates"]),
        beta_start=pd.Series(beta.values(), dtype=float),
        fixed=beta.fixed_names,
        control=control,
        maximum=maximum,
        LLStart=_float(saved["loglik_start"]),
        LL0={name: math.nan for name in diag["ll_components"]},
        LLout={name: _float(v) for name, v in diag["ll_components"].items()},
        se=_series(saved["se"]),
        robse=_series(saved["robse"]),
        varcov=varcov,
        robvarcov=robvarcov,
        corrmat=_corr(varcov),
        robcorrmat=_corr(robvarcov),
        hessianEigenValue=eig,
        eigValue=float(np.max(np.real(eig))) if eig is not None else math.nan,
        code=int(diag["code"]),
        message=diag["message"],
        successfulEstimation=bool(saved["converged"]),
        routine=diag["routine"],
        nIter=diag["iterations"],
        nObs=n_obs,
        nIndivs=int(saved["n_indivs"]),
        nFreeParams=k,
        AIC=-2 * maximum + 2 * k,
        BIC=-2 * maximum + k * math.log(n_obs),
        startTime=datetime.fromisoformat(diag["start_time"]),
        timePre=diag["time_pre"],
        timeEst=diag["time_estimation"],
        timePost=diag["time_post"],
        timeTaken=diag["time_taken"],
        classShares=None if diag["class_shares"] is None else _series(diag["class_shares"]),
    )
