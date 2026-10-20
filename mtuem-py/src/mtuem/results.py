"""The result of fitting a model: estimates, report, files and post-estimation."""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from . import _engine as engine
from ._engine.latent import class_shares as _class_shares
from ._engine.output import __version__, _estimates_table
from .beta import Beta
from .config import Control, EstimationSettings
from .lc_profiled import _PROFILE_STASH, get_profiled_sigma, profile_ic
from .postestimation import cholesky_delta_method, get_corr, values_of_time

_NO_DATA = "{what} needs the estimation data: pass data= to mtuem.load"


class Result:
    """A fitted model. Returned by ``fit`` and by ``mtuem.load``; the model stays reusable.

    ``beta`` is the start ``Beta`` (with its fixed flags, locks and bounds); ``estimates``,
    ``se`` and ``robse`` cover every parameter (fixed ones without standard errors) and
    ``varcov`` and ``robvarcov`` the free ones. ``starts`` holds the multistart attempts, if any
    (None for a candidate skipped or failed), and ``notes`` what the estimation did beyond the
    plain run (an EM fallback, the multistart candidate kept).
    """

    def __init__(
        self,
        model,
        beta: Beta,
        control: Control,
        settings: EstimationSettings,
        estimate: engine.Estimate,
        inputs: engine.Inputs | None = None,
        prob_fn: Callable | None = None,
        starts: list[Result | None] | None = None,
        notes: list[str] | None = None,
    ):
        # The result keeps its own copy, so later edits of the model do not change it.
        self.model = model.model_copy(deep=True)
        self.beta = beta
        self.control = control
        self.settings = settings
        self.starts = starts
        self.notes = list(notes or [])
        self._estimate = estimate
        self._inputs = inputs
        self._prob_fn = prob_fn

    # ------------------------------------------------------------------------------------------
    # Estimates
    # ------------------------------------------------------------------------------------------
    @property
    def estimates(self) -> pd.Series:
        return self._estimate.estimate

    @property
    def se(self) -> pd.Series:
        return self._estimate.se

    @property
    def robse(self) -> pd.Series:
        return self._estimate.robse

    @property
    def varcov(self) -> pd.DataFrame:
        return self._estimate.varcov

    @property
    def robvarcov(self) -> pd.DataFrame:
        return self._estimate.robvarcov

    @property
    def loglik(self) -> float:
        return float(self._estimate.maximum)

    @property
    def loglik_start(self) -> float:
        return float(self._estimate.LLStart)

    @property
    def converged(self) -> bool:
        return bool(self._estimate.successfulEstimation)

    @property
    def n_obs(self) -> int:
        return int(self._estimate.nObs)

    @property
    def n_indivs(self) -> int:
        return int(self._estimate.nIndivs)

    @property
    def n_free(self) -> int:
        return int(self._estimate.nFreeParams)

    @property
    def aic(self) -> float:
        return float(self._estimate.AIC)

    @property
    def bic(self) -> float:
        return float(self._estimate.BIC)

    def __repr__(self) -> str:
        return (
            f"<Result {type(self.model).__name__} '{self.control.model_name}': loglik={self.loglik:.4f}, "
            f"n_free={self.n_free}, converged={self.converged}>"
        )

    # ------------------------------------------------------------------------------------------
    # Report and tables
    # ------------------------------------------------------------------------------------------
    def summary_text(self, print_covar: bool = False, print_corr: bool = False) -> str:
        """The estimation report as text, followed by the notes of the estimation."""
        text = engine.output_text(self._estimate, print_covar, print_corr)
        if self.notes:
            text += "\nNotes:\n" + "".join(f"- {note}\n" for note in self.notes)
        return text

    def summary(self, print_covar: bool = False, print_corr: bool = False) -> None:
        """Print the estimation report."""
        print(self.summary_text(print_covar, print_corr), end="")

    def to_frame(self) -> pd.DataFrame:
        """The estimates table: estimates, classical and robust standard errors and t-ratios."""
        table = _estimates_table(self._estimate)
        table.columns = ["Estimate", "Std.err.", "t-ratio(0)", "Rob.std.err.", "Rob.t-ratio(0)"]
        return table

    # ------------------------------------------------------------------------------------------
    # Post-estimation
    # ------------------------------------------------------------------------------------------
    def _require_data(self, what: str) -> engine.Inputs:
        if self._inputs is None:
            if getattr(self.model, "detached", False):
                raise ValueError(f"{what} needs the model and the estimation data: pass model= and data= to mtuem.load")
            raise ValueError(_NO_DATA.format(what=what))
        return self._inputs

    def delta_method(self, expressions: dict[str, Callable[[dict], float]], varcov: str = "robust") -> pd.DataFrame:
        """Delta method standard errors of functions of the estimates.

        Each expression receives every parameter value in a dict and returns a number.
        ``varcov`` is ``"robust"`` or ``"classical"``.
        """
        if varcov not in ("robust", "classical"):
            raise ValueError(f"varcov must be 'robust' or 'classical', not {varcov!r}")
        return engine.delta_method(self._estimate, expressions, varcov=varcov, silent=True)

    def values_of_time(self, **kwargs) -> dict:
        """Values of leisure and of time assigned to work (Theta-Phi models).

        The budgets and the wage come from the model; keyword arguments go to
        ``mtuem.postestimation.values_of_time`` (for example ``conf_level`` or ``weights``).
        """
        if getattr(self.model, "is_custom", False):
            # No model to take the budgets from: the caller names the columns.
            missing = [k for k in ("Tc", "Ec", "w") if not isinstance(kwargs.get(k), str)]
            if missing:
                raise ValueError(
                    "values_of_time of a custom model needs the column names Tc, Ec and w as keyword "
                    "arguments (the model does not say which columns they are)"
                )
            if self._n_classes():
                raise ValueError(
                    "values_of_time is not available for a custom latent class result; call "
                    "mtuem.postestimation.values_of_time with the class weights instead"
                )
            inputs = self._require_data("values_of_time")
            absent = [kwargs[k] for k in ("Tc", "Ec", "w") if kwargs[k] not in inputs.database.columns]
            if absent:
                raise ValueError(f"values_of_time: the data lack the columns {', '.join(absent)}")
            return values_of_time(self._estimate, inputs, self._prob_fn, **({"silent": True} | kwargs))
        if getattr(self.model, "utility", None) != "theta_phi":
            raise ValueError("values_of_time needs a model with utility='theta_phi'")
        if self.model.elasticity_covariates:
            raise ValueError("values_of_time does not take elasticity_covariates into account")
        inputs = self._require_data("values_of_time")
        db = inputs.database
        args: dict[str, Any] = {"tau": self.model.tau, "silent": True}
        n_classes = self._n_classes()
        if n_classes:
            shared = [n for n in ("Theta", "Phi", "thw") if n in self.model.shared]
            if shared:
                raise ValueError(
                    f"values_of_time needs class-specific Theta, Phi and thw; this model shares {', '.join(shared)}"
                )
            post = self.posterior()[[f"class_{s}" for s in range(1, n_classes + 1)]].to_numpy()
            # Posterior weights of each individual on each of its rows.
            weights = post if len(post) == inputs.nObs else post[inputs.indiv_index]
            args |= {"nClass": n_classes, "weights": weights}
        constants = {}
        for key, name in (("Tc", "Tc"), ("Ec", "Ec"), ("w", "w")):
            v = getattr(self.model, name)
            if isinstance(v, str):
                args[key] = v
            else:
                column = f"__mtuem_{key}"
                constants[column] = float(v)
                args[key] = column
        if constants:
            inputs = SimpleNamespace(database=db.assign(**constants))
        return values_of_time(self._estimate, inputs, self._prob_fn, **(args | kwargs))

    def predict(self) -> pd.DataFrame:
        """Predicted allocations at the estimates, one row per observation.

        The predictions come from the standard likelihood of the model's utility form, also for
        a ``LatentClassMTUEM`` subclass that overrides ``class_loglik``.
        """
        inputs = self._require_data("predict")
        return engine.prediction(self._estimate, self._prob_fn, inputs, silent=True)

    def correlations(self) -> pd.DataFrame:
        """Error standard deviations and correlations implied by a Cholesky covariance, with
        delta method standard errors (robust)."""
        if getattr(self.model, "covariance", None) != "cholesky":
            raise ValueError("correlations needs a model with covariance='cholesky'")
        model = self.model
        equations = model._equations()
        if not self._n_classes():
            names = model._covariance_names()["cholesky"]
            return cholesky_delta_method(self._estimate, names, equations, {"silent": True})
        # One block per class, over the equations the class takes part in.
        frames = []
        for s in range(1, model.n_classes + 1):
            absent = model.absent_times(s)
            active = [i for i, c in enumerate(equations, start=1) if c not in absent]
            names = [model._name(n, s) for n in model._covariance_names(active)["cholesky"]]
            labels = [equations[i - 1] for i in active]
            frame = cholesky_delta_method(self._estimate, names, labels, {"silent": True})
            frame.insert(0, "class", s)
            frames.append(frame)
        return pd.concat(frames, ignore_index=True)

    def residual_covariance(self) -> dict:
        """Covariance (``covar``), correlations (``corr``) and standard deviations (``sigma``) of
        the observed minus the predicted allocations at the estimates, for a single-class model.

        The empirical error covariance, whatever ``covariance`` the model estimated: with
        ``"profiled"`` it is the covariance the likelihood concentrates out.
        """
        if self._n_classes():
            raise ValueError("residual_covariance needs a single-class model; see profiled_covariance")
        inputs = self._require_data("residual_covariance")
        return get_corr(self._estimate, self._prob_fn, inputs)

    def profiled_covariance(self) -> dict:
        """The class covariances a ``covariance="profiled"`` latent class model concentrates out.

        Evaluated at the estimates: ``Sigma``, ``sig`` and ``rho`` (one per class, labelled by
        equation), the class ``shares``, and ``k``, ``AIC`` and ``BIC`` counting the profiled-out
        covariance terms (``n_classes * M * (M + 1) / 2`` for M equations) as free parameters.
        """
        n_classes = self._n_classes()
        if not n_classes or getattr(self.model, "covariance", None) != "profiled":
            raise ValueError("profiled_covariance needs a latent class model with covariance='profiled'")
        self._profiled_stash(self._require_data("profiled_covariance"))
        prof = get_profiled_sigma()
        eqs = self.model._equations()
        out = {
            "Sigma": [pd.DataFrame(m, index=eqs, columns=eqs) for m in prof["Sigma"]],
            "sig": [pd.Series(v, index=eqs) for v in prof["sig"]],
            "rho": [pd.DataFrame(m, index=eqs, columns=eqs) for m in prof["rho"]],
            "shares": pd.Series(prof["shares"], index=[f"class_{s}" for s in range(1, n_classes + 1)]),
        }
        return out | profile_ic(self._estimate, n_classes, len(eqs))

    def _profiled_stash(self, inputs) -> None:
        """Evaluate the profiled likelihood at the estimates, so its stash holds their covariances
        and posteriors (it is left empty when the evaluation is degenerate)."""
        _PROFILE_STASH.clear()
        self._prob_fn(self.estimates.to_dict(), inputs, "estimate")
        if "W" not in _PROFILE_STASH:
            raise ValueError("the profiled likelihood is degenerate at the estimates: no posterior")

    # ------------------------------------------------------------------------------------------
    # Latent classes
    # ------------------------------------------------------------------------------------------
    def _n_classes(self) -> int:
        """The number of classes of a latent class model, 0 for a single-class model."""
        if getattr(self.model, "is_custom", False):
            # A custom model does not declare its classes: the estimation reports their shares.
            shares = self._estimate.classShares
            return 0 if shares is None else len(shares)
        return int(getattr(self.model, "n_classes", 0) or 0)

    def _require_classes(self, what: str) -> int:
        n_classes = self._n_classes()
        if not n_classes:
            raise ValueError(f"{what} needs a latent class model")
        return n_classes

    def posterior(self) -> pd.DataFrame:
        """Posterior class membership probabilities at the estimates.

        One row per individual (``control.indiv_id``), in order of first appearance in the
        data, with columns ``ID, class_1 .. class_S``.
        """
        n_classes = self._require_classes("posterior")
        inputs = self._require_data("posterior")
        columns = [f"class_{s}" for s in range(1, n_classes + 1)]
        if getattr(self.model, "covariance", None) == "profiled":
            # The profiled likelihood mixes row by row (fit rejects panel data for it) and stores
            # its posteriors when an evaluation succeeds.
            self._profiled_stash(inputs)
            post = _PROFILE_STASH["W"]
            out = pd.DataFrame(post, columns=columns)
            out.insert(0, "ID", inputs.ids)
            return out
        post = engine.lc_conditionals(self._estimate, self._prob_fn, inputs)
        return post.rename(columns={f"X{s}": f"class_{s}" for s in range(1, n_classes + 1)})

    def class_shares(self) -> pd.Series:
        """Mean class allocation probabilities at the estimates (``class_1 .. class_S``)."""
        self._require_classes("class_shares")
        shares = self._estimate.classShares
        if shares is None:
            shares = _class_shares(self.estimates.to_dict(), self._require_data("class_shares"))
        return shares.copy()

    # ------------------------------------------------------------------------------------------
    # Files
    # ------------------------------------------------------------------------------------------
    def save(self, folder: Path | None = None) -> dict[str, Path]:
        """Write the report, estimates and covariance files and ``<model_name>_model.json``.

        ``folder`` defaults to ``control.output_dir``. Returns the paths by kind.
        """
        folder = Path(self.control.output_dir if folder is None else folder)
        name = self.control.model_name
        files = engine.save_files(self._estimate, folder, name)
        files["model"] = folder / f"{name}_model.json"
        files["model"].write_text(json.dumps(self._to_json(), indent=2, allow_nan=False), encoding="utf-8")
        return files

    def _to_json(self) -> dict:
        est = self._estimate
        custom = getattr(self.model, "is_custom", False)
        cls = type(self.model)
        eig = est.hessianEigenValue
        beta = self.beta.model_dump(mode="json")
        for param in beta["params"].values():
            # An infinite bound is no bound.
            for side in ("lower", "upper"):
                if param[side] is not None and not math.isfinite(param[side]):
                    param[side] = None
        # The output folder is a path of this machine: a loaded result saves next to its json.
        control = self.control.model_dump(mode="json") | {"output_dir": "."}
        parts = {
            "model": None if custom else self.model.model_dump(mode="json"),
            "beta": beta,
            "control": control,
            "settings": self.settings.model_dump(mode="json"),
        }
        # JSON has no infinity or NaN, and pydantic writes them as null inside untyped fields:
        # check the Python values so nothing is lost silently.
        sources = {"control": self.control, "settings": self.settings}
        if not custom:
            sources["model"] = self.model
        for key, obj in sources.items():
            if not _finite(obj.model_dump()):
                raise ValueError(f"cannot save the result: the {key} holds a number that is not finite")
        if not _finite([p.value for p in self.beta.params.values()]):
            raise ValueError("cannot save the result: the beta holds a number that is not finite")
        return {
            "mtuem_version": __version__,
            "model_class": "mtuem.models:CustomModel" if custom else f"{cls.__module__}:{cls.__qualname__}",
            **parts,
            "estimates": _series(self.estimates),
            "se": _series(self.se),
            "robse": _series(self.robse),
            "varcov": _matrix(self.varcov),
            "robvarcov": _matrix(self.robvarcov),
            "loglik": _num(self.loglik),
            "loglik_start": _num(self.loglik_start),
            "converged": self.converged,
            "n_obs": self.n_obs,
            "n_indivs": self.n_indivs,
            "notes": list(self.notes),
            # What the report needs beyond the estimates.
            "diagnostics": {
                "code": int(est.code),
                "message": est.message,
                "routine": est.routine,
                "iterations": est.nIter if isinstance(est.nIter, str) else int(est.nIter),
                "hessian_eigenvalues": _eigenvalues(eig),
                "ll_components": {k: _num(v) for k, v in est.LLout.items()},
                "class_shares": None if est.classShares is None else _series(est.classShares),
                "bounds": {k: list(v) for k, v in est.bounds.items()},
                "start_time": est.startTime.isoformat(),
                "time_pre": est.timePre,
                "time_estimation": est.timeEst,
                "time_post": est.timePost,
                "time_taken": est.timeTaken,
            },
        }


# ----------------------------------------------------------------------------------------------
# JSON helpers: values that are not finite (NaN standard errors of fixed parameters) are null
# ----------------------------------------------------------------------------------------------
def _num(x) -> float | None:
    x = float(x)
    return x if math.isfinite(x) else None


def _finite(obj) -> bool:
    """False when a float anywhere in nested dicts, lists or arrays is infinite or NaN."""
    if isinstance(obj, dict):
        return all(_finite(v) for v in obj.values())
    if isinstance(obj, (list, tuple, set)):
        return all(_finite(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return obj.dtype.kind not in "fc" or bool(np.all(np.isfinite(obj)))
    if isinstance(obj, (float, np.floating)):
        return math.isfinite(obj)
    return True


def _series(s: pd.Series) -> dict:
    return {str(k): _num(v) for k, v in s.items()}


def _matrix(m: pd.DataFrame) -> dict:
    return {"names": [str(n) for n in m.index], "values": [[_num(v) for v in row] for row in m.to_numpy()]}


def _eigenvalues(eig) -> dict | None:
    if eig is None:
        return None
    eig = np.asarray(eig)
    out = {"real": [_num(v) for v in np.real(eig)]}
    if np.iscomplexobj(eig):
        out["imag"] = [_num(v) for v in np.imag(eig)]
    return out
