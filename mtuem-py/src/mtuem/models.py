"""Declarative MTUEM models: describe the model, get its parameters, fit it to data."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any, Literal

import numpy as np
import pandas as pd
from pydantic import model_validator

from . import _engine as engine
from . import likelihood, likelihood_jax
from .beta import Beta
from .config import Control, EstimationSettings, _Strict
from .results import Result

ProbFn = Callable[[Mapping, Any, str], dict]
LcParsFn = Callable[[Any, Any], dict]

# Settings keys the model builds itself; ``extra`` may not override them.
_BUILT = frozenset({
    "work_times", "free_times", "free_goods", "goods_cost", "Tc", "Ec", "w", "tau", "optimal_tw",
    "work_elasticities", "times_elasticities", "goods_elasticities", "work_elasticity", "work_satiation",
    "times_satiations", "goods_satiations", "sig", "rho", "cholesky",
})
_OMITTED = ("omitted_times", "omitted_goods")


def _numeric(v):
    """Plain numbers as float (as ``Inputs.attach`` does); traced or array values unchanged."""
    return float(v) if isinstance(v, (int, float, np.number)) else v


class MTUEM(_Strict):
    """A single-class MTUEM, described by its data columns, utility form and error structure.

    ``utility`` picks the likelihood (``mtuem.likelihood.theta_phi``, ``alpha_beta`` or
    ``additive_quadratic``) and ``covariance`` the error structure: ``profiled`` (estimated from
    the residuals), ``diagonal`` (``sig{i}``), ``full`` (``sig{i}`` and ``rho{i}{j}``) or
    ``cholesky`` (``chk{i}{j}``, row-major lower triangle). Equations are numbered 1 = work, then
    the free times, then the free goods (without the omitted ones of the additive quadratic form).

    ``Tc``, ``Ec`` and ``w`` are column names or constants; ``goods_cost`` maps goods to prices
    (constants, or column names except for the additive quadratic form), 1 when missing.
    ``elasticity_covariates`` adds linear terms: ``{"Phi": ["female"]}`` creates ``Phi_female``
    and the elasticity becomes ``Phi + female * Phi_female`` row by row. ``extra`` is passed to the
    likelihood settings unchanged (for example ``omitted_times``).
    """

    work_times: list[str]
    free_times: list[str] = []
    free_goods: list[str] = []
    goods_cost: dict[str, float | str] | None = None
    Tc: str | float
    Ec: str | float
    w: str | float
    tau: float = 168
    utility: Literal["theta_phi", "alpha_beta", "additive_quadratic"] = "theta_phi"
    covariance: Literal["profiled", "diagonal", "full", "cholesky"] = "profiled"
    optimal_tw: bool = True
    elasticity_covariates: dict[str, list[str]] = {}
    extra: dict[str, Any] = {}

    @model_validator(mode="after")
    def _check(self) -> MTUEM:
        if len(self.work_times) != 1:
            raise ValueError("work_times must name exactly one column")
        allocations = self.work_times + self.free_times + self.free_goods
        repeated = sorted({c for c in allocations if allocations.count(c) > 1})
        if repeated:
            raise ValueError(
                f"columns named more than once in work_times, free_times and free_goods: {', '.join(repeated)}"
            )
        unknown = [g for g in self.goods_cost or {} if g not in self.free_goods]
        if unknown:
            raise ValueError(f"goods_cost names goods that are not in free_goods: {', '.join(unknown)}")
        clash = sorted(_BUILT.intersection(self.extra))
        if clash:
            raise ValueError(f"extra must not set keys the model builds: {', '.join(clash)}")
        aq = self.utility == "additive_quadratic"
        if aq:
            if self.elasticity_covariates:
                raise ValueError("elasticity_covariates are not available with utility='additive_quadratic'")
            if any(isinstance(v, str) for v in (self.goods_cost or {}).values()):
                raise ValueError(
                    "goods_cost columns are not available with utility='additive_quadratic'; use constants"
                )
            aq_settings = likelihood._with_defaults(self._allocation_settings(), likelihood._AQ_DEFAULTS)
            likelihood._additive_quadratic_equations(aq_settings)  # checks the omitted equations
        elif any(k in self.extra for k in _OMITTED):
            raise ValueError("omitted_times / omitted_goods apply to utility='additive_quadratic' only")
        if self.utility == "alpha_beta" and not self.optimal_tw:
            raise ValueError("optimal_tw=False is not available with utility='alpha_beta'")
        elasticities = [n for names in self._groups().values() for n in names]
        unknown = [k for k in self.elasticity_covariates if k not in elasticities]
        if unknown:
            raise ValueError(
                f"elasticity_covariates names unknown elasticities: {', '.join(unknown)}; "
                f"this model has {', '.join(elasticities)}"
            )
        fixed = [k for k in self.elasticity_covariates if self._base_defaults()[k][1]]
        if fixed:
            raise ValueError(
                f"elasticity_covariates on {', '.join(fixed)}: it is fixed to normalise the model, "
                "so it takes no covariate terms"
            )
        for k, covs in self.elasticity_covariates.items():
            if len(set(covs)) != len(covs):
                raise ValueError(f"elasticity_covariates['{k}'] names a covariate more than once")
            digits = [c for c in covs if c.isdigit()]
            if digits:
                raise ValueError(
                    f"elasticity_covariates['{k}'] names covariates made only of digits ({', '.join(digits)}): "
                    f"{k}_{digits[0]} would read as the class suffix of {k}; rename the column"
                )
        return self

    # ------------------------------------------------------------------------------------------
    # Structure
    # ------------------------------------------------------------------------------------------
    def _allocation_settings(self) -> dict:
        s = {
            "work_times": list(self.work_times),
            "free_times": list(self.free_times),
            "free_goods": list(self.free_goods),
        }
        s.update({k: self.extra[k] for k in _OMITTED if k in self.extra})
        return s

    def _groups(self) -> dict[str, list[str]]:
        """Likelihood settings key -> parameter base names of the utility."""
        th = [f"th{k}" for k in range(1, len(self.free_times) + 1)]
        ph = [f"ph{k}" for k in range(1, len(self.free_goods) + 1)]
        if self.utility == "theta_phi":
            work = {"work_elasticities": ["Theta", "Phi", "thw"]}
        elif self.utility == "alpha_beta":
            work = {"work_elasticities": ["alpha", "beta"]}
        else:
            return {
                "work_elasticity": ["theta_w"],
                "work_satiation": ["beta_w"],
                "times_elasticities": th,
                "times_satiations": [f"beta{k}" for k in range(1, len(self.free_times) + 1)],
                "goods_elasticities": ph,
                "goods_satiations": [f"eta{k}" for k in range(1, len(self.free_goods) + 1)],
            }
        return work | {"times_elasticities": th, "goods_elasticities": ph}

    def _equations(self) -> list[str]:
        """The allocation columns that enter the likelihood, in equation order."""
        omitted = [c for k in _OMITTED for c in self.extra.get(k, [])]
        return [c for c in self.work_times + self.free_times + self.free_goods if c not in omitted]

    def _covariance_names(self) -> dict[str, list[str]]:
        M = len(self._equations())
        sig = [f"sig{i}" for i in range(1, M + 1)]
        rho = [f"rho{i}{j}" for i in range(1, M + 1) for j in range(i + 1, M + 1)]
        chol = [f"chk{i}{j}" for i in range(1, M + 1) for j in range(1, i + 1)]
        return {
            "sig": sig if self.covariance in ("diagonal", "full") else [],
            "rho": rho if self.covariance == "full" else [],
            "cholesky": chol if self.covariance == "cholesky" else [],
        }

    def _base_defaults(self) -> dict[str, tuple[float, bool]]:
        """Utility parameters -> (start value, fixed), in parameter order."""
        g = self._groups()
        out: dict[str, tuple[float, bool]] = {}
        if self.utility == "theta_phi":
            out |= {"Theta": (1.0, True), "Phi": (1.0, False), "thw": (0.0, False)}
        elif self.utility == "alpha_beta":
            out |= {"alpha": (0.25, False), "beta": (0.25, False)}
        else:
            # Location: theta_w = 0. Scale: the first goods satiation (else the first times
            # satiation, else beta_w) is 1, as in the R additive quadratic examples.
            scale = (g["goods_satiations"] or g["times_satiations"] or ["beta_w"])[0]
            out |= {"theta_w": (0.0, True), "beta_w": (1.0, scale == "beta_w")}
            for th, bt in zip(g["times_elasticities"], g["times_satiations"]):
                out |= {th: (0.3, False), bt: (1.0, bt == scale)}
            for ph, eta in zip(g["goods_elasticities"], g["goods_satiations"]):
                out |= {ph: (0.3, False), eta: (1.0, eta == scale)}
            return out
        out |= {n: (0.3, False) for n in g["times_elasticities"] + g["goods_elasticities"]}
        return out

    def _defaults(self) -> dict[str, tuple[float, bool]]:
        out: dict[str, tuple[float, bool]] = {}
        for name, default in self._base_defaults().items():
            out[name] = default
            out |= {f"{name}_{c}": (0.0, False) for c in self.elasticity_covariates.get(name, [])}
        cov = self._covariance_names()
        out |= {n: (10.0, False) for n in cov["sig"]}
        out |= {n: (0.0, False) for n in cov["rho"]}
        diag = {f"chk{i}{i}" for i in range(1, len(self._equations()) + 1)}
        out |= {n: (10.0 if n in diag else 0.0, False) for n in cov["cholesky"]}
        return out

    # ------------------------------------------------------------------------------------------
    # Public description
    # ------------------------------------------------------------------------------------------
    def beta(self) -> Beta:
        """A new ``Beta`` with every parameter of the model and its default start value.

        Theta-Phi fixes ``Theta`` at 1; alpha-beta fixes nothing; the additive quadratic form
        fixes ``theta_w`` at 0 and the first goods satiation at 1. Elasticities start at 0.3,
        covariate terms and correlations at 0, ``sig`` and the Cholesky diagonal at 10.
        """
        defaults = self._defaults()
        values = {k: v for k, (v, _) in defaults.items()}
        return Beta.from_dict(values, fixed=[k for k, (_, fixed) in defaults.items() if fixed])

    def columns(self) -> list[str]:
        """Every data column the model reads."""
        cols = self.work_times + self.free_times + self.free_goods
        cols += [v for v in (self.Tc, self.Ec, self.w) if isinstance(v, str)]
        cols += [v for v in (self.goods_cost or {}).values() if isinstance(v, str)]
        cols += [c for covs in self.elasticity_covariates.values() for c in covs]
        return list(dict.fromkeys(cols))

    def normalization_groups(self) -> dict[str, list[str]]:
        """Parameter names by role, as ``MultiStart.normalization`` needs them.

        Covariate terms are left out. The multistart normalizations are defined for the
        Theta-Phi and alpha-beta forms; the additive quadratic form only lists its covariance
        parameters.
        """
        g = self._groups()
        aq = self.utility == "additive_quadratic"
        return {
            "work_elasticities": [] if aq else list(g["work_elasticities"]),
            "times_elasticities": [] if aq else list(g["times_elasticities"]),
            "goods_elasticities": [] if aq else list(g["goods_elasticities"]),
            **self._covariance_names(),
        }

    # ------------------------------------------------------------------------------------------
    # Likelihood
    # ------------------------------------------------------------------------------------------
    def _settings(self, beta: Mapping, inputs) -> dict:
        """Likelihood settings at the parameter values ``beta``."""
        n = inputs.nObs

        def column(name: str) -> np.ndarray:
            return inputs.matrix([name])[:, 0]

        def context(v):
            return column(v) if isinstance(v, str) else np.full(n, float(v))

        def value(name: str):
            v = _numeric(beta[name])
            for c in self.elasticity_covariates.get(name, []):
                v = v + _numeric(beta[f"{name}_{c}"]) * column(c)
            return v

        costs = self.goods_cost or {}
        s = self._allocation_settings() | {
            "goods_cost": {g: context(costs[g]) if isinstance(costs.get(g), str) else float(costs.get(g, 1.0))
                           for g in self.free_goods},
            "Tc": context(self.Tc), "Ec": context(self.Ec), "w": context(self.w), "tau": self.tau,
        }
        if self.utility != "alpha_beta":
            s["optimal_tw"] = self.optimal_tw
        s |= {key: {n: value(n) for n in names} for key, names in self._groups().items()}
        s |= {key: [_numeric(beta[n]) for n in names] for key, names in self._covariance_names().items() if names}
        s |= {k: v for k, v in self.extra.items() if k not in _OMITTED}
        return s

    def _compile(self, backend: str) -> tuple[ProbFn, LcParsFn | None, bool]:
        """(probabilities function, lc_pars, autodiff) of the model for the engine.

        On ``backend="jax"`` the JAX likelihood is used (autodiff), unless the form has no JAX
        version: then the NumPy one, with a ``FallbackWarning``.
        """
        if backend not in ("jax", "numpy"):
            raise ValueError(f"backend must be 'jax' or 'numpy', not {backend!r}")
        autodiff = backend == "jax" and likelihood_jax.require_jax_form(self.utility)
        loglik = getattr(likelihood_jax if autodiff else likelihood, self.utility)
        settings = self._settings

        def prob_fn(beta: Mapping, inputs, functionality: str = "estimate") -> dict:
            return {"model": loglik(settings(beta, inputs), functionality, inputs)}

        return prob_fn, None, autodiff

    # ------------------------------------------------------------------------------------------
    # Estimation
    # ------------------------------------------------------------------------------------------
    def _validate_data(self, data: pd.DataFrame, control: Control | None = None) -> None:
        """Raise ``ValueError`` naming the columns that are missing, non-numeric or not finite."""
        cols = self.columns()
        missing = [c for c in cols if c not in data.columns]
        if control is not None and control.indiv_id and control.indiv_id not in data.columns:
            missing.append(control.indiv_id)
        if missing:
            raise ValueError(f"the data lack the columns: {', '.join(missing)}")
        if control is not None and control.indiv_id and data[control.indiv_id].isna().any():
            raise ValueError(f"the individual id column {control.indiv_id} has missing values")
        non_numeric = [
            c for c in cols
            if not (pd.api.types.is_numeric_dtype(data[c]) or pd.api.types.is_bool_dtype(data[c]))
        ]
        if non_numeric:
            raise ValueError(f"the model needs numeric columns, these are not: {', '.join(non_numeric)}")
        not_finite = [c for c in cols if not np.all(np.isfinite(data[c].to_numpy(dtype=float, na_value=np.nan)))]
        if not_finite:
            raise ValueError(f"columns with NaN or infinite values: {', '.join(not_finite)}")

    def fit(
        self,
        data: pd.DataFrame,
        beta: Beta,
        control: Control,
        settings: EstimationSettings | None = None,
    ) -> Result:
        """Estimate the model on ``data`` from the start values ``beta``; returns a new ``Result``.

        ``beta`` must hold exactly the parameters of ``self.beta()``; its fixed flags and bounds
        apply. ``settings`` defaults to ``EstimationSettings()``. The model and ``beta`` are not
        changed, so the model can be fitted again (for example over a grid of start values).
        """
        settings = EstimationSettings() if settings is None else settings
        self._check_beta(beta)
        if settings.method == "em":
            raise NotImplementedError("method='em' is not available for this model yet")
        if settings.multistart is not None:
            raise NotImplementedError("multistart is not available for this model yet")
        self._validate_data(data, control)
        beta = beta.model_copy(deep=True)
        silent = not settings.verbose
        prob_fn, lc_pars, autodiff = self._compile(settings.backend)
        inputs = engine.make_inputs(
            data, beta.values(), beta.fixed_names, control, lc_pars, autodiff=autodiff, silent=silent
        )
        est = engine.estimate(
            beta.values(),
            beta.fixed_names,
            prob_fn,
            inputs,
            max_iterations=settings.max_iterations,
            max_function_evals=settings.max_function_evals,
            hessian=settings.resolved_hessian(),
            saddle_escapes=settings.saddle_escapes,
            silent=silent,
            bounds=beta.bounds(),
        )
        return Result(self, beta, control, settings, est, inputs, prob_fn)

    def _check_beta(self, beta: Beta) -> None:
        expected = self.beta().names
        missing = [n for n in expected if n not in beta]
        extra = [n for n in beta.names if n not in expected]
        if missing or extra:
            parts = [f"missing: {', '.join(missing)}" if missing else "", f"extra: {', '.join(extra)}" if extra else ""]
            raise ValueError(
                "beta does not hold the parameters of this model (" + "; ".join(p for p in parts if p) + ")"
            )
