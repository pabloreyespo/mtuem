"""Declarative MTUEM models: describe the model, get its parameters, fit it to data."""

from __future__ import annotations

from collections.abc import Callable, Collection, Mapping
from types import SimpleNamespace
from typing import Any, Literal

import numpy as np
import pandas as pd
from pydantic import model_validator

from . import _engine as engine
from . import likelihood, likelihood_jax
from .beta import Beta
from .config import Control, EstimationSettings, _Strict
from .lc_profiled import mtuem_lc_profiled_likelihood
from .multistart import estimate_once, run_multistart
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


def _fit(model, data: pd.DataFrame, beta: Beta, control: Control, settings: EstimationSettings | None) -> Result:
    """Fit ``model``: check ``beta`` and ``data``, build the inputs once, then run the multistart
    when ``settings.multistart`` is set, else one estimation by ``settings.method``."""
    settings = EstimationSettings() if settings is None else settings
    model._check_beta(beta)
    model._validate_data(data, control)
    prob_fn, lc_pars, autodiff = model._compile(settings.backend)
    ms = settings.multistart
    if lc_pars is None or getattr(model, "covariance", None) == "profiled":
        # EM needs the class components; the profiled likelihood mixes the classes inside.
        what = "a latent class model" if lc_pars is None else "class likelihoods, not covariance='profiled'"
        if settings.method == "em":
            raise ValueError(f"method='em' needs {what}")
        if ms is not None and ms.first_em:
            raise ValueError(f"MultiStart.first_em needs {what}")
    beta = beta.model_copy(deep=True)
    silent = not settings.verbose
    inputs = engine.make_inputs(
        data, beta.values(), beta.fixed_names, control, lc_pars, autodiff=autodiff, silent=silent
    )
    if ms is not None:
        return run_multistart(model, beta, control, settings, inputs, prob_fn)
    est, notes = estimate_once(beta.values(), beta, prob_fn, inputs, settings, silent)
    return Result(model, beta, control, settings, est, inputs, prob_fn, notes=notes)


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

    def _covariance_terms(self, equations: list[int] | None = None) -> dict[str, list[tuple[str, int, int]]]:
        """Covariance parameters as (name, i, j) over the equation numbers ``equations`` (all by default).

        The names keep the numbers of the full model, so a subset (the active equations of a
        latent class) gives the terms among those equations only.
        """
        eqs = list(range(1, len(self._equations()) + 1)) if equations is None else list(equations)
        sig = [(f"sig{i}", i, i) for i in eqs]
        rho = [(f"rho{i}{j}", i, j) for i in eqs for j in eqs if j > i]
        chol = [(f"chk{i}{j}", i, j) for i in eqs for j in eqs if j <= i]
        return {
            "sig": sig if self.covariance in ("diagonal", "full") else [],
            "rho": rho if self.covariance == "full" else [],
            "cholesky": chol if self.covariance == "cholesky" else [],
        }

    def _covariance_names(self, equations: list[int] | None = None) -> dict[str, list[str]]:
        return {key: [t[0] for t in terms] for key, terms in self._covariance_terms(equations).items()}

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
    def _values(self, beta: Mapping, inputs, suffix: str = "", shared: Collection[str] = ()) -> dict:
        """Parameter base name -> value at ``beta``, elasticities with their covariate terms applied.

        A base name ``n`` reads ``beta[n + suffix]`` (``beta[n]`` when ``n`` is in ``shared``).
        Covariate terms are also kept under their own names.
        """
        values = {n: _numeric(beta[n if n in shared else n + suffix]) for n in self._defaults()}
        for name, covs in self.elasticity_covariates.items():
            v = values[name]
            for c in covs:
                v = v + values[f"{name}_{c}"] * inputs.matrix([c])[:, 0]
            values[name] = v
        return values

    def _likelihood_settings(self, values: Mapping, inputs, free_times: list[str] | None = None) -> dict:
        """Likelihood settings at the parameter values ``values`` (base name -> value).

        ``free_times`` keeps only those free times (a latent class that does not take part in
        the others): the other equations, with their elasticities and covariance terms, leave
        the likelihood. Only the Theta-Phi and alpha-beta forms take a subset.
        """
        n = inputs.nObs
        kept = list(self.free_times) if free_times is None else list(free_times)

        def context(v):
            return inputs.matrix([v])[:, 0] if isinstance(v, str) else np.full(n, float(v))

        groups = self._groups()
        if self.utility != "additive_quadratic":
            groups["times_elasticities"] = [
                th for th, t in zip(groups["times_elasticities"], self.free_times, strict=True) if t in kept
            ]
        dropped = [t for t in self.free_times if t not in kept]
        active = [i for i, c in enumerate(self._equations(), start=1) if c not in dropped]
        costs = self.goods_cost or {}
        s = self._allocation_settings() | {
            "free_times": [t for t in self.free_times if t in kept],
            "goods_cost": {g: context(costs[g]) if isinstance(costs.get(g), str) else float(costs.get(g, 1.0))
                           for g in self.free_goods},
            "Tc": context(self.Tc), "Ec": context(self.Ec), "w": context(self.w), "tau": self.tau,
        }
        if self.utility != "alpha_beta":
            s["optimal_tw"] = self.optimal_tw
        s |= {key: {n: values[n] for n in names} for key, names in groups.items()}
        s |= {key: [values[n] for n in names] for key, names in self._covariance_names(active).items() if names}
        s |= {k: v for k, v in self.extra.items() if k not in _OMITTED}
        return s

    def _settings(self, beta: Mapping, inputs) -> dict:
        """Likelihood settings at the parameter values ``beta``."""
        return self._likelihood_settings(self._values(beta, inputs), inputs)

    def _backend(self, backend: str) -> bool:
        """True when ``backend`` runs this utility form on JAX (autodiff); warns on a fallback."""
        if backend not in ("jax", "numpy"):
            raise ValueError(f"backend must be 'jax' or 'numpy', not {backend!r}")
        return backend == "jax" and likelihood_jax.require_jax_form(self.utility)

    def _compile(self, backend: str) -> tuple[ProbFn, LcParsFn | None, bool]:
        """(probabilities function, lc_pars, autodiff) of the model for the engine.

        On ``backend="jax"`` the JAX likelihood is used (autodiff), unless the form has no JAX
        version: then the NumPy one, with a ``FallbackWarning``.
        """
        autodiff = self._backend(backend)
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
        apply. ``settings`` defaults to ``EstimationSettings()``: ``method="em"`` (latent class
        models) runs EM then BFGS, with a BFGS fallback recorded in ``Result.notes``;
        ``multistart`` estimates from random starts and keeps the best (``Result.starts``). The
        model and ``beta`` are not changed, so the model can be fitted again.
        """
        return _fit(self, data, beta, control, settings)

    def _check_beta(self, beta: Beta) -> None:
        expected = self.beta().names
        missing = [n for n in expected if n not in beta]
        extra = [n for n in beta.names if n not in expected]
        if missing or extra:
            parts = [f"missing: {', '.join(missing)}" if missing else "", f"extra: {', '.join(extra)}" if extra else ""]
            raise ValueError(
                "beta does not hold the parameters of this model (" + "; ".join(p for p in parts if p) + ")"
            )


class Data:
    """The data as a class likelihood reads it.

    ``data["Tf1"]`` is a column and ``data[["Tf1", "Tf2"]]`` an N x k matrix, both as arrays of
    ``xp``: ``jax.numpy`` on the JAX backend, ``numpy`` otherwise. ``inputs`` holds the engine
    inputs that the likelihood functions of ``mtuem.likelihood`` and ``mtuem.likelihood_jax``
    take.
    """

    def __init__(self, inputs, xp):
        self.inputs = inputs
        self.xp = xp

    def __getitem__(self, key: str | list[str]):
        if isinstance(key, str):
            return self.xp.asarray(self.inputs.matrix([key])[:, 0])
        return self.xp.asarray(self.inputs.matrix(list(key)))

    @property
    def n_obs(self) -> int:
        return int(self.inputs.nObs)


class LatentClassMTUEM(MTUEM):
    """A latent class MTUEM: ``n_classes`` classes, each an ``MTUEM`` with its own parameters.

    Parameters of class ``s`` carry the suffix ``_s`` (``Theta_1``, ``chk21_2``), except the
    base names listed in ``shared``, which every class uses unsuffixed. ``participation`` has one
    row per class over ``free_times`` (all True by default): a class that does not take part in
    a free time has no equation for it, so its elasticity and its covariance terms leave that
    class's likelihood. Those parameters stay in ``beta()``, fixed at 0 and locked, and the
    covariance terms keep the equation numbers of the full model (``chk43_1``), as the R scripts
    of the activity participation article do.

    Class allocation is a logit with utility ``asc_s + sum_c c * x_{c.lower()}_s`` over the
    ``allocation`` columns; the ``reference_class`` has its ``asc`` and ``x_*`` fixed at 0.
    With ``covariance="profiled"`` the class covariances are profiled out of the likelihood
    (Theta-Phi form, the same equations in every class).

    Override ``class_loglik`` to change the likelihood of a class.
    """

    n_classes: int
    participation: list[list[bool]] | None = None
    allocation: list[str] = []
    reference_class: int = 1
    shared: list[str] = []

    @model_validator(mode="after")
    def _check_classes(self) -> LatentClassMTUEM:
        S = self.n_classes
        if S < 2:
            raise ValueError(f"n_classes must be at least 2, not {S}; use MTUEM for one class")
        if self.participation is not None:
            if len(self.participation) != S:
                raise ValueError(f"participation needs one row per class ({S}), not {len(self.participation)}")
            K = len(self.free_times)
            if any(len(row) != K for row in self.participation):
                raise ValueError(f"participation rows need one entry per free time ({K}: {', '.join(self.free_times)})")
        if not 1 <= self.reference_class <= S:
            raise ValueError(f"reference_class must be a class number from 1 to {S}, not {self.reference_class}")
        if len(set(self.allocation)) != len(self.allocation):
            raise ValueError("allocation names a covariate more than once")
        lowered: dict[str, str] = {}
        for c in self.allocation:
            if c.lower() in lowered:
                raise ValueError(
                    f"allocation covariates {lowered[c.lower()]} and {c} both give the parameter x_{c.lower()}_s"
                )
            lowered[c.lower()] = c
        absent = any(self.absent_times(s) for s in range(1, S + 1))
        if self.covariance == "profiled":
            if self.utility != "theta_phi":
                raise ValueError("covariance='profiled' latent class models need utility='theta_phi'")
            if absent:
                raise ValueError(
                    "classes that do not take part in every free time are not available with "
                    "covariance='profiled': the profiled latent class likelihood needs the same equations "
                    "in every class"
                )
        if absent and self.utility == "additive_quadratic":
            raise ValueError("classes that do not take part in every free time are not available with "
                             "utility='additive_quadratic'")
        names = self._defaults()
        if len(set(self.shared)) != len(self.shared):
            raise ValueError("shared names a parameter more than once")
        unknown = [n for n in self.shared if n not in names]
        if unknown:
            raise ValueError(
                f"shared names unknown parameters: {', '.join(unknown)}; the class parameters are {', '.join(names)}"
            )
        for s in range(1, S + 1):
            reasons = self._absent_names(s)
            for n in self.shared:
                if n in reasons:
                    raise ValueError(f"shared names {n}, but {reasons[n]}")
        return self

    # ------------------------------------------------------------------------------------------
    # Structure
    # ------------------------------------------------------------------------------------------
    def _classes(self) -> range:
        return range(1, self.n_classes + 1)

    def absent_times(self, s: int) -> list[str]:
        """The free times class ``s`` (1-based) does not take part in."""
        if s not in self._classes():
            raise ValueError(f"class must be a number from 1 to {self.n_classes}, not {s}")
        if self.participation is None:
            return []
        return [t for t, takes_part in zip(self.free_times, self.participation[s - 1]) if not takes_part]

    def _active_times(self, s: int) -> list[str]:
        absent = self.absent_times(s)
        return [t for t in self.free_times if t not in absent]

    def _absent_names(self, s: int) -> dict[str, str]:
        """Base names of the parameters class ``s`` does not use -> the reason."""
        absent = self.absent_times(s)
        if not absent:
            return {}
        out: dict[str, str] = {}
        for th, t in zip(self._groups()["times_elasticities"], self.free_times, strict=True):
            if t in absent:
                out[th] = f"class {s} does not take part in {t}"
                out |= {f"{th}_{c}": out[th] for c in self.elasticity_covariates.get(th, [])}
        eqs = self._equations()
        for terms in self._covariance_terms().values():
            for name, i, j in terms:
                hit = [eqs[k - 1] for k in sorted((i, j)) if eqs[k - 1] in absent]
                if hit:
                    out[name] = f"class {s} does not take part in {hit[0]}"
        return out

    def _allocation_names(self, s: int) -> list[str]:
        return [f"x_{c.lower()}_{s}" for c in self.allocation]

    def _name(self, base: str, s: int) -> str:
        return base if base in self.shared else f"{base}_{s}"

    # ------------------------------------------------------------------------------------------
    # Public description
    # ------------------------------------------------------------------------------------------
    def beta(self) -> Beta:
        """A new ``Beta`` with every parameter of the model and its default start value.

        Each class gets the parameters and defaults of ``MTUEM.beta()`` with the suffix ``_s``
        (``Theta_s`` fixed at 1 in the Theta-Phi form), then ``asc_s`` and ``x_{c}_s`` at 0.
        Parameters of free times a class does not take part in are fixed at 0 and locked.
        The ``reference_class`` has its ``asc`` and ``x_*`` fixed (not locked). ``shared``
        parameters come last, without suffix.
        """
        base = self._defaults()
        values: dict[str, float] = {}
        fixed: list[str] = []
        locked: dict[str, str] = {}

        def add(name: str, value: float, is_fixed: bool) -> None:
            values[name] = value
            if is_fixed:
                fixed.append(name)

        for s in self._classes():
            absent = self._absent_names(s)
            for name, (value, is_fixed) in base.items():
                if name in self.shared:
                    continue
                if name in absent:
                    add(f"{name}_{s}", 0.0, True)
                    locked[f"{name}_{s}"] = absent[name]
                else:
                    add(f"{name}_{s}", value, is_fixed)
            reference = s == self.reference_class
            for name in [f"asc_{s}", *self._allocation_names(s)]:
                add(name, 0.0, reference)
        for name, (value, is_fixed) in base.items():
            if name in self.shared:
                add(name, value, is_fixed)
        start = Beta.from_dict(values, fixed=fixed)
        return Beta(params=start.params, locked=locked)

    def columns(self) -> list[str]:
        """Every data column the model reads, allocation covariates included."""
        return list(dict.fromkeys(super().columns() + self.allocation))

    def normalization_groups(self) -> dict[str, list[str]]:
        """Parameter names by role, as ``MTUEM.normalization_groups``, for every class.

        Class parameters come class by class, then the shared ones.
        """
        groups = super().normalization_groups()
        return {
            key: [f"{n}_{s}" for s in self._classes() for n in names if n not in self.shared]
            + [n for n in names if n in self.shared]
            for key, names in groups.items()
        }

    # ------------------------------------------------------------------------------------------
    # Likelihood
    # ------------------------------------------------------------------------------------------
    def class_loglik(self, s: int, params: SimpleNamespace, data: Data):
        """Per-row log-likelihood of class ``s`` (1-based).

        ``params`` holds the values of class ``s`` under the base names (``Theta``, ``th1``,
        ``chk21``, covariate terms such as ``Phi_female``); an elasticity with covariate terms
        already holds its per-row value ``Phi + female * Phi_female``. ``data`` gives the columns
        and the array module ``data.xp``. The default is the standard likelihood over the
        equations of the free times class ``s`` takes part in.

        Override it to customize a class likelihood. On the JAX backend it runs under
        ``jax.jit``: compute with ``data.xp`` and no Python branching on parameter values.
        """
        settings = self._likelihood_settings(vars(params), data.inputs, self._active_times(s))
        loglik = getattr(likelihood if data.xp is np else likelihood_jax, self.utility)
        return loglik(settings, "estimate", data.inputs)

    def _compile(self, backend: str) -> tuple[ProbFn, LcParsFn | None, bool]:
        """(probabilities function, lc_pars, autodiff) of the latent class model for the engine.

        The probabilities function returns ``class_1`` .. ``class_S`` (from ``class_loglik``)
        and their mixture ``model``; with ``covariance="profiled"`` only ``model``.
        """
        autodiff = self._backend(backend)
        profiled = self.covariance == "profiled"
        if profiled and type(self).class_loglik is not LatentClassMTUEM.class_loglik:
            raise ValueError(
                "covariance='profiled' mixes the classes inside its likelihood, so it does not use "
                "class_loglik; choose another covariance to customize class_loglik"
            )
        if autodiff:
            import jax.numpy as xp

            alloc, mix, mix_profiled = likelihood_jax.class_alloc, likelihood_jax.lc, likelihood_jax.lc_profiled
        else:
            xp = np
            alloc, mix, mix_profiled = engine.class_alloc, engine.lc, mtuem_lc_profiled_likelihood
        loglik = getattr(likelihood_jax if autodiff else likelihood, self.utility)
        classes = list(self._classes())

        def lc_pars(b, inputs) -> dict:
            utilities = {}
            for s in classes:
                v = getattr(b, f"asc_{s}")
                for c, name in zip(self.allocation, self._allocation_names(s)):
                    v = v + getattr(b, c) * getattr(b, name)
                utilities[f"class_{s}"] = v
            return {"pi_values": alloc({"utilities": utilities})}

        def prob_fn(beta: Mapping, inputs, functionality: str = "estimate") -> dict:
            pis = inputs.attach(beta).pi_values
            values = {s: self._values(beta, inputs, f"_{s}", self.shared) for s in classes}
            if profiled:
                settings = [self._likelihood_settings(values[s], inputs) for s in classes]
                return {"model": mix_profiled({"class_settings": settings, "classProb": pis}, functionality, inputs)}
            P = {}
            if functionality == "prediction":
                # Predicted allocations come from the standard likelihood (numpy frames).
                for s in classes:
                    settings = self._likelihood_settings(values[s], inputs, self._active_times(s))
                    P[f"class_{s}"] = loglik(settings, "prediction", inputs)
                P["model"] = engine.lc({"inClassProb": dict(P), "classProb": pis}, inputs, functionality)
                return P
            data = Data(inputs, xp)
            for s in classes:
                P[f"class_{s}"] = self.class_loglik(s, SimpleNamespace(**values[s]), data)
            P["model"] = mix({"inClassProb": dict(P), "classProb": pis}, inputs, functionality)
            return P

        return prob_fn, lc_pars, autodiff

    # ------------------------------------------------------------------------------------------
    # Estimation
    # ------------------------------------------------------------------------------------------
    def _validate_data(self, data: pd.DataFrame, control: Control | None = None) -> None:
        super()._validate_data(data, control)
        indiv_id = control.indiv_id if control is not None else None
        if self.covariance == "profiled" and indiv_id and data[indiv_id].duplicated().any():
            raise ValueError(
                "covariance='profiled' mixes the classes row by row, so it does not take panel data "
                f"(several rows per {indiv_id})"
            )

    def _check_beta(self, beta: Beta) -> None:
        super()._check_beta(beta)
        for name, reason in self.beta().locked.items():
            if not beta[name].fixed:
                raise ValueError(f"{name} must stay fixed: {reason}")
