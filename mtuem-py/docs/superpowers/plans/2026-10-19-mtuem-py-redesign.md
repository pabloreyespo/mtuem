# mtuem-py Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Apollo-shaped API of `mtuem-py` with a declarative, pydantic-configured, statsmodels-style API, with JAX as the default engine, without losing functionality or precision, and prove it by replicating the activity-participation article.

**Architecture:** The numeric core (likelihoods, BFGS optimizer, Hessian, EM, output formatting) moves from `mtuem/apollo.py` into a private `mtuem/_engine/` package with neutral names and keeps its contract: a probabilities function `fn(beta: dict, inputs, functionality) -> dict` plus optional `lc_pars`. Declarative models (`MTUEM`, `LatentClassMTUEM`) and the hand-written `CustomModel` all compile to that contract, so estimation code is shared. `fit` returns a `Result` that owns output, saving and post-estimation.

**Tech Stack:** Python 3.14 (uv), numpy, scipy, pandas, pyarrow, pydantic v2, jax (x64), pytest. R 4.4.3 with apollo 0.3.9, arrow and mtuem for reference fixtures.

**Spec:** `docs/superpowers/specs/2026-10-19-mtuem-py-redesign-design.md`

## Global Constraints

- Python via uv only, 3.14 pinned; run everything with `uv run`.
- Dependencies: numpy, scipy, pandas, pyarrow, pydantic, jax (with `jax_enable_x64`). JAX is required, not an extra.
- No public function, class, module, file name or printed header contains "apollo" (case-insensitive). Docs may cite Apollo (Hess and Palma) as the method reproduced.
- Plain names under the package (`mtuem.MTUEM`, `mtuem.load`), no `mtuem_` prefix on new public API.
- The parameter object is `Beta`; parameter names follow `name_class` (`Theta_1`, `chk21_2`, `x_female_2`); single-class models have no suffix.
- One optimizer: the engine's BFGS (BHHH seeding, saddle escape). No `bgw`, no `estimationRoutine`.
- Tests first; 80% coverage minimum. Existing R parity tolerances are kept: LL `rel=1e-10, abs=1e-8` at fixed points, maximum `abs=1e-3`, estimates `rtol=1e-4, atol=1e-4`, SEs `rtol=2e-3` (NumPy backend); JAX backend SEs `rtol=5e-3`.
- No em dashes, no emojis in code, comments or docs.
- Slow tests carry `@pytest.mark.slow` and are skipped unless `--runslow` is passed.

## Review Focus

1. A DataFrame with NaN or non-numeric values in a column the model uses: `fit` raises `ValueError` naming the columns (Task 7).
2. A latent class whose participation pattern excludes every free activity: the class is a 1-equation model (work only) and fits (Task 9).
3. `beta.free("th2_1")` or `beta.set(th2_1=...)` on a parameter of an activity the class does not take part in: `ValueError` naming the class and activity (Task 9).
4. Panel data (several rows per `indiv_id`): `result.posterior()` has one row per individual, in order of first appearance (Task 9).
5. `mtuem.load` of a result whose model class is not importable (for example the example subclass): `ImportError` telling the user to pass `model_class=` (Task 8).

---

### Task 1: Branch, dependencies and configuration objects

**Files:**
- Modify: `pyproject.toml` (move jax to dependencies, add pydantic, pytest marker `slow`)
- Create: `tests/conftest.py` (`--runslow` option)
- Create: `src/mtuem/config.py`
- Test: `tests/test_config.py`

**Interfaces:**
- Produces: `Control(model_name: str = "model", description: str = "", indiv_id: str | None = None, output_dir: Path = Path("."))`; `MultiStart(n: int = 100, seed: int | None = None, lower: dict[str, float] | None = None, upper: dict[str, float] | None = None, first_em: bool = False, em_iter_max: int = 5, non_saddle: bool = True, normalization: Literal["theta_phi", "alpha_beta"] | None = None)`; `EstimationSettings(method: Literal["bfgs", "em"] = "bfgs", fallback: Literal["bfgs"] | None = "bfgs", hessian: Literal["analytic", "numeric", "none"] | None = None, max_iterations: int = 200, max_function_evals: int = 100000, em_max_iterations: int = 100, em_tolerance: float = 1e-5, saddle_escapes: int = 3, backend: Literal["jax", "numpy"] = "jax", multistart: MultiStart | None = None, verbose: bool = True)`. All are pydantic `BaseModel` with `extra="forbid"`; `hessian=None` resolves to `"analytic"` on jax and `"numeric"` on numpy via `EstimationSettings.resolved_hessian() -> str`; `"analytic"` with `backend="numpy"` is a validation error.

- [ ] **Step 1:** In the `mtuem` repo (`..`), create branch `py-redesign` and commit the current `mtuem-py` tree as the baseline (excluding `.venv`, `__pycache__`, `output/`), so the redesign diff is reviewable.
- [ ] **Step 2: Write failing tests** `test_control_defaults`, `test_settings_reject_unknown_field` (`EstimationSettings(maxIterations=5)` raises `ValidationError`), `test_settings_reject_nonpositive_iterations`, `test_numpy_backend_rejects_analytic_hessian`, `test_resolved_hessian_defaults`, `test_settings_json_round_trip` (`EstimationSettings.model_validate_json(s.model_dump_json()) == s`, with a nested `MultiStart`).
- [ ] **Step 3:** Run `uv run pytest tests/test_config.py -v`; expect import failure.
- [ ] **Step 4:** Implement `config.py`; update `pyproject.toml`; `uv sync`.
- [ ] **Step 5:** Run the new tests and the whole suite; expect PASS (old suite unchanged).
- [ ] **Step 6:** Commit `feat: pydantic configuration objects, jax required`.

### Task 2: Beta

**Files:**
- Create: `src/mtuem/beta.py`
- Test: `tests/test_beta.py`

**Interfaces:**
- Produces: `Param(value: float, fixed: bool = False, lower: float | None = None, upper: float | None = None)`; `Beta(params: dict[str, Param])` with `fix(names: str | list[str], value: float | None = None) -> Beta`, `free(names) -> Beta`, `set(**values: float) -> Beta`, `lower(names, bound: float) -> Beta`, `upper(names, bound: float) -> Beta`, `classmethod from_dict(values: Mapping[str, float], fixed: Iterable[str] = ()) -> Beta`, `names -> list[str]`, `free_names -> list[str]`, `fixed_names -> list[str]`, `values() -> dict[str, float]`, `bounds() -> dict[str, tuple[float | None, float | None]]`, `__getitem__(name) -> Param`, `__contains__`. Methods mutate and return `self` (chainable). A `locked: dict[str, str]` private attribute maps names that must stay fixed to a reason (used by Task 9).
- Name resolution `_resolve(name) -> list[str]`: exact match, else every key matching `^{re.escape(name)}_\d+$`; no match raises `KeyError` listing close names (`difflib.get_close_matches`).

- [ ] **Step 1: Write failing tests:** `test_fix_base_name_covers_all_classes` (`Theta` fixes `Theta_1`, `Theta_2`), `test_fix_exact_name`, `test_base_name_does_not_prefix_match` (`chk2` raises `KeyError`, does not touch `chk21_1`), `test_unknown_name_suggests` (message contains `Phi_1`), `test_lower_bound_replaces_constraint`, `test_start_outside_bounds_raises` (`ValueError` on `lower` above current value), `test_from_dict_marks_fixed`, `test_json_round_trip`, `test_locked_parameter_cannot_be_freed` (`ValueError` containing the lock reason).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement `beta.py`.
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit `feat: Beta parameter object`.

### Task 3: Move the engine under neutral names

**Files:**
- Create: `src/mtuem/_engine/__init__.py`, `_engine/inputs.py`, `_engine/latent.py`, `_engine/optimize.py`, `_engine/loglike.py`, `_engine/estimate.py`, `_engine/em.py`, `_engine/output.py`, `_engine/postestimation.py` (moved bodies of `apollo.py` sections)
- Delete: `src/mtuem/apollo.py`
- Modify: `src/mtuem/autodiff.py`, `multistart.py`, `postestimation.py`, `lc_profiled.py`, `likelihood.py` (imports, names)
- Modify: all tests importing `mtuem.apollo` (mechanical rename only)

**Interfaces:**
- Produces (private): `Inputs(database, beta: Mapping, fixed: list[str], control: Control, lc_pars: Callable | None, autodiff: bool)` (former `ApolloInputs`, always works in logs; `.attach`, `.matrix`, `.panel`, `.nIndivs`, `.nObs`, `.indiv_index`, `.first_rows`, `.ids`); `Estimate` dataclass (former `ApolloModel`, same fields minus `apollo_` prefixes: `beta_start`, `fixed`, `control`, `estimationRoutine` renamed `routine`); `estimate(beta, fixed, prob_fn, inputs, *, max_iterations, max_function_evals, hessian, saddle_escapes, silent, bounds=None) -> Estimate`; `lc_em(beta, fixed, prob_fn, inputs, *, max_iterations, tolerance, post_em: int, estimate_kwargs: dict, bounds=None) -> Estimate`; `class_alloc`, `lc`, `lc_conditionals`, `prediction`, `delta_method(est, expressions, varcov="robust", silent=False)`, `output_text(est, print_covar, print_corr) -> str`, `save_files(est, folder: Path, name: str) -> dict[str, Path]`.
- Converged runs report code 0 and message "successful convergence"; the BGW code table is removed. The report header reads "mtuem" with the package version, never "Apollo".

- [ ] **Step 1:** Rename imports and call sites in the existing tests to the private engine (`from mtuem import _engine as engine`); change `assert model.code == 4` to `== 0`. Run; expect FAIL (module missing).
- [ ] **Step 2:** Split `apollo.py` into the modules above, renaming only identifiers and strings; drop `estimationRoutine`, `bgw`, `workInLogs` (always logs), `noValidation`, `noDiagnostics`, `nCores`.
- [ ] **Step 3:** Run the full suite; expect every previously passing test to PASS with unchanged tolerances.
- [ ] **Step 4:** Commit `refactor: move estimation engine to mtuem._engine`.

### Task 4: Parameter bounds in BFGS and EM

**Files:**
- Modify: `src/mtuem/_engine/optimize.py`, `_engine/estimate.py`, `_engine/em.py`
- Modify: `tests/reference/make_reference.R` (add case `eq1_bounded`)
- Create: `tests/reference/eq1_bounded_ref.json`
- Test: `tests/test_bounds.py`

**Interfaces:**
- Consumes: `estimate(..., bounds: dict[str, tuple[float | None, float | None]] | None)`.
- Produces: `maximise_bounded(fun_vec, x0, lower: np.ndarray, upper: np.ndarray, **maximise_kwargs) -> OptResult`, the adaptive log-barrier of R's `constrOptim` (used by maxLik for Apollo inequality constraints):

```text
mu = 1e-4; outer_eps = 1e-5; outer_iterations = 100
g(x) = stacked (x_i - lower_i) and (upper_i - x_i) for finite bounds; require g(x0) > 0
repeat k = 1..outer_iterations:
    g_old = g(x_k)
    R(x) = f(x) + mu * sum(g_old * log(g(x)) - A x)      # A = d g / d x (rows of +-1)
    x_{k+1} = maximise(R) from x_k                         # -inf outside g > 0
    stop when |R(x_{k+1}) - R(x_k)| < outer_eps * (|R(x_k)| + outer_eps)
```

- [ ] **Step 1:** Add `eq1_bounded` to `make_reference.R` (eq1 with Apollo constraint `"thw > <value that binds>"`; pick the bound as the unconstrained `thw` estimate plus 0.1 of its absolute value) and run `Rscript tests/reference/make_reference.R eq1_bounded` from `mtuem-py`; check the JSON is written.
- [ ] **Step 2: Write failing tests:** `test_unbinding_bounds_match_unbounded` (estimates equal within `1e-8`), `test_binding_bound_matches_r` (estimate `rtol=1e-4`, LL `abs=1e-3` vs `eq1_bounded_ref.json`), `test_start_on_bound_raises`, `test_em_respects_bounds` (a 2-class `lc_cholesky` with `Phi_1` lower bound above its free optimum ends at or above the bound).
- [ ] **Step 3:** Run; expect FAIL.
- [ ] **Step 4:** Implement `maximise_bounded`; `estimate` and the EM M-step call it when any bound applies to a free parameter.
- [ ] **Step 5:** Run the bounds tests and the full suite; expect PASS.
- [ ] **Step 6:** Commit `feat: parameter bounds via adaptive barrier`.

### Task 5: Public likelihood functions and JAX ports of alpha-beta and additive quadratic

**Files:**
- Modify: `src/mtuem/likelihood.py` (rename public functions), `src/mtuem/autodiff.py` -> rename to `src/mtuem/likelihood_jax.py`
- Test: `tests/test_likelihood_jax.py` (extend from `test_autodiff.py`)

**Interfaces:**
- Produces: `mtuem.likelihood.theta_phi`, `.alpha_beta`, `.additive_quadratic`, `.residuals` (former `mtuem_likelihood`, `mtuem_ab_likelihood`, `additive_quadratic_mtuem_likelihood`, `mtuem_residuals`; same `(settings, functionality, inputs)` signature); `mtuem.likelihood_jax.theta_phi`, `.alpha_beta`, `.additive_quadratic`, `.class_alloc`, `.lc` returning per-row log-likelihoods. `mtuem.likelihood_jax.PORTED: frozenset[str]` lists the forms that have a JAX version; a form missing from it falls back to NumPy with `warnings.warn(..., FallbackWarning)` where `class FallbackWarning(UserWarning)` lives in `mtuem/_engine/__init__.py`.

- [ ] **Step 1: Write failing tests:** `test_alpha_beta_jax_matches_numpy` and `test_additive_quadratic_jax_matches_numpy` (log-likelihood equal `rel=1e-10` on `maed` at three parameter vectors, for profiled, `sig`, `sig`+`rho` and `cholesky` errors where the NumPy version supports them), `test_jax_gradient_matches_central_differences` (`rtol=1e-5`) for both forms.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Port both forms with `jnp`, mirroring the NumPy degeneracy guards (`-inf` where NumPy returns zero likelihood). If a form cannot be ported exactly, leave it out of `PORTED` and keep its test as `xfail` with the reason.
- [ ] **Step 4:** Run; expect PASS (or documented xfail).
- [ ] **Step 5:** Commit `feat: jax alpha-beta and additive quadratic likelihoods`.

### Task 6: JAX profiled latent class likelihood

**Files:**
- Modify: `src/mtuem/lc_profiled.py`, `src/mtuem/likelihood_jax.py`
- Test: `tests/test_lc_profiled_jax.py`

**Interfaces:**
- Produces: `mtuem.likelihood_jax.lc_profiled(settings, functionality, inputs)` with the NumPy signature; the inner fixed point (posterior weights -> weighted class covariances -> posteriors) runs with `jax.lax.while_loop` for the value, and the derivative uses implicit differentiation (`jax.custom_vjp` on the fixed point: solve `(I - dT/dz)^T v = cotangent`). Adds `"lc_profiled"` to `PORTED` only if the tests pass.

- [ ] **Step 1: Write failing tests:** `test_lc_profiled_jax_value_matches_numpy` (`rel=1e-10` at the `lc_profiled_ref.json` start, test and estimate points), `test_lc_profiled_jax_gradient_matches_central_differences` (`rtol=1e-5`).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run. If the gradient test cannot pass, mark both `xfail` with the reason, keep `lc_profiled` out of `PORTED`, and add `test_lc_profiled_falls_back_with_warning` (`pytest.warns(FallbackWarning)` when fitting on `backend="jax"`).
- [ ] **Step 5:** Commit `feat: jax profiled latent class likelihood` (or `...fallback`).

### Task 7: Declarative single-class `MTUEM`

**Files:**
- Create: `src/mtuem/models.py`
- Test: `tests/test_models_single.py`

**Interfaces:**
- Produces: `class MTUEM(pydantic.BaseModel)` with fields `work_times: list[str]`, `free_times: list[str] = []`, `free_goods: list[str] = []`, `goods_cost: dict[str, float | str] | None = None`, `Tc: str | float`, `Ec: str | float`, `w: str | float`, `tau: float = 168`, `utility: Literal["theta_phi", "alpha_beta", "additive_quadratic"] = "theta_phi"`, `covariance: Literal["profiled", "diagonal", "full", "cholesky"] = "profiled"`, `optimal_tw: bool = True`, `elasticity_covariates: dict[str, list[str]] = {}`, `extra: dict[str, Any] = {}` (passed through to the likelihood settings, for example `omitted_times`). Methods: `beta() -> Beta`; `columns() -> list[str]` (every data column used); `normalization_groups() -> dict[str, list[str]]` (for `MultiStart.normalization`); `_compile(backend: str) -> tuple[ProbFn, LcParsFn | None, bool]` (prob function, lc_pars, autodiff flag); `fit(data: pd.DataFrame, beta: Beta, control: Control, settings: EstimationSettings | None = None) -> Result` (Task 8 implements `Result`).
- Parameter names: work `Theta`, `Phi`, `thw` (theta-phi; alpha-beta and additive quadratic use the keys their likelihood settings read); times elasticities `th{k}` and goods `ph{k}` (1-based position in `free_times` / `free_goods`); `diagonal` adds `sig{i}`; `full` adds `sig{i}` and `rho{i}{j}` (i < j); `cholesky` adds `chk{i}{j}` (i >= j, row-major lower triangle); equations are numbered 1 = work, then free times, then goods. Covariate terms: `{elasticity}_{covariate}`, entering linearly. Defaults: `Theta=1` fixed, `Phi=1`, `thw=0`, elasticities `0.3`, `sig`/`chk` diagonal `10`, others `0`.
- `fit` validates data: missing columns, non-numeric or NaN values in `columns()` raise `ValueError` naming them (Review Focus 1).

- [ ] **Step 1: Write failing tests:** `test_beta_names_match_r_convention` (eq3 cholesky spec produces `Theta, Phi, thw, th1, ph1, chk11, chk21, chk22, chk31, chk32, chk33`), `test_compiled_loglik_matches_reference` (for every single-class case in `tests/reference/*_ref.json`, log-likelihood at the R start, test and estimate points equals the reference `rel=1e-10, abs=1e-8`, both backends; the reference names are mapped with a rename dict in `tests/cases.py`), `test_elasticity_covariates_enter_linearly`, `test_fit_rejects_nan_column`, `test_fit_rejects_missing_column`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement `models.py` (`fit` may raise `NotImplementedError` until Task 8; its data validation is tested through a private `_validate_data`).
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit `feat: declarative single-class MTUEM`.

### Task 8: `fit`, `Result`, saving and loading

**Files:**
- Create: `src/mtuem/results.py`, `src/mtuem/io.py`
- Modify: `src/mtuem/models.py`, `src/mtuem/postestimation.py` (functions take the engine `Estimate` plus inputs; renamed `values_of_time`, `cholesky_delta_method`, `get_corr`)
- Test: `tests/test_results.py`, rewrite `tests/test_r_parity.py` single-class parts on the public API

**Interfaces:**
- Produces: `class Result` with attributes `model`, `beta: Beta` (start), `control`, `settings`, `estimates: pd.Series`, `se`, `robse: pd.Series`, `varcov`, `robvarcov: pd.DataFrame`, `loglik: float`, `loglik_start: float`, `converged: bool`, `n_obs`, `n_indivs`, `n_free: int`, `aic`, `bic: float`, `starts: list[Result] | None`, `_estimate` (engine object); methods `summary(print_covar=False, print_corr=False) -> None`, `summary_text(...) -> str`, `save(folder: Path | None = None) -> dict[str, Path]` (engine files plus `<model_name>_model.json`), `to_frame() -> pd.DataFrame` (columns `Estimate, Std.err., t-ratio(0), Rob.std.err., Rob.t-ratio(0)`), `values_of_time(**kwargs) -> dict`, `predict() -> pd.DataFrame`, `delta_method(expressions: dict[str, Callable[[dict], float]], varcov="robust") -> pd.DataFrame`, `correlations() -> pd.DataFrame` (cholesky models). `mtuem.load(path: Path, data: pd.DataFrame | None = None, model_class: type | None = None) -> Result`; methods needing data raise `ValueError("... pass data= to mtuem.load")` when loaded without it.
- `_model.json` holds `{"mtuem_version", "model_class": "module:qualname", "model": model.model_dump(), "beta", "control", "settings", "estimates", "se", "robse", "varcov", "robvarcov", "loglik", "loglik_start", "converged", "n_obs", "n_indivs"}`.

- [ ] **Step 1: Write failing tests:** `test_fit_single_class_matches_r` (rewritten `test_estimation_matches_r` on `MTUEM.fit`, numpy backend tolerances; jax backend with SE `rtol=5e-3`), `test_summary_header_has_no_apollo`, `test_save_writes_engine_files_and_json`, `test_load_round_trip` (estimates, varcov equal; `values_of_time` works with `data=`), `test_load_without_data_methods_raise`, `test_load_unimportable_class_raises_with_hint` (Review Focus 5), `test_to_frame_columns`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement `MTUEM.fit` (single BFGS path; EM, fallback and multistart arrive in Task 10), `Result`, `io.load`.
- [ ] **Step 4:** Run the suite; expect PASS.
- [ ] **Step 5:** Commit `feat: Result object, save and load`.

### Task 9: `LatentClassMTUEM` and the `class_loglik` hook

**Files:**
- Modify: `src/mtuem/models.py`, `src/mtuem/results.py` (`posterior`, class shares)
- Test: `tests/test_models_latent.py`; rewrite the latent class parts of `tests/test_r_parity.py`

**Interfaces:**
- Produces: `class LatentClassMTUEM(MTUEM)` adding `n_classes: int`, `participation: list[list[bool]] | None = None` (one row per class over `free_times`; default all True), `allocation: list[str] = []`, `reference_class: int = 1`, `shared: list[str] = []` (base names without class suffix, for example the six `chk` names of `lc_shared`). Method `class_loglik(self, s: int, params: SimpleNamespace, data: Data) -> Array` returns the per-row log-likelihood of class `s` (1-based) built from the standard likelihood for the class's active equations; override it to customize. `Data` wraps the engine inputs: `data["Tf1"] -> Array`, `data[["Tf1", "Tf2"]] -> Array (N x k)`, attribute `xp` is `jax.numpy` or `numpy` for the active backend. `absent_times(s: int) -> list[str]`. `Result.posterior() -> pd.DataFrame` (columns `ID, class_1..class_S`).
- `beta()` builds per-class names with suffix `_s`, fixes `Theta_s`, fixes at 0 and locks (Beta `locked`) the elasticities and covariance terms of absent activities, fixes `asc` and `x_*` of `reference_class`, and puts `shared` names without suffix. `covariance="profiled"` uses the profiled latent class likelihood.

- [ ] **Step 1: Write failing tests:** `test_lc_cholesky_matches_r_at_r_optimum` and `test_lc_reaches_at_least_r_likelihood` for `lc_profiled`, `lc_shared` (ported from `test_r_parity.py`), `test_compiled_lc_loglik_matches_reference` (both backends, `rel=1e-10`), `test_absent_activity_parameters_are_locked` (Review Focus 3), `test_class_without_free_activities_fits` (Review Focus 2), `test_posterior_one_row_per_individual_on_panel_data` (Review Focus 4, duplicate each `maed` row twice under the same id), `test_class_loglik_override_is_used` (a subclass subtracting a constant shifts the LL by `N * const` per class mix as expected).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the suite; expect PASS.
- [ ] **Step 5:** Commit `feat: declarative latent class MTUEM with class_loglik hook`.

### Task 10: Estimation orchestration (EM, fallback, multistart)

**Files:**
- Modify: `src/mtuem/models.py` (shared `_fit` used by all models), `src/mtuem/multistart.py`
- Test: `tests/test_orchestration.py`; port `tests/test_ported_customMultiStart.py`

**Interfaces:**
- Consumes: `EstimationSettings`, `MultiStart`, engine `estimate`, `lc_em`.
- Produces: `method="em"` runs `lc_em(post_em=2)` (EM then BFGS); a raised `FloatingPointError`, `ValueError` or `LinAlgError`, or a non-converged result, triggers `fallback="bfgs"` from the same start, recorded in `Result.notes: list[str]`; `method="em"` on a single-class model raises `ValueError`. `multistart` runs `run_multistart(model, data, beta, control, settings) -> tuple[Result, list[Result | None]]` with exactly the semantics of `customMultiStart` (uniform draws in `[lower, upper]`, default start plus or minus 1; optional EM warmup; normalization groups from `model.normalization_groups()`; keep best converged non-saddle); `Result.starts` holds every attempt. No candidate qualifying raises `RuntimeError` with the attempts attached as `.attempts`.

- [ ] **Step 1: Write failing tests:** `test_em_then_bfgs_matches_direct_bfgs_optimum` (`lc_cholesky`, LL `abs=1e-3`), `test_fallback_runs_when_em_fails` (monkeypatched EM raising), `test_em_on_single_class_raises`, `test_multistart_same_seed_same_result`, `test_multistart_keeps_best_non_saddle`, `test_multistart_all_fail_raises_with_attempts`, plus the ported `customMultiStart` tests rewritten on the new API.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the suite; expect PASS.
- [ ] **Step 5:** Commit `feat: em, fallback and multistart from EstimationSettings`.

### Task 11: `CustomModel` (hand-written likelihood)

**Files:**
- Modify: `src/mtuem/models.py`
- Test: `tests/test_custom_model.py`

**Interfaces:**
- Produces: `class CustomModel` with `__init__(self, probabilities: Callable[[dict, Inputs, str], dict], lc_pars: Callable | None = None, backend: Literal["jax", "numpy"] = "numpy")` and the same `fit(data, beta, control, settings) -> Result`. `probabilities` returns a dict with `"model"` (and `"class_s"` components for latent class), per-row or per-individual log-likelihoods. The inputs object exposes `attach(beta) -> SimpleNamespace` (columns, parameters, lc_pars). `settings.backend` must match the declared backend, else `ValueError`.

- [ ] **Step 1: Write failing tests:** `test_custom_model_reproduces_declarative_eq3` (same estimates as `MTUEM` eq3 cholesky, `rtol=1e-8`), `test_custom_latent_class_with_lc_pars` (reproduces `lc_cholesky` LL at the R optimum), `test_custom_backend_mismatch_raises`.
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit `feat: CustomModel for hand-written likelihoods`.

### Task 12: Public surface, examples and README

**Files:**
- Modify: `src/mtuem/__init__.py`, `README.md`, `examples/maed_*.py` (rewrite on the new API), every remaining test still using old names
- Test: `tests/test_public_api.py`

**Interfaces:**
- `mtuem.__all__` = `Control, EstimationSettings, MultiStart, Beta, Param, MTUEM, LatentClassMTUEM, CustomModel, Result, load, likelihood, likelihood_jax, FallbackWarning, DATASETS, load_dataset, maed`, plus the `get_*` utility functions renamed only where they carry `mtuem_`.

- [ ] **Step 1: Write failing tests:** `test_no_public_name_mentions_apollo` (walk `mtuem` and its public submodules; no name, module or file under `src/mtuem` contains "apollo", case-insensitive), `test_examples_run` (each `examples/maed_*.py` runs in a subprocess with `uv run python` and exits 0; mark `slow` for the latent class ones).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Update exports, rewrite examples and README (usage with the new API; a section crediting Apollo as the method the engine reproduces; the Biogeme note kept).
- [ ] **Step 4:** Run the full suite plus `uv run pytest --cov=mtuem --cov-report=term-missing`; expect PASS and coverage at least 80%.
- [ ] **Step 5:** Commit `feat: new public API, examples and README`.

### Task 13: `ActivityParticipationMTUEM` example

**Files:**
- Create: `examples/activity_participation.py`
- Test: `tests/test_activity_participation.py`

**Interfaces:**
- Produces: `class ActivityParticipationMTUEM(LatentClassMTUEM)` with fields `penalization: float`, `likelihood_floor: float = 0.0`; `class_loglik` returns `super().class_loglik(...) - penalization * data[absent].sum(axis=1)`, then `xp.maximum(ll, log(likelihood_floor))` when the floor is positive. Also `def article_model(work_times, free_times, participation, allocation, penalization, likelihood_floor=0.0) -> ActivityParticipationMTUEM` and `def article_beta(model) -> Beta` reproducing `get_apollo_beta()` / `get_apollo_fixed()` of `R/models/hobbies_cholesky.R` (2 classes) and `R/models/hobbies_rest.R` (any number of classes). A `__main__` block fits the 2015 model at penalization 5 on the fixture data.

- [ ] **Step 1: Write failing tests:** `test_penalty_matches_r_formula` (class LL equals the unpenalized class LL minus `penalization * hours of absent activities`, rowwise, `rel=1e-12`), `test_floor_applies_in_logs`, `test_article_beta_matches_r_start_values` (names, values and fixed set equal the fixture JSON of Task 14 for 2015 and 2022-rest).
- [ ] **Step 2:** Run; expect FAIL.
- [ ] **Step 3:** Implement with commented code meant for learning.
- [ ] **Step 4:** Run; expect PASS.
- [ ] **Step 5:** Commit `docs: activity participation example`.

Note: Task 13's last test reads the Task 14 fixtures; implement Task 14 Steps 1-2 first if running strictly in order.

### Task 14: Article reference fixtures (R)

**Files:**
- Create: `tests/reference/make_article_reference.R`
- Create: `tests/reference/article/<run>_p<k>.parquet` and `<run>_p<k>.json` for `2015_p5`, `2015-new_p5`, `2022_p5`, `2022-bigcities_p4`, `2022-rest_p5`

**Interfaces:**
- Produces per target JSON: `run, penalization, likelihood_floor, warm_started, indiv_id, times_internal, act_mix, socioecon, beta_start, fixed, r_estimate, r_se, r_robse, r_maximum, r_minutes`. The parquet is `prepared_data$database` restricted to the columns the model uses plus the id.
- The script runs from the article root (`Rscript <mtuem-py>/tests/reference/make_article_reference.R <run> <penalization>`), sets `study_year`, sources the same files as `scripts/<run>/experiments/fit_penalization_grid.R` (rest model and bigcities location handled as there), calls `prepare_hobbies_data`, the matching `initialise_*`, `get_apollo_beta()`, `get_apollo_fixed()`, and reads R results from `output/<run>/models/penalization[/conce_valpo]/corner_hobbies_p<k>_{estimates.csv,model.rds,convergence.csv}`. `likelihood_floor` and `warm_started` come from the convergence CSV when present.

- [ ] **Step 1:** Write the script; run it for the five targets.
- [ ] **Step 2:** Verify each JSON parses, `beta_start` keys equal `r_estimate` keys, and the parquet row count equals R's `nObs` from `_output.txt`.
- [ ] **Step 3:** Commit `test: article reference fixtures`.

### Task 15: Article replication (slow)

**Files:**
- Create: `tests/test_article_replication.py`

**Interfaces:**
- Consumes: Task 13 `article_model`, `article_beta`; Task 14 fixtures.

- [ ] **Step 1: Write the tests:** `test_article_replication[run-backend]` for the five targets and `backend in ("jax", "numpy")`, `@pytest.mark.slow`: build the model, start from `beta_start`, `EstimationSettings(method="em", fallback="bfgs", max_iterations=500, backend=backend)`; assert estimates `rtol=1e-4, atol=1e-4` and SEs and robust SEs at the backend's tolerance against the R values. On failure the message lists the LL difference and the ten largest parameter differences. Write a CSV summary (target, backend, LL Python, LL R, max abs estimate diff, max rel SE diff, minutes) to `output/replication/`.
- [ ] **Step 2:** Run `uv run pytest tests/test_article_replication.py --runslow -v` (hours). Expect PASS; any target landing on a different optimum is reported as a failure, never loosened.
- [ ] **Step 3:** Commit `test: activity participation article replication`.
