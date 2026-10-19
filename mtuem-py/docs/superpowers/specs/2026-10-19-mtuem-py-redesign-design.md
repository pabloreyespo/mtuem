# mtuem (Python) redesign: a tailored, Python-native API

Date: 2026-10-19
Status: draft for review

## Intent

`mtuem-py` is a working Python port of the R package `mtuem`. Because Apollo has no Python version,
the port carries its own estimation engine in `mtuem.apollo` and keeps Apollo's names and calling
conventions (`apollo_beta`, `apollo_fixed`, `apollo_control`, `apollo_probabilities`,
`validate_inputs`, `customMultiStart`, ...). Model scripts look like Apollo scripts.

The package only serves microeconomic time-use (and expenditure) models, so most of what Apollo
asks the user to write by hand can be derived from a short model description. Pablo wants:

1. No public name that says "apollo". Keep only the parts of the engine these models need, under
   mtuem names.
2. The dictionaries a user fills in today (control, beta and fixed, estimation settings) become
   typed pydantic objects that are imported, filled and passed to estimation.
3. A more Python-natural way of working (sklearn/statsmodels style), with every design decision
   approved by Pablo (the decisions below are his unless marked "builder's call").
4. No loss of functionality or numerical precision, and the results of the article
   `articles/mtuem-lc-activity-participation` must still be reproducible.

## Goals

- A declarative model API: the user describes the model (times, covariance, classes,
  participation patterns, allocation covariates) and the package builds the likelihood,
  parameter names, Cholesky terms and class allocation.
- Typed configuration with pydantic: `Control`, `Beta`, `EstimationSettings`, `MultiStart`.
- A result object (statsmodels style) with output, saving, loading and post-estimation methods.
- JAX as the default engine (exact derivatives), NumPy kept as an alternative engine.
- An override hook and a fully hand-written likelihood path, so anything possible today stays
  possible.
- Replication of the article's selected models against the saved R results.

## Non-goals

- Backward compatibility with the current `mtuem.apollo` API. It is replaced, not kept beside the
  new one (Pablo's decision; nothing outside `mtuem-py` depends on it, it is untracked and at 0.1.0).
- Changing the R package `mtuem`.
- Reimplementing the article's data preparation (`R/data/enut_data.R`, `hobbies_config()`) in
  Python. The replication uses the databases exactly as R prepared them (see Replication).
- Replicating the committed-time article (`articles/2026-mtuem-lc-committed-time`). Its code never
  ran on this library; the design must make that model expressible, but there is no R result to
  match.

## Anti-goals (would be a failure even if it runs)

- Any drift in log-likelihood or estimates against the R reference results beyond the tolerances
  of the current tests.
- A public function, class, module, file name or printed header that says "apollo". The docs may
  still cite Apollo (Hess and Palma) as the method the engine reproduces.
- A declarative API that cannot express a model the current port can estimate.

## Decisions taken with Pablo

### Naming

- Plain names under the package: `mtuem.estimate`-style, no `mtuem_` prefix on the new API
  (`mtuem.MTUEM`, `mtuem.Control`, `mtuem.load`, ...).
- The parameter object keeps the name **beta**: class `Beta`, created with `model.beta()`.
- Parameter names keep the R convention `name_class` (`Theta_1`, `chk21_2`, `x_female_2`) so that
  outputs line up with the R runs.

### Configuration objects (pydantic)

- `Control` replaces `apollo_control`: model name, description, individual id column, output
  directory, seed, and the remaining control flags the engine uses today.
- `Beta` replaces `apollo_beta` + `apollo_fixed` + constraint strings. One entry per parameter with
  start value, fixed flag and lower/upper bound. Helpers:

  ```python
  beta = model.beta()                       # every name the spec needs, with sensible starts
  beta.fix("Theta", 1)                      # base name: all classes (Theta_1, Theta_2)
  beta.fix(["th2_1", "chk31_1"], 0)         # single entries
  beta.set(Phi_1=0.8, Phi_2=0.9)            # start values
  beta.lower("Phi", 1e-10)                  # replaces "Phi_1 > 0 + 1e-10"
  beta = mtuem.Beta.from_dict(old_apollo_beta, fixed=old_apollo_fixed)
  ```

  Names that do not exist in the spec are rejected with a clear error.
- `EstimationSettings` replaces `estimate_settings` and the hand-chained EM/fallback/multistart
  logic:

  ```python
  settings = mtuem.EstimationSettings(
      method="em",            # "em" or "bfgs"; EM refines with BFGS at the end
      fallback="bfgs",        # used when EM fails (replaces the R tryCatch)
      hessian="analytic",     # "analytic" (exact, JAX) or "numeric" (numDeriv port)
      max_iterations=500,
      backend="jax",          # default; "numpy" is the alternative
      multistart=mtuem.MultiStart(n=20, seed=123),   # None = off
  )
  ```

  `MultiStart` mirrors today's `customMultiStart` exactly (builder's call): number of candidates
  drawn within the start values plus or minus one, seed, optional EM warmup (`first_em`,
  `em_iter_max`), theta-phi or alpha-beta normalization, and it keeps the single best non-saddle
  optimum. Every attempt is available on the result for inspection.
- Bounds (`beta.lower`, `beta.upper`) are a new optimizer feature, since the current port has no
  constraint support. Builder's call: reproduce Apollo's handling of inequality constraints so
  results match R runs that used them, and apply the same treatment in the EM M-steps. None of
  the replication targets uses constraints.

### Model API (sklearn-like description, statsmodels-like fit)

```python
model = mtuem.LatentClassMTUEM(
    work_times=["Tw"], free_times=["Tf1", "Tf2", "Tf3"],
    Tc="Tc", Ec="Ec", w="w", tau=168,
    covariance="cholesky", n_classes=2,
    participation=[[True, False, True], [True, True, True]],
    allocation=["female", "mayor45", "university"],
)
beta = model.beta()
result = model.fit(df, beta, control, settings)   # returns a new result; model stays reusable
result.summary(); result.save(folder); result.values_of_time(); result.posterior()
```

- `fit` returns a separate result object (statsmodels style). The model is a reusable spec, so a
  penalization grid is a loop that collects results.
- The package builds per-class parameters, Cholesky terms for the active equations of each class,
  fixes the parameters of absent activities, and builds the class allocation utilities (class 1
  as reference by default).

### Custom models

- `LatentClassMTUEM.class_loglik(self, s, params, data)` returns the log-likelihood of class `s`.
  Overriding it (like `forward` in a PyTorch module) customizes one class's likelihood while
  estimation, EM, output and post-estimation keep working.
- A fully hand-written likelihood function remains available for models the declarative API does
  not cover (the successor of `apollo_probabilities`).
- `ActivityParticipationMTUEM` (participation penalty `ll - penalization * hours of absent
  activities`, plus the optional likelihood floor for exact reproduction of old R runs) lives in
  `examples/` as a commented, runnable example. It is not part of the core package.

### Estimation engine

- JAX is a required dependency and the default backend. NumPy stays available via
  `backend="numpy"` and keeps the numDeriv-style numerical derivatives of today.
- One optimizer only: the engine's BFGS with BHHH seeding and saddle escape. The `bgw` label and
  the `estimationRoutine` switch are removed (in the current port both labels ran the same
  algorithm).
- The alpha-beta and additive quadratic likelihoods have no JAX version today; they follow the
  same rule as the profiled model below (builder's call, applying Pablo's rule): port to JAX,
  else fall back to NumPy with a warning.
- The profiled latent class likelihood is ported to JAX where possible (implicit differentiation
  of its fixed point). Whatever cannot be ported falls back to NumPy automatically with a warning.
- The numeric core (likelihood formulas, optimizer, Hessian routine, EM, output formatting) is
  moved, not rewritten; the R reference tests are rewritten against the new API with the same
  tolerances.

### Outputs

Agreed by default (Pablo did not object):

- `result.summary()` prints today's Apollo-layout report, with an mtuem header.
- `result.save(folder)` writes the same files as today (`_output.txt`, `_estimates.csv`,
  `_covariance.csv`, ...) plus `_model.json` with spec, `Beta`, `Control`, settings and estimates.
- `mtuem.load(path)` restores a full result from it (replaces `apollo_loadModel`).
- `result.to_frame()` returns the estimates table as a pandas DataFrame.

### Replication targets

Refit in Python from the same default starting values as the R scripts (`get_apollo_beta()`),
with the full EM then maximum likelihood sequence, and compare estimates and standard errors
(classical and robust) against the saved R results in
`articles/mtuem-lc-activity-participation/output/<run>/models/penalization/`. Post-estimation
quantities (implied sigma and rho, values of time, posterior class shares) are out of scope for
the replication tests. All targets use the engine's single BFGS optimizer (Pablo: "stick with
just BFGS"). The 2015 R result came from Apollo's `bgw`; at the same optimum the engine's
standard errors are the same estimator (inverse Hessian, BHHH meat for robust), so 2015 is
compared like the others.

| Run | Penalization | Model |
|---|---|---|
| 2015 | 5 | `hobbies_cholesky`, 2 classes |
| 2015-new | 5 | `hobbies_cholesky`, 2 classes |
| 2022 | 5 | `hobbies_cholesky`, 2 classes |
| 2022-bigcities | 4 | `hobbies_cholesky` with `conce`, `valpo` allocation (`penalization/conce_valpo/`) |
| 2022-rest | 5 (builder's call, see below) | `hobbies_rest`, four classes, rest Tf5 as free time |

These are slow tests (marked, run on demand). The R fits took 10 to 15 minutes (2015, 2015-new,
2022-bigcities), 37 minutes (2022) and 244 minutes (2022-rest, 4 classes, 82 free parameters), so
the suite may take hours; that is accepted. If a refit lands on a different optimum than R, the
test fails and the difference is reported, never loosened away. Standard errors under the JAX
backend are compared at a slightly looser tolerance than under NumPy, since R uses numerical
derivatives.

## Constraints

- Python via uv only, 3.14 pinned (`uv python pin`), run with `uv run`.
- Dependencies: numpy, scipy, pandas, pyarrow, pydantic, jax (with x64 enabled).
- Tests written first; 80% coverage minimum.
- No em dashes in docs and comments; no emojis.

## What the current code has that must survive

Every feature below exists in the port today and must remain reachable, either declaratively or
through the override/hand-written path:

- Likelihoods: Theta-Phi MTUEM (`mtuem_likelihood`), alpha-beta (`mtuem_ab_likelihood`),
  additive quadratic, profiled latent class (`mtuem_lc_profiled_likelihood`), residuals.
- Error structures: profiled, `sig`, `sig` + `rho`, `cholesky`; covariance shared across classes
  (`lc_shared` reference case).
- Goods equations (`free_goods`, `goods_cost`), `optimal_tw`.
- Latent class: class allocation, EM, posteriors/conditionals, prediction.
- Post-estimation: values of time (with posterior weights), correlations, Cholesky delta method,
  generic delta method, profiled sigma and information criteria.
- Multistart (`customMultiStart`), the bundled datasets (`maed`, `enut_i`, `enut_ii`,
  `enut_bigcities` and raw variants).

## Left to the builder (builder's calls, Pablo can override)

- **2022-rest penalization:** p4 and p5 R results both exist and no analysis script names one;
  the builder targets p5, consistent with the other 2022 runs, and also checks p4 if cheap.
- **Data for replication:** an R script (like `tests/reference/make_reference.R`) exports each
  target's prepared database and the R estimates to parquet/JSON fixtures, so Python tests do not
  depend on R at run time and do not reimplement the article's data preparation.
- **Covariates in elasticities:** the declarative spec accepts linear covariate terms on
  elasticities (for example `covariates={"Phi": ["female"], "thw": ["female"]}`, creating
  `Phi_female_1`, ...), as in `apollo_timeuse/examples/enut_1eq_sleep_covariates.r`. Derived
  columns (dummies) are prepared by the user in the DataFrame beforehand.
- **Model family:** a single-class `MTUEM` and `LatentClassMTUEM`, with `utility="thph" | "ab" |
  "additive_quadratic"` and `covariance="profiled" | "diagonal" | "full" | "cholesky"`, plus a
  `shared=[...]` option for parameters common to all classes. Exact names and argument layout are
  the builder's.
- Module layout, private function names, the exact pydantic field names, and how the hand-written
  likelihood is passed to `fit`.
- Porting the existing examples (`examples/maed_*.py`) to the new API and updating the README.
