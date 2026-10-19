# mtuem (Python)

A Python port of the R package [mtuem](https://github.com/pabloreyespo/mtuem), which estimates
Microeconomic Time-Use-Expenditure Models (Jara-Diaz et al., 2008). It covers the same likelihoods,
latent class variants, multi-start module, values of time and datasets as the R version, and prints
and saves results in the same layout as Apollo.

The R package only supplies likelihood functions and relies on Apollo for estimation. Apollo has no
Python version, so `mtuem.apollo` is a small engine with the parts of Apollo these models need:
`validate_inputs`, `estimate`, `lc_em`, `model_output`, `save_output`, `delta_method`, `prediction`,
`class_alloc`, `lc` and `lc_conditionals`. Model scripts keep Apollo's structure, so porting an R
script is mostly mechanical (compare `examples/*.py` with `../examples/*.r`).

## Why not Biogeme

Biogeme builds the likelihood one observation at a time from its expression language. Several MTUEM
specifications need statistics of the whole sample at every evaluation: the profiled covariance is
the covariance of the current residuals, and the profiled latent class model runs an inner fixed
point over all observations. Neither fits Biogeme's row-wise model, so the engine is built on
numpy and scipy instead, which also allows reproducing Apollo's output layout and standard errors
exactly.

## Installation

```bash
cd mtuem-py
uv sync
uv run python examples/maed_1eq.py
```

Requires Python 3.14 (pinned with `uv python pin`), numpy, scipy, pandas and pyarrow. JAX is
optional, for exact derivatives: `uv sync --extra jax`.

## Usage

```python
import mtuem
from mtuem import apollo

database = mtuem.maed()
database = database[database["EcI"] > 0]

apollo_beta = {"Theta": 1, "Phi": 1, "thw": 0}
apollo_fixed = ["Theta"]
apollo_inputs = apollo.validate_inputs(
    database, apollo_beta, apollo_fixed,
    {"modelName": "maed-1eq", "indivID": "PeID", "outputDirectory": "output"},
)

def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
    b = apollo_inputs.attach(apollo_beta)          # database columns + parameters
    settings = {
        "work_times": ["Tw"],
        "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
        "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
    }
    P = {"model": mtuem.mtuem_likelihood(settings, functionality, apollo_inputs)}
    return apollo.prepare_prob(P, apollo_inputs, functionality)

model = apollo.estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
apollo.model_output(model)
apollo.save_output(model)                          # <modelName>_output.txt, _estimates.csv, ...
vot = mtuem.mtuem_values_of_time(model, apollo_inputs, Ec="EcI")
```

Differences from the R calling conventions:

| R | Python |
|---|---|
| `apollo_attach(apollo_beta, apollo_inputs)` | `b = apollo_inputs.attach(apollo_beta)`, then `b.Theta`, `b.Tc`, ... |
| likelihoods read `apollo_inputs` from the calling frame | pass it as the third argument |
| settings are R lists | settings are dicts (lists for positional values such as `cholesky`) |
| `apollo_lcPars <- function(apollo_beta, apollo_inputs)` | `def apollo_lcPars(b, apollo_inputs)`, `b` holds parameters and columns |
| `apollo_deltaMethod` with string expressions | `apollo.delta_method(model, {"expression": {"VoL": lambda p: ...}})` |
| `customMultiStart(...)` uses the global R seed | pass `seed=` |
| `best_model` is `NA` when no candidate qualifies | `best_model` is `None` |

## Working in logs (default)

R's mtuem returns densities, `exp(ll)`, and Apollo takes their log again. A density below about
`exp(-745)` is 0 in double precision, so a single badly fitted observation makes the whole
log-likelihood `-Inf`. In mtuem-py, `apollo_control$workInLogs` (Apollo's own option) defaults to
`TRUE`. Every likelihood function (`mtuem_likelihood`, `mtuem_ab_likelihood`,
`additive_quadratic_mtuem_likelihood`, `mtuem_lc_profiled_likelihood`) then returns
log-likelihoods. `apollo.lc` mixes classes with log-sum-exp, and `apollo.class_alloc` carries
exact log-probabilities. Nothing is exponentiated, so nothing underflows.

On ENUT II, two workers with absurd wages make the density-mode likelihood `-Inf` at the starting
values, the same failure as in R. In logs, the 5,452-worker model estimates normally. Where both
modes are finite they agree to about 1e-14.

Consequence for model code: if you modify a component yourself, do it in logs. For example a
participation penalty becomes `ll - penalization * hours` instead of
`L * exp(-penalization * hours)`. `apollo_inputs.work_in_logs` tells you which mode you are in.
Set `workInLogs = False` in `apollo_control` to get densities, as in R. Inputs that are not
built by `validate_inputs` (a bare namespace with a `database`) also get densities.

## Exact derivatives with JAX (optional)

By default, gradients, individual scores and Hessians are numerical, as in Apollo. With the
`jax` extra (`uv sync --extra jax`), a model written with `mtuem.autodiff` gets exact ones.
Use `ad.mtuem_loglik`, `ad.class_alloc` and `ad.lc` in place of `mtuem.mtuem_likelihood`,
`apollo.class_alloc` and `apollo.lc`, and pass `autodiff=True` to `validate_inputs`:

```python
from mtuem import autodiff as ad

def apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate"):
    b = apollo_inputs.attach(apollo_beta)
    settings = {...}                                   # the same settings as before
    P = {"model": ad.mtuem_loglik(settings, functionality, apollo_inputs)}   # log-likelihood
    return apollo.prepare_prob(P, apollo_inputs, functionality)

apollo_inputs = apollo.validate_inputs(database, apollo_beta, apollo_fixed, control, autodiff=True)
model = apollo.estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
```

Everything else (`lc_em`, `customMultiStart`, output, prediction, values of time) works
unchanged. The functions always return log-likelihoods (the autodiff backend requires `workInLogs = TRUE`). Model code must use `jax.numpy` and must not branch on parameter values
with Python `if`.

Covered: the Theta-Phi likelihood with profiled, `sig`, `sig` + `rho` or `cholesky` errors,
class allocation and latent class mixtures (including panel data). Not covered:
`mtuem_lc_profiled_likelihood` (its inner fixed point is not differentiable as written), the
alpha-beta and the additive quadratic likelihoods.

On the 2022-bigcities participation model at penalization 4 (2 classes, 33 parameters, 2,361
workers), all three reach LL -28,778.5583 with standard errors within 0.15% of each other:

| | Time |
|---|---|
| R (Apollo, EM + BGW) | 637 s |
| mtuem-py, numerical derivatives (EM + BFGS) | 234 s |
| mtuem-py, autodiff (BFGS, no EM needed) | 3.9 s |

On small models (a few parameters, under a thousand rows) the numerical path is just as fast,
because JIT compilation takes one to three seconds.

## Functions

Likelihoods: `mtuem_likelihood` (Theta-Phi), `mtuem_ab_likelihood` (alpha-beta),
`additive_quadratic_mtuem_likelihood`, `mtuem_lc_profiled_likelihood`, `mtuem_residuals`.

Post-estimation: `mtuem_values_of_time`, `mtuem_get_corr`, `mtuem_cholesky_deltaMethod`,
`cholesky_to_cov`, `mtuem_get_profiled_sigma`, `mtuem_profile_ic`.

Estimation: `customMultiStart`, plus the `mtuem.apollo` engine.

Allocation helpers: `get_tw_*`, `get_ti_*`, `get_xi_*`, `get_values_of_time_*` for `thph`, `albe`
and `additive_quadratic`, `get_additive_quadratic_essentials`, `get_cond_err`.

Data: `mtuem.maed()` and `mtuem.load_dataset(name)` for `maed`, `enut.i`, `enut.i.raw`, `enut.ii`,
`enut.ii.raw`, `enut.bigcities` and `enut.bigcities.raw`, converted from the R package's `.rda`
files by `data-raw/convert_data.py`.

## Agreement with the R package

`tests/test_r_parity.py` compares against results produced by the R package
(`tests/reference/make_reference.R`, R 4.4.3, Apollo 0.3.9, mtuem 1.2.2) for the eight README
models:

- Log-likelihoods agree to about 1e-11 at the starting values, at a perturbed vector and at the
  R estimates, for all eight models.
- For the five single-class models, estimation from the README starting values reaches the same
  optimum: estimates agree to about 1e-5, classical and robust standard errors to about 0.01%.
- Latent class standard errors at the R optimum agree within 0.03%.
- Values of time, delta method results, residual correlations, predictions and the profiled
  covariances match.
- The estimates table prints character for character as `apollo_modelOutput` prints it.

From the README starting values, a single latent class estimation is path dependent in both
languages. The profiled and shared covariance models reach a higher likelihood in Python than
in R (R stops with false and singular convergence). The per-class Cholesky model can stop at a
different local optimum. This is why the examples use `customMultiStart`.

## Where it deliberately differs from Apollo

- `workInLogs` defaults to `TRUE` (Apollo: `FALSE`), so likelihood functions return
  log-likelihoods. See "Working in logs".
- `estimationRoutine = "bgw"` runs a BFGS line-search optimiser whose curvature is seeded with
  the BHHH matrix, as BGW does. It reports BGW's codes (4 is relative function convergence), so
  `customMultiStart`'s acceptance rule works unchanged. `"bfgs"` runs the same optimiser and
  reports code 0 on success, as maxLik does.
- When the optimiser stops at a saddle point (for example the symmetric solution of a latent class
  model started from identical classes), `estimate` steps along the direction of positive
  curvature and re-estimates, up to `saddleEscapes` times (default 3). Apollo escapes such
  points only by numerical noise.
- When `numDeriv::hessian`'s first step (10% of each parameter) leaves the region where the
  likelihood is positive (for example a correlation of -0.96 stepped past -1), the Hessian is
  recomputed with relative steps of 1e-4. Apollo falls back to maxLik's Hessian in that case.
- Parameters are scaled by `max(|start|, 1)`, not `|start|`, so near-zero starting values (common
  after EM) do not freeze those parameters.
- `LL(start)` is the log-likelihood at the starting values. Apollo prints 0 for these models,
  because its validation step uses the `validate` functionality, which returns 1 for every row.
- `BIC` uses the number of rows. Apollo prints `-Inf` because it counts zero modelled outcomes
  for a continuous density.
- `lc_em` maximises the expected complete-data log-likelihood over all parameters at once, so,
  unlike `apollo_lcEM`, it accepts parameters shared across classes.
- `apollo.prediction` works for latent class models built with `apollo.lc`, returning per-class
  allocations and class allocation probabilities. `apollo_prediction` fails on them in R.
- Like the R version, `get_values_of_time_thph` and `get_values_of_time_albe` use a fixed
  168-hour week in the leisure denominator, whatever `tau` is.

## Tests

```bash
uv run pytest
```
