# mtuem (Python)

Estimate Microeconomic Time-Use-Expenditure Models (Jara-Diaz et al., 2008) in Python. This is
the Python version of the R package [mtuem](https://github.com/pabloreyespo/mtuem): the same
likelihoods (Theta-Phi, alpha-beta, additive quadratic), error structures, latent class models,
multistart, values of time and datasets, with its own estimation engine.

You describe the model (its time and expenditure columns, utility form, error covariance,
classes and class allocation), take its parameters with `model.beta()`, and call `model.fit`,
which returns a `Result` with the report, the saved files and the post-estimation quantities.

## Installation

mtuem needs Python 3.14 and is managed with [uv](https://docs.astral.sh/uv/):

```bash
cd mtuem-py
uv python pin 3.14
uv sync
uv run python examples/maed_1eq.py
```

The dependencies (numpy, scipy, pandas, pyarrow, pydantic and JAX) are installed by `uv sync`.
JAX is required: it is the default backend.

## A single-class model

```python
import mtuem

data = mtuem.maed()
data = data[data["EcI"] > 0]

model = mtuem.MTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"],
    Tc="Tc", Ec="EcI", w="w", tau=168,
    utility="theta_phi",          # or "alpha_beta", "additive_quadratic"
    covariance="cholesky",        # or "profiled", "diagonal", "full"
)

beta = model.beta()               # every parameter, with default start values
print(beta.names)                 # Theta, Phi, thw, th1, ph1, chk11, chk21, chk22, chk31, ...
beta.set(th1=0.5, ph1=0.5, chk11=100, chk22=100, chk33=100)   # start values
beta.fix("Theta", 1)              # model.beta() already fixes Theta at 1 for identification

control = mtuem.Control(model_name="maed-3eq", indiv_id="PeID", output_dir="output")
result = model.fit(data, beta, control)

result.summary()                  # the estimation report
result.to_frame()                 # estimates, standard errors and t-ratios as a DataFrame
result.values_of_time()           # VoL and VTAW with delta method standard errors
result.correlations()             # sds and correlations implied by the Cholesky factor
result.predict()                  # predicted allocations
result.save()                     # report, csv files and maed-3eq_model.json in output/
```

Parameters follow the names of the R scripts: equations are numbered 1 = work, then the free
times, then the goods; `th1` is the elasticity of the first free time, `ph1` of the first good,
`sig2` the standard deviation of equation 2, `rho12` a correlation and `chk21` an entry of the
Cholesky factor (row by row over the lower triangle). `elasticity_covariates={"Phi": ["female"]}`
adds the parameter `Phi_female`, so that the elasticity is `Phi + female * Phi_female`.

`model.fit` does not change the model or `beta`: fit the same model again with other settings
or data, for example over a grid of specifications, and collect the results.

## Latent class models

```python
model = mtuem.LatentClassMTUEM(
    work_times=["Tw"], free_times=["Tf1"], free_goods=["Ef1"],
    Tc="Tc", Ec="EcI", w="w",
    covariance="cholesky",
    n_classes=2,
    allocation=["female", "older45"],     # class membership logit, class 1 is the reference
)
beta = model.beta()     # Theta_1, Phi_1, ..., chk33_1, asc_1, x_female_1, ..., Theta_2, ...
beta.set(chk22=20)      # a base name sets every class: chk22_1 and chk22_2
result = model.fit(data, beta, control, mtuem.EstimationSettings(method="em"))

result.class_shares()
result.posterior()                        # posterior class membership, one row per individual
result.values_of_time()                   # by class, weighted by the posteriors
```

Parameters of class `s` carry the suffix `_s`. `shared=["chk11", "chk21", ...]` makes
parameters common to every class (no suffix). `participation` says which free times each class
takes part in: for a class that does not take part in a free time, its elasticity and covariance
terms are fixed at 0 and locked. With `covariance="profiled"` the class covariances are
concentrated out of the likelihood (`result.profiled_covariance()` reports them, with
information criteria that count them).

To change the likelihood of a class, subclass and override `class_loglik(s, params, data)`;
estimation, EM, output, posteriors and values of time keep working (predictions use the standard
likelihood). On the JAX backend write it with `data.xp` (`jax.numpy`).

## A hand-written likelihood

For models the declarative classes do not cover, `CustomModel` takes a probabilities function
that returns the log-likelihood of each row under `"model"`:

```python
from mtuem import likelihood_jax

def probabilities(params, inputs, functionality="estimate"):
    b = inputs.attach(params)          # data columns and parameters as attributes
    settings = {
        "work_times": ["Tw"],
        "work_elasticities": {"Theta": b.Theta, "Phi": b.Phi, "thw": b.thw},
        "Tc": b.Tc, "Ec": b.EcI, "w": b.w, "tau": 168,
    }
    return {"model": likelihood_jax.theta_phi(settings, functionality, inputs)}

model = mtuem.CustomModel(probabilities, backend="jax", columns=["Tw", "Tc", "EcI", "w"])
beta = mtuem.Beta.from_dict({"Theta": 1, "Phi": 1, "thw": 0}, fixed=["Theta"])
result = model.fit(data, beta, control)
```

A latent class custom model also returns its classes (`"class_1"`, ...) and passes
`lc_pars=` (needed for EM). The likelihood functions are `mtuem.likelihood_jax` (JAX) and
`mtuem.likelihood` (NumPy): `theta_phi`, `alpha_beta`, `additive_quadratic`, the profiled latent
class `lc_profiled`, and `class_alloc` and `lc` for class allocation and mixtures.

## Estimation settings

```python
settings = mtuem.EstimationSettings(
    method="em",              # "bfgs" (default) or "em": EM, then BFGS from the EM estimates
    fallback="bfgs",          # BFGS from the start values when EM fails (recorded in result.notes)
    backend="jax",            # "jax" (default) or "numpy"
    hessian=None,             # "analytic" on JAX, "numeric" on NumPy by default
    max_iterations=500,
    multistart=mtuem.MultiStart(n=20, seed=123, normalization="theta_phi"),
)
```

**JAX or NumPy.** With `backend="jax"` gradients, individual scores and Hessians are exact
(automatic differentiation) and the likelihood is compiled. With `backend="numpy"` they are
numerical, as in Apollo (a port of `numDeriv`), which reproduces the standard errors of the R
runs most closely. A likelihood form without a JAX version falls back to NumPy with a
`mtuem.FallbackWarning`. JIT compilation takes one to three seconds, so on small models both
backends take about as long; on larger latent class models JAX is much faster.

**Bounds.** `beta.lower("Phi", 1e-10)` and `beta.upper("rho12", 1)` bound parameters (base
names apply to every class). The optimizer handles them with an adaptive barrier, as Apollo
does with inequality constraints, and finishes with an active-set step, so optima on a bound
are found; the EM steps respect them too.

**EM and multistart.** `method="em"` needs class likelihoods (a latent class model, not
`covariance="profiled"`). `MultiStart` draws `n` starting points uniformly within the start
values plus or minus one (`lower=` and `upper=` change the range of a parameter), optionally
runs a few EM iterations on each (`first_em=True`, `em_iter_max`), and keeps the best converged
optimum that is not a saddle point. `result.starts` holds the result of every candidate.

## Saving and loading

`result.save()` writes `<model_name>_output.txt`, `_estimates.csv`, the covariance and
correlation files, and `<model_name>_model.json` with the model, `Beta`, `Control`, settings and
estimates. `mtuem.load` restores the result:

```python
result = mtuem.load("output/maed-3eq_model.json")             # report, estimates, delta method
result = mtuem.load("output/maed-3eq_model.json", data=data)  # also values of time, predictions
```

A subclass of `MTUEM` is restored by its import path, or by `model_class=`. A `CustomModel`
is not saved: pass it as `model=` together with `data=`.

## Examples

`examples/maed_*.py` are the Python versions of the examples of the R package: a single
equation model, three-equation models with each error structure (profiled, `sig`, `sig` and
`rho`, Cholesky), and two-class latent class models with per-class, profiled and shared
covariances. Run them from `mtuem-py` with `uv run python examples/maed_1eq.py`.

## Datasets

`mtuem.maed()` and `mtuem.load_dataset(name)` for `maed`, `enut.i`, `enut.i.raw`, `enut.ii`,
`enut.ii.raw`, `enut.bigcities` and `enut.bigcities.raw` (`mtuem.DATASETS`), converted from the
R package's `.rda` files by `data-raw/convert_data.py`.

## Credit: Apollo

The R package estimates its models with [Apollo](https://www.apollochoicemodelling.com/)
(Hess and Palma, 2019, Journal of Choice Modelling 32, 100170). Apollo has no Python version,
so the estimation engine of this package reproduces the parts of Apollo these models use:
maximum likelihood with BFGS seeded by the BHHH matrix, classical and robust covariance
matrices, EM for latent class models, class allocation and posteriors, the delta method,
inequality constraints and the layout of the estimation report. The results match those of the
R package with Apollo (see below). Credit for the method goes to Stephane Hess and David Palma.

Where the engine deliberately differs from Apollo:

- Likelihoods work in logs, so a badly fitted observation cannot make the log-likelihood
  underflow to `-Inf` (R's mtuem returns densities and Apollo takes their log).
- When the optimizer stops at a saddle point (for example the symmetric solution of a latent
  class model started from identical classes), it steps along the direction of positive
  curvature and re-estimates, up to `saddle_escapes` times (default 3).
- Parameters are scaled by `max(|start|, 1)`, not `|start|`, so near-zero start values (common
  after EM) do not freeze those parameters.
- When the first step of the numerical Hessian leaves the region where the likelihood is
  finite, the Hessian is recomputed with relative steps of 1e-4.
- `LL(start)` is the log-likelihood at the start values and `BIC` uses the number of rows.
- EM maximizes over all parameters at once, so it accepts parameters shared across classes.
- Predictions work for latent class models.

## Why not Biogeme

Biogeme builds the likelihood one observation at a time from its expression language. Several MTUEM
specifications need statistics of the whole sample at every evaluation: the profiled covariance is
the covariance of the current residuals, and the profiled latent class model runs an inner fixed
point over all observations. Neither fits Biogeme's row-wise model, so the engine is built on
numpy, scipy and JAX instead, which also allows reproducing Apollo's output layout and standard
errors exactly.

## Agreement with the R package

`tests/test_r_parity.py` compares against results of the R package
(`tests/reference/make_reference.R`, R 4.4.3, Apollo 0.3.9, mtuem 1.2.2) for eight models:
log-likelihoods agree to about 1e-11 at fixed points, single-class estimates from the R start
values to about 1e-5 and standard errors to about 0.01% (NumPy backend). Values of time, delta
method results, residual correlations, predictions and profiled covariances match, and the
estimates table prints as Apollo prints it. From the R start values, a single latent class
estimation is path dependent in both languages, which is why the latent class examples use a
multistart.

## Tests

```bash
uv run pytest               # the default suite
uv run pytest --runslow     # also the slow tests: examples, replications and long fits
```
