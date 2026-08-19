---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# mtuem

<!-- badges: start -->
<!-- badges: end -->

The **mtuem** package implements the likelihood functions for estimating Microeconomic Time-Use-Expenditure Models (MTUEM) within the [Apollo](https://www.apollochoicemodelling.com/) choice modelling framework. The models are based on Jara-Diaz et al. (2008), which maximizes a Cobb-Douglas utility function over time-use and goods consumption subject to time, monetary, and technical constraints.

## Installation

```r
install.packages("devtools")
devtools::install_github("https://github.com/pabloreyespo/mtuem")
```

## Model Formulation Guide

### The MTUEM Framework

The MTUEM derives optimal allocations of time and expenditure from utility maximization. The model has two budget constraints:

- **Time budget**: Total time equals the time budget (typically 168 hours per week), split into committed time (Tc), work time (Tw), and freely allocated time.
- **Monetary budget**: Total expenditure equals labor income minus committed expenses (Ec), split into freely allocated goods.

### Two Parameterizations

The package offers two equivalent parameterizations:

| | Theta-Phi (`mtuem_likelihood`) | Alpha-Beta (`mtuem_ab_likelihood`) |
|---|---|---|
| Work elasticities | `Theta`, `Phi`, `thw` | `alpha`, `beta` |
| Identification | Fix `Theta = 1` | No fixed parameters needed |
| Use when | Direct interpretation of time value | Alternative formulation |

### Equations and the "Leave One Out" Rule

- **work_times** (required): Column name(s) for work time, typically `c("Tw")`.
- **free_times** (optional): Column names for freely allocated time activities. You must leave at least one out for identification. The omitted category is recovered as the residual from the time budget.
- **free_goods** (optional): Column names for freely allocated goods. You must leave at least one out for identification. The omitted category is recovered as the residual from the monetary budget.

### Error Structure

The likelihood assumes normally distributed errors. Four specification options exist:

1. **Profiled covariance** (no `sig`, no `rho`): The error covariance is profiled out of the likelihood, recomputed from the current residuals at every evaluation. Simplest approach, valid only for single-equation models.
2. **Explicit sigma, profiled correlations** (`sig` specified, no `rho`): Standard deviations are estimated as parameters. Correlations are profiled out. Valid for multi-equation models.
3. **Full estimation** (both `sig` and `rho` specified): Both standard deviations and correlations are estimated as parameters.
4. **Cholesky parameterisation** (`cholesky` specified): The covariance is estimated through its lower-triangular Cholesky factor. Positive definite by construction and free of the bounds that constrain `rho`, so this is the **recommended** option for multi-equation models. Recover standard deviations and correlations with `mtuem_cholesky_deltaMethod()`.

### Retrieving the Covariance Matrix

`mtuem_get_corr(model, apollo_probabilities, apollo_inputs)` returns the covariance, correlation and standard deviations implied by the residuals. Under the Cholesky parameterisation the covariance is instead recovered from the estimated factor, either by rebuilding `L` and forming `L %*% t(L)`, or with `mtuem_cholesky_deltaMethod()` when standard errors are needed.

## Values of Time

The values of time follow from the estimates and the observed budgets,

$$VoL = \frac{w\,T_w - E_c}{\tau - T_w - T_c}\cdot\frac{\Theta}{\Phi}, \qquad VTAW = \frac{w\,T_w - E_c}{T_w}\cdot\frac{\theta_w}{\Phi}$$

`mtuem_values_of_time(model, apollo_inputs, ...)` reports them. It rebuilds the optimal work time from the estimates, averages each leading fraction over the sample, and applies `apollo_deltaMethod()` to the parameter ratios that remain. It prints one line per quantity with the mean and its confidence interval, for `VoL`, `VTAW` and the observed wage rate, and returns the same numbers in `$summary` alongside the per-observation values in `$individual`. Pass `silent = TRUE` to suppress the printed report.

For a latent class model, pass `nClass` and `apollo_probabilities`. Each class is averaged with its own class membership weights, so observations the model assigns elsewhere barely count towards it. The weights default to the posterior probabilities, recovered from `apollo_lcConditionals()` or, when `apollo_lc` was bypassed, from `mtuem_get_profiled_sigma()`; pass `weights` to use the prior shares or any other weighting instead.

`apollo_prediction()` also returns `VoL` and `VTAW` per individual alongside the optimal allocations, for single-class models.

---

## Example 1: Single Equation (Work Only)

The simplest model estimates only the work equation. No free times or free goods are specified. Each example first runs a single estimation from the declared starting values, then re-estimates with `customMultiStart()`, which searches over randomly drawn starting values and returns the best non-saddle optimum in `$best_model`, and closes with the values of time.


``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-1eq",
  modelDescr = "Single equation MTUEM (work only)",
  indivID = "PeID",
  outputDirectory = "output"
)

est_set <- list(
  writeIter = FALSE,
  maxIterations = 500,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

apollo_beta <- c(Theta = 1, Phi = 1, thw = 0)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)

  mtuem_settings <- list(
    work_times = c("Tw"),
    work_elasticities = work_elasticities,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: same model, best of several random starts
multi_start <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw")
  )
)
model <- multi_start$best_model
apollo_modelOutput(model)

### Values of time. mtuem_values_of_time rebuilds the optimal work time from the
### estimates and reports the sample means of
###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi,  equivalently VoL - w
### with delta method standard errors, printing one line per quantity with its
### confidence interval. vot$summary and vot$individual hold the same numbers.
vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")

### Predicted allocations
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

---

## Example 2: Three Equations with Different Error Structures

Three equations: work (Tw), one free time (Tf1) and one free good (Ef1), in four error structure specifications. The Cholesky specification (2d) is the recommended one.

### 2a: Profiled Covariance Matrix

No `sig` or `rho` parameters. The covariance is profiled out of the likelihood and recomputed from the current residuals at every evaluation. Cheapest in parameter count, but the error structure carries no standard errors.


``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-profiled",
  modelDescr = "Three equations with profiled covariance",
  indivID = "PeID",
  outputDirectory = "output"
)

est_set <- list(
  writeIter = FALSE,
  maxIterations = 500,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Only structural parameters; no sig or rho
apollo_beta <- c(Theta = 1, Phi = 1, thw = 0,
                 theta_1 = 0.5, phi_1 = 0.5)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("Tf1"),
    free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: same model, best of several random starts
multi_start <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("theta_1"),
    goods_elasticities = c("phi_1")
  )
)
model <- multi_start$best_model
apollo_modelOutput(model)

### Error covariance implied by the residuals
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs)
print(corr_result$covar)
print(corr_result$corr)

### Values of time. mtuem_values_of_time rebuilds the optimal work time from the
### estimates and reports the sample means of
###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi,  equivalently VoL - w
### with delta method standard errors, printing one line per quantity with its
### confidence interval. vot$summary and vot$individual hold the same numbers.
vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")

### Predicted allocations
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

### 2b: Explicit Standard Deviations, Profiled Correlations

One estimated `sig` per equation. Correlations are profiled out.


``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-sig",
  modelDescr = "Three equations with explicit sigma, profiled correlations",
  indivID = "PeID",
  outputDirectory = "output"
)

est_set <- list(
  writeIter = FALSE,
  maxIterations = 500,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Structural parameters + one sigma per equation
apollo_beta <- c(Theta = 1, Phi = 1, thw = 0,
                 theta_1 = 0.5, phi_1 = 0.5,
                 sig_1 = 100, sig_2 = 100, sig_3 = 100)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)
  sig <- list(sig_1 = sig_1, sig_2 = sig_2, sig_3 = sig_3)

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("Tf1"),
    free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    sig = sig,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: same model, best of several random starts
multi_start <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("theta_1"),
    goods_elasticities = c("phi_1"),
    sig = c("sig_1", "sig_2", "sig_3")
  )
)
model <- multi_start$best_model
apollo_modelOutput(model)

### Error covariance implied by the residuals
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs)
print(corr_result$covar)
print(corr_result$corr)

### Values of time. mtuem_values_of_time rebuilds the optimal work time from the
### estimates and reports the sample means of
###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi,  equivalently VoL - w
### with delta method standard errors, printing one line per quantity with its
### confidence interval. vot$summary and vot$individual hold the same numbers.
vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")

### Predicted allocations
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

### 2c: Full Estimation (Sigma and Rho)

Both `sig` and `rho` are estimated. For 3 equations the `rho` parameters run across the upper triangle: `rho_12`, `rho_13`, `rho_23`.


``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-full",
  modelDescr = "Three equations with full covariance estimation",
  indivID = "PeID",
  outputDirectory = "output"
)

est_set <- list(
  writeIter = FALSE,
  maxIterations = 500,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Structural + sigma + rho parameters
### rho ordering for 3 equations: rho_12, rho_13, rho_23 (horizontal upper triangle)
apollo_beta <- c(
  Theta = 1, 
  Phi = 1, 
  thw = 0,
  theta_1 = 0.5, 
  phi_1 = 0.5,
  sig_1 = 100, 
  sig_2 = 100, 
  sig_3 = 100,
  rho_12 = 0, rho_13 = 0, rho_23 = 0)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)
  sig <- list(sig_1 = sig_1, sig_2 = sig_2, sig_3 = sig_3)
  rho <- list(rho_12 = rho_12, rho_13 = rho_13, rho_23 = rho_23)

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("Tf1"),
    free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    sig = sig,
    rho = rho,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: same model, best of several random starts
multi_start <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("theta_1"),
    goods_elasticities = c("phi_1"),
    sig = c("sig_1", "sig_2", "sig_3"),
    rho = c("rho_12", "rho_13", "rho_23")
  )
)
model <- multi_start$best_model
apollo_modelOutput(model)

### Error covariance implied by the residuals
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs)
print(corr_result$covar)
print(corr_result$corr)

### Values of time. mtuem_values_of_time rebuilds the optimal work time from the
### estimates and reports the sample means of
###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi,  equivalently VoL - w
### with delta method standard errors, printing one line per quantity with its
### confidence interval. vot$summary and vot$individual hold the same numbers.
vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")

### Predicted allocations
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

### 2d: Cholesky Parameterisation (recommended)

The covariance is estimated through its lower-triangular Cholesky factor (`chk11`, `chk21`, `chk22`, `chk31`, `chk32`, `chk33`), given row-wise over the lower triangle. This keeps the matrix positive definite by construction and avoids the bounds that constrain `rho` in 2c.

The estimates are the entries of the factor, not the covariance itself, so they have to be decoded. Rebuilding `L` and forming `Sigma = L %*% t(L)` gives the covariance, `sqrt(diag(Sigma))` the standard deviations and `cov2cor(Sigma)` the correlations; `mtuem_cholesky_deltaMethod()` returns the same quantities with standard errors. Both are shown in the script.


``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-cholesky",
  modelDescr = "Three equations with Cholesky covariance estimation",
  indivID = "PeID",
  outputDirectory = "output"
)

est_set <- list(
  writeIter = FALSE,
  maxIterations = 500,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Cholesky parameters are given row-wise over the lower triangle
apollo_beta <- c(
  Theta = 1,
  Phi = 1,
  thw = 0,
  theta_1 = 0.5,
  phi_1 = 0.5,
  chk11 = 100,
  chk21 = 0, chk22 = 100,
  chk31 = 0, chk32 = 0, chk33 = 100
)
apollo_fixed <- c("Theta")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(Theta = Theta, Phi = Phi, thw = thw)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)
  chk <- list(chk11 = chk11, chk21 = chk21, chk22 = chk22,
              chk31 = chk31, chk32 = chk32, chk33 = chk33)

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("Tf1"),
    free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    cholesky = chk,
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

cholesky_names <- c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33")
eq_names <- c("Tw", "Tf1", "Ef1")

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: same model, best of several random starts
multi_start <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("theta_1"),
    goods_elasticities = c("phi_1"),
    cholesky = cholesky_names
  )
)
model <- multi_start$best_model
apollo_modelOutput(model)

# -----------------------------------------------------------------------------
# Decoding the Cholesky factor
# -----------------------------------------------------------------------------
# The estimates chk11...chk33 are the non-zero entries of the lower triangular
# L, read row by row. The covariance is Sigma = L %*% t(L), the standard
# deviations are sqrt(diag(Sigma)) and the correlations are cov2cor(Sigma).

cholesky_to_cov <- function(estimates, cholesky_names) {
  M <- (-1 + sqrt(1 + 8 * length(cholesky_names))) / 2
  L <- matrix(0, M, M)
  L[cbind(rep(1:M, 1:M), sequence(1:M))] <- estimates[cholesky_names]
  L %*% t(L)
}

L_entries <- model$estimate[cholesky_names]
Sigma <- cholesky_to_cov(L_entries, cholesky_names)
dimnames(Sigma) <- list(eq_names, eq_names)

print(Sigma)
print(sqrt(diag(Sigma)))
print(cov2cor(Sigma))

### Same quantities with standard errors, via the delta method
print(mtuem_cholesky_deltaMethod(model, cholesky_names, eq_names))

### Values of time. mtuem_values_of_time rebuilds the optimal work time from the
### estimates and reports the sample means of
###   VoL  = ((w * Tw - Ec) / (tau - Tw - Tc)) * Theta / Phi
###   VTAW = ((w * Tw - Ec) / Tw) * thw / Phi,  equivalently VoL - w
### with delta method standard errors, printing one line per quantity with its
### confidence interval. vot$summary and vot$individual hold the same numbers.
vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")

### Predicted allocations
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

---

## Example 3: Latent Class Model (2 Classes, 3 Equations each)

A latent class model with 2 classes, each with a work equation (Tw), a free time equation (Tf1) and a free good equation (Ef1). Class membership is modelled with covariates. The three variants differ in how much the error covariance costs in parameters.

| Variant | Covariance parameters | Single estimation |
|---|---|---|
| 3a Cholesky | 12 (6 per class) | `apollo_lcEM()` |
| 3b Profiled | 0 | `apollo_estimate()`, `apollo_lc` is bypassed |
| 3c Shared | 6 | `apollo_estimate()`, shared parameters block `apollo_lcEM` |

`apollo_prediction()` does not work for these models: it looks for `inClassProb[[1]]$chosen`, which a continuous density component does not provide. Call `mtuem_likelihood(..., functionality = "prediction")` per class instead, or rebuild the allocations and values of time from the estimates as each script does at the end.

### 3a: Per-Class Cholesky Covariance

Estimated with an EM warm start through `apollo_lcEM()`, then with `customMultiStart(..., first_em = TRUE)`, which warm starts each candidate the same way. The per-class Cholesky factors are decoded class by class at the bottom.


``` r
library(mtuem)
library(apollo)

est_set <- list(
  writeIter = FALSE,
  silent = FALSE,
  maxIterations = 500,
  scaleHessian = FALSE,
  scaleAfterConvergence = FALSE,
  hessianRoutine = "numDeriv"
)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-lc-2class-3eq",
  modelDescr = "Latent class model with 2 classes and 3 equations each",
  indivID = "PeID",
  outputDirectory = "output",
  noValidation = TRUE,
  noDiagnostics = TRUE
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

apollo_beta <- c(
  # Class 1 structural
  Theta_1 = 1,
  Phi_1 = 1,
  thw_1 = 0,
  th1_1 = 0.3,
  ph1_1 = 0.3,
  # Class 1 error structure
  chk11_1 = 10,
  chk21_1 = 0, chk22_1 = 20,
  chk31_1 = 0, chk32_1 = 0, chk33_1 = 20,
  # Class 1 membership
  asc_1 = 0,
  x_female_1 = 0,
  x_older45_1 = 0,
  # Class 2 structural
  Theta_2 = 1,
  Phi_2 = 1,
  thw_2 = 0,
  th1_2 = 0.3,
  ph1_2 = 0.3,
  # Class 2 error structure
  chk11_2 = 10,
  chk21_2 = 0, chk22_2 = 20,
  chk31_2 = 0, chk32_2 = 0, chk33_2 = 20,
  # Class 2 membership
  asc_2 = 0,
  x_female_2 = 0,
  x_older45_2 = 0
)

apollo_fixed <- c(
  "Theta_1", "Theta_2",
  "asc_1", "x_female_1", "x_older45_1"
)

### Class membership parameters
apollo_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["Theta"]] <- list(Theta_1, Theta_2)
  lcpars[["Phi"]] <- list(Phi_1, Phi_2)
  lcpars[["thw"]] <- list(thw_1, thw_2)
  lcpars[["th1"]] <- list(th1_1, th1_2)
  lcpars[["ph1"]] <- list(ph1_1, ph1_2)
  lcpars[["chk11"]] <- list(chk11_1, chk11_2)
  lcpars[["chk21"]] <- list(chk21_1, chk21_2)
  lcpars[["chk22"]] <- list(chk22_1, chk22_2)
  lcpars[["chk31"]] <- list(chk31_1, chk31_2)
  lcpars[["chk32"]] <- list(chk32_1, chk32_2)
  lcpars[["chk33"]] <- list(chk33_1, chk33_2)

  V <- list()
  V[["class_1"]] <- asc_1 + (female * x_female_1) + (older45 * x_older45_1)
  V[["class_2"]] <- asc_2 + (female * x_female_2) + (older45 * x_older45_2)

  classAlloc_settings <- list(
    classes = c(class_1 = 1, class_2 = 2),
    utilities = V
  )
  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  for (s in 1:length(pi_values)) {
    work_elasticities <- list(
      Theta = Theta[[s]],
      Phi = Phi[[s]],
      thw = thw[[s]]
    )

    times_elasticities <- list(
      th1 = th1[[s]]
    )

    goods_elasticities <- list(
      ph1 = ph1[[s]]
    )

    chk <- list(
      chk11 = chk11[[s]],
      chk21 = chk21[[s]],
      chk22 = chk22[[s]],
      chk31 = chk31[[s]],
      chk32 = chk32[[s]],
      chk33 = chk33[[s]]
    )

    mtuem_settings <- list(
      work_times = c("Tw"),
      free_times = c("Tf1"),
      free_goods  = c("Ef1"),
      work_elasticities = work_elasticities,
      times_elasticities = times_elasticities,
      goods_elasticities = goods_elasticities,
      cholesky = chk,
      Tc = Tc, Ec = EcI, w = w, tau = 168,
      componentName = paste0("class_", s)
    )

    P[[paste0("class_", s)]] <- mtuem_likelihood(mtuem_settings, functionality)
  }

  lc_settings <- list(inClassProb = P, classProb = pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

cholesky_names_1 <- c("chk11_1", "chk21_1", "chk22_1", "chk31_1", "chk32_1", "chk33_1")
cholesky_names_2 <- c("chk11_2", "chk21_2", "chk22_2", "chk31_2", "chk32_2", "chk33_2")
eq_names <- c("Tw", "Tf1", "Ef1")

### Single estimation with an EM warm start
model <- apollo_lcEM(apollo_beta, apollo_fixed,
                     apollo_probabilities, apollo_inputs,
                     lcEM_settings = list(EMmaxIterations = 10),
                     estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: best of several random starts, each warm started
### with apollo_lcEM (first_em = TRUE)
multi_starts <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  first_em = TRUE,
  em_iter_max = 10,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("th1"),
    goods_elasticities = c("ph1"),
    cholesky = c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33"),
    covs = c("asc", "x_female", "x_older45")
  ),
  nClass = 2
)

### best_model is NA when no candidate reaches a non-saddle optimum
if (!is.list(multi_starts$best_model)) {
  stop("No candidate converged to a non-saddle optimum. Increase nCandidates, ",
       "widen apolloBetaMin/apolloBetaMax, or set non_saddle = FALSE to inspect ",
       "the saddle solutions in multi_starts$models.")
}
model <- multi_starts$best_model
apollo_modelOutput(model)

# -----------------------------------------------------------------------------
# Decoding the Cholesky factor
# -----------------------------------------------------------------------------
# The estimates chk11...chk33 are the non-zero entries of the lower triangular
# L, read row by row. The covariance is Sigma = L %*% t(L), the standard
# deviations are sqrt(diag(Sigma)) and the correlations are cov2cor(Sigma).

cholesky_to_cov <- function(estimates, cholesky_names) {
  M <- (-1 + sqrt(1 + 8 * length(cholesky_names))) / 2
  L <- matrix(0, M, M)
  L[cbind(rep(1:M, 1:M), sequence(1:M))] <- estimates[cholesky_names]
  L %*% t(L)
}

for (s in 1:2) {
  cholesky_names <- get(paste0("cholesky_names_", s))
  Sigma <- cholesky_to_cov(model$estimate[cholesky_names], cholesky_names)
  dimnames(Sigma) <- list(eq_names, eq_names)
  cat("
--- Class", s, "---
")
  print(Sigma)
  print(sqrt(diag(Sigma)))
  print(cov2cor(Sigma))
  ### Same quantities with standard errors
  print(mtuem_cholesky_deltaMethod(model, cholesky_names, eq_names))
}

### apollo_prediction does not work for these latent class models: it looks for
### inClassProb[[1]]$chosen, which a continuous density component does not
### provide. Call mtuem_likelihood(..., functionality = "prediction") per class
### instead.

### Values of time by class. apollo_prediction cannot be used, so
### mtuem_values_of_time rebuilds the optimal work time from the estimates and
### averages the values of time with the posterior class membership weights,
### which it recovers from the model. The delta method is applied per class and
### the means are printed with their confidence intervals.
vot <- mtuem_values_of_time(model, apollo_inputs, apollo_probabilities,
                            nClass = 2, Ec = "EcI")
```

### 3b: Profiled Covariance

`mtuem_lc_profiled_likelihood()` concentrates the class-specific covariance out of the likelihood with an inner fixed point, so `apollo_beta` holds only structural and class membership parameters. It already returns the class mixture, so `apollo_lc` is bypassed and `apollo_lcEM` has no class structure to attach to: estimation goes through `apollo_estimate()` and `first_em` must be `FALSE`. Use `mtuem_profile_ic()` so the information criteria count the profiled covariance parameters.


``` r
library(mtuem)
library(apollo)

est_set <- list(
  writeIter = FALSE,
  silent = FALSE,
  maxIterations = 500,
  scaleHessian = FALSE,
  scaleAfterConvergence = FALSE,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-lc-2class-profiled-cov",
  modelDescr = "Latent class model with 2 classes and a profiled-out error covariance",
  indivID = "PeID",
  outputDirectory = "output",
  noValidation = TRUE,
  noDiagnostics = TRUE
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

apollo_beta <- c(
  # Class 1 structural
  Theta_1 = 1,
  Phi_1 = 1,
  thw_1 = 0,
  th1_1 = 0.3,
  ph1_1 = 0.3,
  # Class 1 membership
  asc_1 = 0,
  x_female_1 = 0,
  x_older45_1 = 0,
  # Class 2 structural
  Theta_2 = 1,
  Phi_2 = 1,
  thw_2 = 0,
  th1_2 = 0.3,
  ph1_2 = 0.3,
  # Class 2 membership
  asc_2 = 0,
  x_female_2 = 0,
  x_older45_2 = 0
)

apollo_fixed <- c(
  "Theta_1", "Theta_2",
  "asc_1", "x_female_1", "x_older45_1"
)

apollo_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["Theta"]] <- list(Theta_1, Theta_2)
  lcpars[["Phi"]] <- list(Phi_1, Phi_2)
  lcpars[["thw"]] <- list(thw_1, thw_2)
  lcpars[["th1"]] <- list(th1_1, th1_2)
  lcpars[["ph1"]] <- list(ph1_1, ph1_2)

  V <- list()
  V[["class_1"]] <- asc_1 + (female * x_female_1) + (older45 * x_older45_1)
  V[["class_2"]] <- asc_2 + (female * x_female_2) + (older45 * x_older45_2)

  classAlloc_settings <- list(
    classes = c(class_1 = 1, class_2 = 2),
    utilities = V
  )
  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  class_settings <- list()
  for (s in 1:length(pi_values)) {
    work_elasticities <- list(
      Theta = Theta[[s]],
      Phi = Phi[[s]],
      thw = thw[[s]]
    )

    times_elasticities <- list(
      th1 = th1[[s]]
    )

    goods_elasticities <- list(
      ph1 = ph1[[s]]
    )

    class_settings[[s]] <- list(
      work_times = c("Tw"),
      free_times = c("Tf1"),
      free_goods = c("Ef1"),
      goods_cost = list(Ef1 = 1),
      work_elasticities = work_elasticities,
      times_elasticities = times_elasticities,
      goods_elasticities = goods_elasticities,
      Tc = Tc, Ec = EcI, w = w, tau = 168,
      componentName = paste0("class_", s)
    )
  }

  lc_settings <- list(
    class_settings = class_settings,
    classProb = pi_values
  )

  P[["model"]] <- mtuem_lc_profiled_likelihood(lc_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Single estimation from the starting values above
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: best of several random starts (first_em = FALSE
### because apollo_lc is bypassed, so there is no EM structure to warm start)
multi_starts <- customMultiStart(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  first_em = FALSE,
  verbose = TRUE,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("th1"),
    goods_elasticities = c("ph1"),
    covs = c("asc", "x_female", "x_older45")
  ),
  nClass = 2
)

### best_model is NA when no candidate reaches a non-saddle optimum
if (!is.list(multi_starts$best_model)) {
  stop("No candidate converged to a non-saddle optimum. Increase nCandidates, ",
       "widen apolloBetaMin/apolloBetaMax, or set non_saddle = FALSE to inspect ",
       "the saddle solutions in multi_starts$models.")
}
model <- multi_starts$best_model
apollo_modelOutput(model)

### Re-evaluate at the converged estimates so the profiled covariance matches
### the final parameter values before retrieving it
invisible(apollo_probabilities(model$estimate, apollo_inputs, functionality = "estimate"))

profiled <- mtuem_get_profiled_sigma()
cat("\n--- Profiled-out error covariance (Class 1) ---\n")
cat("sig:\n"); print(profiled$sig[[1]])
cat("rho:\n"); print(profiled$rho[[1]])
cat("\n--- Profiled-out error covariance (Class 2) ---\n")
cat("sig:\n"); print(profiled$sig[[2]])
cat("rho:\n"); print(profiled$rho[[2]])
cat("\nClass shares:\n")
print(profiled$shares)

### AIC/BIC must count the profiled-out covariance parameters, which Apollo
### never sees: nClass * nEq * (nEq + 1) / 2
print(mtuem_profile_ic(model, nClass = 2, nEq = 3))

### Values of time by class. apollo_prediction cannot be used, so
### mtuem_values_of_time rebuilds the optimal work time from the estimates and
### averages the values of time with the posterior class membership weights,
### which it recovers from the model. The delta method is applied per class and
### the means are printed with their confidence intervals.
vot <- mtuem_values_of_time(model, apollo_inputs, apollo_probabilities,
                            nClass = 2, Ec = "EcI")
```

### 3c: Shared Covariance

Both classes share a single lower-Cholesky factor. The shared parameters are declared without a class suffix in `apollo_beta` and are not listed in `apollo_lcPars`, so the class loop references them without indexing by class. `apollo_lcEM()` cannot be used because it does not handle parameters shared across classes.


``` r
library(mtuem)
library(apollo)

est_set <- list(
  writeIter = FALSE,
  silent = FALSE,
  maxIterations = 500,
  scaleHessian = FALSE,
  scaleAfterConvergence = FALSE,
  estimationRoutine = "bgw",
  hessianRoutine = "numDeriv"
)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-lc-2class-shared-cov",
  modelDescr = "Latent class model with 2 classes sharing a single error covariance",
  indivID = "PeID",
  outputDirectory = "output",
  noValidation = TRUE,
  noDiagnostics = TRUE
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Parameters for 2 classes, with a SINGLE shared covariance structure
apollo_beta <- c(
  # Class 1 structural
  Theta_1 = 1,
  Phi_1 = 1,
  thw_1 = 0,
  th1_1 = 0.3,
  ph1_1 = 0.3,
  # Class 1 membership
  asc_1 = 0,
  x_female_1 = 0,
  x_older45_1 = 0,
  # Class 2 structural
  Theta_2 = 1,
  Phi_2 = 1,
  thw_2 = 0,
  th1_2 = 0.3,
  ph1_2 = 0.3,
  # Class 2 membership
  asc_2 = 0,
  x_female_2 = 0,
  x_older45_2 = 0,
  # Shared error structure
  chk11 = 6,
  chk21 = 0, chk22 = 7,
  chk31 = 0, chk32 = 0, chk33 = 30
)

apollo_fixed <- c(
  "Theta_1", "Theta_2",
  "asc_1", "x_female_1", "x_older45_1"
)

### Class membership parameters
apollo_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["Theta"]] <- list(Theta_1, Theta_2)
  lcpars[["Phi"]] <- list(Phi_1, Phi_2)
  lcpars[["thw"]] <- list(thw_1, thw_2)
  lcpars[["th1"]] <- list(th1_1, th1_2)
  lcpars[["ph1"]] <- list(ph1_1, ph1_2)

  V <- list()
  V[["class_1"]] <- asc_1 + (female * x_female_1) + (older45 * x_older45_1)
  V[["class_2"]] <- asc_2 + (female * x_female_2) + (older45 * x_older45_2)

  classAlloc_settings <- list(
    classes = c(class_1 = 1, class_2 = 2),
    utilities = V
  )
  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  chk <- list(
    chk11 = chk11,
    chk21 = chk21,
    chk22 = chk22,
    chk31 = chk31,
    chk32 = chk32,
    chk33 = chk33
  )

  for (s in 1:length(pi_values)) {
    work_elasticities <- list(
      Theta = Theta[[s]],
      Phi = Phi[[s]],
      thw = thw[[s]]
    )

    times_elasticities <- list(
      th1 = th1[[s]]
    )

    goods_elasticities <- list(
      ph1 = ph1[[s]]
    )

    mtuem_settings <- list(
      work_times = c("Tw"),
      free_times = c("Tf1"),
      free_goods = c("Ef1"),
      goods_cost = list(Ef1 = 1),
      work_elasticities = work_elasticities,
      times_elasticities = times_elasticities,
      goods_elasticities = goods_elasticities,
      cholesky = chk,
      Tc = Tc, Ec = EcI, w = w, tau = 168,
      componentName = paste0("class_", s)
    )

    P[[paste0("class_", s)]] <- mtuem_likelihood(mtuem_settings, functionality)
  }

  lc_settings <- list(inClassProb = P, classProb = pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

cholesky_names <- c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33")
eq_names <- c("Tw", "Tf1", "Ef1")

### Single estimation from the starting values above. apollo_lcEM cannot be
### used here: it does not handle parameters shared across classes.
model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         estimate_settings = est_set)
apollo_modelOutput(model)

### Starting values module: best of several random starts
multi_starts <- customMultiStart(
  apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 5
  ),
  estimate_settings = est_set,
  first_em = FALSE,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("th1"),
    goods_elasticities = c("ph1"),
    cholesky = cholesky_names,
    covs = c("asc", "x_female", "x_older45")
  ),
  nClass = 2
)

### best_model is NA when no candidate reaches a non-saddle optimum
if (!is.list(multi_starts$best_model)) {
  stop("No candidate converged to a non-saddle optimum. Increase nCandidates, ",
       "widen apolloBetaMin/apolloBetaMax, or set non_saddle = FALSE to inspect ",
       "the saddle solutions in multi_starts$models.")
}
model <- multi_starts$best_model
apollo_modelOutput(model)

# -----------------------------------------------------------------------------
# Decoding the Cholesky factor
# -----------------------------------------------------------------------------
# The estimates chk11...chk33 are the non-zero entries of the lower triangular
# L, read row by row. The covariance is Sigma = L %*% t(L), the standard
# deviations are sqrt(diag(Sigma)) and the correlations are cov2cor(Sigma).

cholesky_to_cov <- function(estimates, cholesky_names) {
  M <- (-1 + sqrt(1 + 8 * length(cholesky_names))) / 2
  L <- matrix(0, M, M)
  L[cbind(rep(1:M, 1:M), sequence(1:M))] <- estimates[cholesky_names]
  L %*% t(L)
}

Sigma <- cholesky_to_cov(model$estimate[cholesky_names], cholesky_names)
dimnames(Sigma) <- list(eq_names, eq_names)

print(Sigma)
print(sqrt(diag(Sigma)))
print(cov2cor(Sigma))

### Same quantities with standard errors
print(mtuem_cholesky_deltaMethod(model, cholesky_names, eq_names))

### apollo_prediction does not work for these latent class models: it looks for
### inClassProb[[1]]$chosen, which a continuous density component does not
### provide. Call mtuem_likelihood(..., functionality = "prediction") per class
### instead.

### Values of time by class. apollo_prediction cannot be used, so
### mtuem_values_of_time rebuilds the optimal work time from the estimates and
### averages the values of time with the posterior class membership weights,
### which it recovers from the model. The delta method is applied per class and
### the means are printed with their confidence intervals.
vot <- mtuem_values_of_time(model, apollo_inputs, apollo_probabilities,
                            nClass = 2, Ec = "EcI")
```
