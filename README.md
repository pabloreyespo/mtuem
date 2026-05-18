
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mtuem

<!-- badges: start -->

<!-- badges: end -->

The **mtuem** package implements the likelihood functions for estimating
Microeconomic Time-Use-Expenditure Models (MTUEM) within the
[Apollo](https://www.apollochoicemodelling.com/) choice modelling
framework. The models are based on Jara-Diaz et al. (2008), which
maximizes a Cobb-Douglas utility function over time-use and goods
consumption subject to time, monetary, and technical constraints.

## Installation

``` r
install.packages("devtools")
devtools::install_github("https://github.com/pabloreyespo/mtuem")
```

## Model Formulation Guide

### The MTUEM Framework

The MTUEM derives optimal allocations of time and expenditure from
utility maximization. The model has two budget constraints:

- **Time budget**: Total time equals the time budget (typically 168
  hours per week), split into committed time (Tc), work time (Tw), and
  freely allocated time.
- **Monetary budget**: Total expenditure equals labor income minus
  committed expenses (Ec), split into freely allocated goods.

### Two Parameterizations

The package offers two equivalent parameterizations:

|  | Theta-Phi (`mtuem_likelihood`) | Alpha-Beta (`mtuem_ab_likelihood`) |
|----|----|----|
| Work elasticities | `Theta`, `Phi`, `thw` | `alpha`, `beta` |
| Identification | Fix `Theta = 1` | No fixed parameters needed |
| Use when | Direct interpretation of time value | Alternative formulation |

### Equations and the “Leave One Out” Rule

- **work_times** (required): Column name(s) for work time, typically
  `c("Tw")`.
- **free_times** (optional): Column names for freely allocated time
  activities. You must leave at least one out for identification. The
  omitted category is recovered as the residual from the time budget.
- **free_goods** (optional): Column names for freely allocated goods.
  You must leave at least one out for identification. The omitted
  category is recovered as the residual from the monetary budget.

### Error Structure

The likelihood assumes normally distributed errors. Three specification
options exist:

1.  **Inferred covariance** (no `sig`, no `rho`): The error covariance
    is inferred from residuals. Simplest approach, valid only for
    single-equation models.
2.  **Explicit sigma, inferred correlations** (`sig` specified, no
    `rho`): Standard deviations are estimated as parameters.
    Correlations are inferred from residuals. Valid for multi-equation
    models.
3.  **Full estimation** (both `sig` and `rho` specified): Both standard
    deviations and correlations are estimated as parameters.

### Retrieving the Covariance Matrix

After estimation, use
`mtuem_get_corr(model, apollo_probabilities, apollo_inputs, vars)` where
`vars` is a character vector of the time-use-expenditure column names.
Alternatively, call the likelihood with `functionality = "get_covar"`.

------------------------------------------------------------------------

## Example 1: Single Equation (Work Only)

The simplest model estimates only the work equation. No free times or
free goods are specified.

``` r
library(mtuem)
library(apollo)

### Initialise
apollo_initialise()
apollo_control <- list(
  modelName = "maed-1eq",
  modelDescr = "Single equation MTUEM (work only)",
  indivID = "PeID"
)

### Load data
database <- arrow::read_parquet("data/maed.parquet.gzip")
database <- database[database$EcI > 0, ]

### Parameters: Theta fixed to 1 for identification
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

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)

### Predictions
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)
```

------------------------------------------------------------------------

## Example 2: Three Equations with Different Error Structures

This example shows a model with 3 equations: work (Tw), one free time
(Tf1), and one free good (Ef1). We demonstrate three error structure
specifications.

### 2a: Inferred Covariance Matrix

No `sig` or `rho` parameters. The covariance is inferred from residuals.

``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-inferred",
  modelDescr = "Three equations with inferred covariance",
  indivID = "PeID"
)

database <- arrow::read_parquet("data/maed.parquet.gzip")
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

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)

### Retrieve covariance matrix after estimation
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs,
                               vars = c("Tw", "Tf1", "Ef1"))
corr_result$covar   # Covariance matrix
corr_result$corr    # Correlation matrix
```

### 2b: Explicit Standard Deviations, Inferred Correlations

Estimate `sig` parameters for each equation. Correlations are inferred
from residuals.

``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-sig",
  modelDescr = "Three equations with explicit sigma, inferred correlations",
  indivID = "PeID"
)

database <- arrow::read_parquet("data/maed.parquet.gzip")
database <- database[database$EcI > 0, ]

### Structural parameters + one sigma per equation
apollo_beta <- c(Theta = 1, Phi = 1, thw = 0,
                 theta_1 = 0.5, phi_1 = 0.5,
                 sig_1 = 5, sig_2 = 3, sig_3 = 10)
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

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)
```

### 2c: Full Estimation (Sigma and Rho)

Estimate both standard deviations and correlations. For 3 equations, the
upper triangle of the correlation matrix has 3 elements: rho\[1,2\],
rho\[1,3\], rho\[2,3\]. These must be ordered horizontally.

``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-full",
  modelDescr = "Three equations with full covariance estimation",
  indivID = "PeID"
)

database <- arrow::read_parquet("data/maed.parquet.gzip")
database <- database[database$EcI > 0, ]

### Structural + sigma + rho parameters
### rho ordering for 3 equations: rho_12, rho_13, rho_23 (horizontal upper triangle)
apollo_beta <- c(Theta = 1, Phi = 1, thw = 0,
                 theta_1 = 0.5, phi_1 = 0.5,
                 sig_1 = 5, sig_2 = 3, sig_3 = 10,
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

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)

### Retrieve the estimated covariance matrix
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs,
                               vars = c("Tw", "Tf1", "Ef1"))
print(corr_result$covar)
print(corr_result$corr)
```

------------------------------------------------------------------------

## Example 3: Latent Class Model (1 Class, 1 Equation)

This example shows the latent class framework with a single class and
one equation. This serves as a stepping stone to multi-class models.

``` r
library(mtuem)
library(apollo)

apollo_initialise()
apollo_control <- list(
  modelName = "maed-lc-1class",
  modelDescr = "Latent class model with 1 class and 1 equation",
  indivID = "PeID",
  noValidation = TRUE
)

database <- arrow::read_parquet("data/maed.parquet.gzip")
database <- database[database$EcI > 0, ]

### Parameters for class 1 only
apollo_beta <- c(
  Phi_1 = 1,
  thw_1 = 0,
  asc_1 = 0
)
apollo_fixed <- c("asc_1")

### Class membership parameters (required by apollo_lc even for 1 class)
apollo_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["Phi"]] <- list(Phi_1)
  lcpars[["thw"]] <- list(thw_1)

  V <- list()
  V[["class_1"]] <- asc_1

  classAlloc_settings <- list(
    classes = c(class_1 = 1),
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
      Theta = 1,
      Phi = Phi[[s]],
      thw = thw[[s]]
    )

    mtuem_settings <- list(
      work_times = c("Tw"),
      work_elasticities = work_elasticities,
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

### Estimate (use customMultiStart for latent class models)
source("mtuem/R/customMultiStart.R")
multi_starts <- customMultiStart(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  customMultistart_settings = list(
    apolloBetaMax = apollo_beta + 1,
    apolloBetaMin = apollo_beta - 1,
    nCandidates = 20
  ),
  estimate_settings = list(
    writeIter = FALSE,
    silent = TRUE,
    maxIterations = 500,
    estimationRoutine = "bfgs",
    hessianRoutine = "maxLik"
  ),
  first_em = TRUE,
  em_iter_max = 3,
  verbose = FALSE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c(),
    goods_elasticities = c(),
    sig = c(),
    rho = c()
  ),
  nClass = 1
)

apollo_modelOutput(multi_starts$best_model)
```

------------------------------------------------------------------------

## Post-Estimation

### Values of Time

After estimation, compute the Value of Leisure Time (VoL) and Value of
Time Assigned to Work (VTAW):

``` r
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)

Tw <- pred$Tw
w  <- database$w
Tc <- database$Tc
Ec <- database$Ec

cteVoL  <- mean((w * Tw - Ec) / (168 - Tw - Tc))
cteVTAW <- mean((w * Tw - Ec) / Tw)

delta <- apollo_deltaMethod(model, list(
  expression = c(
    VoL  = paste0(cteVoL,  "*Theta/Phi"),
    VTAW = paste0(cteVTAW, "*thw/Phi")
  )
))

delta
```

### Covariance Matrix Retrieval

``` r
### Using mtuem_get_corr (recommended)
result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs,
                          vars = c("Tw", "Tf1", "Ef1"))
result$covar
result$corr

### Using get_covar functionality directly
covar_result <- apollo_probabilities(apollo_beta, apollo_inputs, functionality = "get_covar")
```
