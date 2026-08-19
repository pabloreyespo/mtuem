rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# Example 2a: three equations (Tw, Tf1, Ef1) with the covariance profiled out.
# No sig or rho parameters: the covariance is recomputed from the residuals at
# every likelihood evaluation, so it carries no standard errors. Prefer
# maed_3eq_cholesky.r when the error structure itself is of interest.

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
