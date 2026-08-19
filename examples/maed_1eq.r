rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# Example 1: single equation (work only). Theta is fixed to 1 for identification.

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
