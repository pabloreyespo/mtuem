rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# Latent class model (2 classes, 3 equations each: Tw, Tf1, Ef1) with the error
# covariance PROFILED OUT: mtuem_lc_profiled_likelihood concentrates the
# class-specific covariance out of the likelihood with an inner fixed point, so
# apollo_beta holds only structural and class membership parameters (0 covariance
# parameters instead of the 12 in maed_lc_2class_cholesky.r).
#
# The profiled likelihood already returns the class mixture, so apollo_lc is
# bypassed and apollo_lcEM cannot be used: estimation goes through
# apollo_estimate and customMultiStart runs with first_em = FALSE.

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
