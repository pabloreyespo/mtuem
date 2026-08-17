rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# =============================================================================
# Example: Latent Class Model (2 Classes, 3 Equations each)
# =============================================================================
# This example shows a latent class model with 2 classes, each with:
#   - 1 work equation (Tw)
#   - 2 free time equations (Tf1, Tf2)
#   - Full covariance estimation (sig + rho per class)
#   - Class membership modeled with covariates (female, older45)
# Uses customMultiStart for EM initialization and multi-start optimization.
# =============================================================================

est_set <- list(
  writeIter = FALSE,
  silent = FALSE,
  maxIterations = 500,
  scaleHessian = FALSE,
  scaleAfterConvergence = FALSE,
  estimationRoutine = "bgw",
  hessianRoutine = "maxLik"
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

### Parameters for 2 classes
### Class 1: Theta_1 (fixed), Phi_1, thw_1, th1_1, th2_1, sig/rho params
### Class 2: Theta_2 (fixed), Phi_2, thw_2, th1_2, th2_2, sig/rho params
### Class membership: asc_1 (fixed), x_female_1 (fixed), x_older45_1 (fixed)
apollo_beta <- c(
  # Class 1 structural
  Theta_1 = 1,
  Phi_1 = 1,
  thw_1 = 0,
  th1_1 = 0.3,
  th2_1 = 0.3,
  # Class 1 error structure (3 equations: Tw, Tf1, Tf2)
  chk11_1 = 10,
  chk21_1 = 0, chk22_1 = 5,
  chk31_1 = 0, chk32_1 = 0, chk33_1 = 5,
  # Class 1 membership
  asc_1 = 0,
  x_female_1 = 0,
  x_older45_1 = 0,
  # Class 2 structural
  Theta_2 = 1,
  Phi_2 = 1,
  thw_2 = 0,
  th1_2 = 0.3,
  th2_2 = 0.3,
  # Class 2 error structure
  chk11_2 = 10,
  chk21_2 = 0, chk22_2 = 5,
  chk31_2 = 0, chk32_2 = 0, chk33_2 = 5,
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
  lcpars[["th2"]] <- list(th2_1, th2_2)
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
      th1 = th1[[s]],
      th2 = th2[[s]]
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
      free_times = c("Tf1", "Tf2"),
      work_elasticities = work_elasticities,
      times_elasticities = times_elasticities,
      cholesky = chk,
      Tc = Tc, Ec = EcI, w = w, tau = 168,
      componentName = paste0("class_", s)
    )

    # Direct assignment: apollo_lc accepts numeric vectors for single-component classes
    P[[paste0("class_", s)]] <- mtuem_likelihood(mtuem_settings, functionality)
  }

  lc_settings <- list(inClassProb = P, classProb = pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

### Estimate using customMultiStart
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
  estimate_settings = est_set,
  first_em = TRUE,
  em_iter_max = 10,
  verbose = TRUE,
  non_saddle = TRUE,
  normalization = list(
    normalization = "theta_phi",
    work_elasticities = c("Theta", "Phi", "thw"),
    times_elasticities = c("th1", "th2"),
    goods_elasticities = c(),
    cholesky = c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33"),
    covs = c("asc", "x_female", "x_older45")
  ),
  nClass = 2
)

apollo_modelOutput(multi_starts$best_model)

### Predictions
pred <- apollo_prediction(multi_starts$best_model, apollo_probabilities, apollo_inputs)
head(pred[["model"]])

### Retrieve standard deviations and correlations for Class 1 and Class 2 via Delta Method
cat("\n--- Delta Method Results for Class 1 Sigmas and Rhos ---\n")
cholesky_names_1 <- c("chk11_1", "chk21_1", "chk22_1", "chk31_1", "chk32_1", "chk33_1")
eq_names_1 <- c("Tw_c1", "Tf1_c1", "Tf2_c1")
delta_results_1 <- mtuem_cholesky_deltaMethod(multi_starts$best_model, cholesky_names_1, eq_names_1)
print(delta_results_1)

cat("\n--- Delta Method Results for Class 2 Sigmas and Rhos ---\n")
cholesky_names_2 <- c("chk11_2", "chk21_2", "chk22_2", "chk31_2", "chk32_2", "chk33_2")
eq_names_2 <- c("Tw_c2", "Tf1_c2", "Tf2_c2")
delta_results_2 <- mtuem_cholesky_deltaMethod(multi_starts$best_model, cholesky_names_2, eq_names_2)
print(delta_results_2)
