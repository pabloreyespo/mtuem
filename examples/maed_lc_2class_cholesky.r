rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# Latent class model, 2 classes, 3 equations each (Tw, Tf1, Ef1), with a full
# per-class covariance parameterised through its lower-Cholesky factor:
# 6 covariance parameters per class, 12 in total. Class membership uses
# covariates (female, older45).
#
# Reference specification. maed_lc_2class_shared.r shares one Cholesky factor
# across classes (6 parameters) and maed_lc_2class_profiled.r concentrates the
# covariance out entirely (0 parameters).

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
vot <- mtuem_values_zof_time(model, apollo_inputs, apollo_probabilities,
                            nClass = 2, Ec = "EcI")
