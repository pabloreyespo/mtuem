rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# Example 2d (recommended): three equations (Tw, Tf1, Ef1) with the covariance
# estimated through its lower-triangular Cholesky factor. Positive definite by
# construction and free of the bounds that constrain rho.

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
