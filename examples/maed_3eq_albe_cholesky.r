rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# =============================================================================
# Example 2c: Three Equations with Full Covariance Estimation
# =============================================================================
# 3 equations: work (Tw), one free time (Tf1), one free good (Ef1).
# Both sig and rho parameters are estimated.
# For 3 equations, rho parameters are ordered horizontally across the
# upper triangle: rho_12, rho_13, rho_23.
# =============================================================================

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-full",
  modelDescr = "Three equations with full covariance estimation",
  indivID = "PeID",
  outputDirectory = "output"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### Structural + sigma + rho parameters
### rho ordering for 3 equations: rho_12, rho_13, rho_23 (horizontal upper triangle)
apollo_beta <- c(
  alpha = 0.25,
  beta = 0.25,
  theta_1 = 0.125, 
  phi_1 = 0.125,
  chk11 = 100,
  chk21 = 0, chk22 = 100,
  chk31 = 0, chk32 = 0, chk33 = 100
  )
apollo_fixed <- c()

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  work_elasticities <- list(alpha = alpha, beta = beta)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)
  chk <- list(chk11 = chk11, chk21 = chk21, chk22 = chk22, chk31 = chk31, chk32 = chk32, chk33 = chk33)

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

  P[["model"]] <- mtuem_ab_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)

### Retrieve the estimated covariance matrix
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs,
                               vars = c("Tw", "Tf1", "Ef1"))
cat("\n--- Covariance Matrix ---\n")
print(corr_result$covar)
cat("\n--- Correlation Matrix ---\n")
print(corr_result$corr)
