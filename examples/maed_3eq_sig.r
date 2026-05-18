rm(list = ls())
library(apollo)
library(mtuem)
set.seed(42)

# =============================================================================
# Example 2b: Three Equations with Explicit Standard Deviations
# =============================================================================
# 3 equations: work (Tw), one free time (Tf1), one free good (Ef1).
# sig parameters are estimated (one per equation).
# Correlations are inferred from residuals.
# =============================================================================

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-sig",
  modelDescr = "Three equations with explicit sigma, inferred correlations",
  indivID = "PeID",
  outputDirectory = "output"
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

model <- apollo_estimate(apollo_beta, apollo_fixed,
                          apollo_probabilities, apollo_inputs)
apollo_modelOutput(model)

### Retrieve covariance matrix
corr_result <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs,
                               vars = c("Tw", "Tf1", "Ef1"))
cat("\n--- Covariance Matrix ---\n")
print(corr_result$covar)
cat("\n--- Correlation Matrix ---\n")
print(corr_result$corr)
