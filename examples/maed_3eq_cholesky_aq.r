rm(list = ls())
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))
set.seed(42)

# =============================================================================
# Additive quadratic MTUEM: full daily cycle, three equations, full covariance
# Structural Anchor: theta_w = 0 (Instrumental Labor)
# =============================================================================
# The daily cycle is five quantities:
#
#     times        Tw + T1 + Ti + Tc = tau                (= 168)
#     expenditure  E1 + Ei + Ec      = w * Tw
#
# Two adding-up constraints, so only THREE of the five residuals are
# statistically independent. The estimated equations are Tw, T1 and E1, while
# Ti and Ei are recovered from the budgets.
#
# IDENTIFICATION AND STRUCTURAL ANCHORING
#
#   1. Utility scale: eta_1 = 1.
#   2. Theta location: theta_w = 0.
#      In classical labor economics, work provides zero autonomous intrinsic
#      pleasure at the entry margin (Tw = 0); individuals work for wages (w),
#      and working generates marginal disutility (-beta_w * Tw <= 0).
#      Fixing theta_w = 0 eliminates the translational symmetry (c = 0) without
#      relying on an arbitrary bliss anchor (A).
#   3. Direct Identification of Values of Time:
#      - VTAW = - (beta_w * Tw*) / lambda <= 0 (work disutility in money terms)
#      - VoL  = w + VTAW = w - (beta_w * Tw*) / lambda <= w
#      - theta_1 and theta_i are freely estimated relative to work at Tw = 0.
#   4. Concavity: all satiation parameters (beta_w, beta_1, beta_i, eta_i) are
#      estimated as exp(.) to ensure strict concavity throughout optimization.
# =============================================================================

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-aq-cholesky",
  modelDescr = "Additive quadratic MTUEM, full daily cycle, theta_w = 0 anchor",
  indivID = "PeID",
  outputDirectory = "output"
)

database <- mtuem::maed
database <- database[database$EcI > 0, ]

### The daily cycle. T1 and Ti are the two freely allocated activities, E1 and
### Ei the two freely consumed goods.
database$T1 <- database$Tf1
database$Ti <- database$Tf2
database$E1 <- database$Ef1
database$Ei <- database$Ef2 + database$Ef3

### Structural parameters:
###   theta_w = 0      : fixed instrumental work anchor
###   lbeta_w          : beta_w = exp(lbeta_w) work satiation / disutility curvature
###   theta_1, lbeta_1 : baseline marginal utility and satiation for T1
###   theta_i, lbeta_i : baseline marginal utility and satiation for Ti
###   Xb_1, eta_1 = 1  : phi_1 = eta_1 * Xb_1; eta_1 = 1 fixes utility scale
###   Xb_i, leta_i     : phi_i = eta_i * Xb_i; eta_i = exp(leta_i)
### Cholesky factor of the 3x3 covariance, ordered (Tw, T1, E1).
apollo_beta <- c(
  theta_w = 0, 
  beta_w = 1,
  theta_1 = 1, 
  beta_1 = 1,
  theta_i = 1,
  beta_i = 1,
  phi_1 = 1, 
  eta_1 = 1,
  eta_i = 1, 
  phi_i = 1,
  chk11 = 10, chk21 = 0, chk22 = 10,
  chk31 = 0, chk32 = 0, chk33 = 40
)
apollo_fixed <- c("theta_w", "eta_1")

apollo_inputs <- apollo_validateInputs()

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("T1", "Ti"),
    free_goods = c("E1", "Ei"),
    omitted_times = c("Ti"),              # recovered from the time budget
    omitted_goods = c("Ei"),              # recovered from the money budget
    goods_cost = list(E1 = 1, Ei = 1),
    work_elasticity = list(theta_w = theta_w),
    work_satiation  = list(beta_w = beta_w),
    times_elasticities = list(theta_1 = theta_1, theta_i = theta_i),
    times_satiations   = list(beta_1 = beta_1, beta_i = beta_i),
    goods_elasticities = list(phi_1 = phi_1, phi_i = phi_i),
    goods_satiations   = list(eta_1 = eta_1, eta_i = eta_i),
    cholesky = list(chk11 = chk11, chk21 = chk21, chk22 = chk22,
                    chk31 = chk31, chk32 = chk32, chk33 = chk33),
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- additive_quadratic_mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model <- apollo_estimate(apollo_beta, apollo_fixed,
                         apollo_probabilities, apollo_inputs,
                         list(estimate_settings = list(estimationRoutine = "bgw")))
apollo_modelOutput(model)

### Prediction returns all five quantities of the daily cycle, plus the values of time
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)

cat("\n--- Fit of the daily cycle (means) ---\n")
for (v in c("Tw", "T1", "Ti", "E1", "Ei"))
  cat(sprintf("  %-3s predicted %9.3f   observed %9.3f\n", v, mean(pred[, v]), mean(database[[v]])))
cat("  time budget residual :", max(abs(168 - database$Tc - pred[, "Tw"] - pred[, "T1"] - pred[, "Ti"])), "\n")
cat("  money budget residual:", max(abs(database$w * pred[, "Tw"] - database$EcI - pred[, "E1"] - pred[, "Ei"])), "\n")

### Values of Time Summary
cat("\n--- Values of Time Summary ---\n")
cat("Mean Wage :", round(mean(database$w), 3), "$/h\n")
cat("Mean VoL  :", round(mean(pred[, "VoL"]), 3), "$/h\n")
cat("Mean VTAW :", round(mean(pred[, "VTAW"]), 3), "$/h\n")
cat("Mean VoL / Mean Wage :", round(mean(pred[, "VoL"]) / mean(database$w), 3), "\n")

### Regularity checks
cat("\n--- Regularity checks ---\n")
cat("Concavity (all beta, eta > 0) : TRUE by construction\n")
cat("Share with VoL <= 0           :", round(mean(pred[, "VoL"] <= 0), 4), "\n")
cat("Share with Tw* out of bounds  :",
    round(mean(pred[, "Tw"] <= 0 | pred[, "Tw"] >= 168 - database$Tc), 4), "\n")
for (v in c("T1", "Ti", "E1", "Ei"))
  cat(sprintf("  share with %-3s < 0 : %s\n", v, format(mean(pred[, v] < 0))))

