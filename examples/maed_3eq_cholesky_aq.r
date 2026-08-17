rm(list = ls())
library(apollo)
library(ggplot2)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))
set.seed(42)

# =============================================================================
# Additive quadratic MTUEM: full daily cycle, three equations, full covariance
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
# IDENTIFICATION (see "Quadratic additive.md", section "Identification rules")
#
# The model has three exact invariances. Adding a multiple of a binding budget
# to the utility function changes it by a constant on the feasible set, so it
# cannot move the optimum, but it moves that budget's multiplier one for one.
#
#   I1 scale     : multiply every parameter by a > 0. Harmless, VoL = mu/lambda
#                  is invariant.
#   I2 theta loc.: add c to theta_w and to every theta_i. Moves mu to mu + c and
#                  leaves lambda alone, so it moves VoL and VTAW. Not harmless.
#   I3 phi loc.  : add P_j*d to every phi_j and -w*d to theta_w. Moves lambda.
#                  Only admissible if theta_w carries the wage, so it is blocked
#                  by an exclusion restriction rather than by a normalization.
#
#   R1 scale             : eta_1 = 1.
#   R2 theta location    : theta_w = beta_w * A_w, i.e. anchor the bliss point of
#                          work at A_w hours. Estimated here at A_w = 0, which
#                          fixes theta_w = 0; any other anchor is the same fit
#                          reached by I2 with c = beta_w * A_w, so it is applied
#                          after estimation instead of re-estimating.
#   R3 exclusions        : the wage must not enter theta_w and P_j must not enter
#                          phi_j, otherwise I3 becomes live and the level of
#                          lambda needs an assumption too. Other covariates are
#                          fine.
#   R4 complete free sets: S_theta, B, S_phi and H are literal sums over the
#                          categories supplied, so a category left out is deleted
#                          from the model, not absorbed into an aggregate. Ti and
#                          Ei stay in the free sets and only their equations are
#                          dropped, through omitted_times / omitted_goods.
#
# Identified: all beta and eta, all contrasts theta_k - theta_l, all phi, the
# allocations, and VoL - VTAW = w. NOT identified: the level of mu, hence the
# levels of VoL and VTAW, which are conditional on A_w and must be reported with
# it. The log-likelihood is flat in A_w, so it cannot be estimated.
#
# VALUES OF TIME. With theta_w = beta_w * A_w and lambda > 0,
#
#     VTAW = beta_w (A_w - Tw*) / lambda        VTAW < 0  <=>  Tw* > A_w
#     VoL  = w + VTAW                           mu   > 0  <=>  VTAW > -w
#
# A_w is the bliss point of work: the hours a respondent would choose if work
# paid nothing. Work is a marginal bad above it and a marginal good below it,
# which is the quadratic counterpart of the sign of theta_w in the Cobb-Douglas
# models, and the data decide who is on which side through Tw*. Raising A_w buys
# time scarcity (mu > 0) and spends the negative value of work, monotonically:
# the sweep at the end of this script prints the whole trade-off.
#
# CONCAVITY. beta_i > 0 and eta_j > 0 are required. For work the exact condition
# is weaker, beta_w + 1/B + w^2/H > 0, the denominator of Tw*. The satiations are
# estimated freely here and checked after estimation; parameterize them as exp(.)
# if a run on other data drifts negative. Do NOT impose Tw* < theta_w/beta_w and
# do NOT estimate theta_w as exp(.): both force VTAW > 0.
# =============================================================================

apollo_initialise()
apollo_control <- list(
  modelName = "maed-3eq-aq-cholesky",
  modelDescr = "Additive quadratic MTUEM, full daily cycle, work bliss anchor",
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

tau <- 168

### Reporting anchor, R2. A_w = 0 reads as "work has no intrinsic reward even at
### the first hour". A_w = 33.3 reproduces the VoL level of the Cobb-Douglas
### model on this sample. Changing it does not change the fit, only the level of
### the values of time and the share of the sample with VTAW < 0.
A_w <- 0

### Structural parameters:
###   theta_w = 0      : work bliss anchor at A_w = 0 (R2)
###   beta_w           : work satiation, the curvature of the disutility of work
###   theta_1, beta_1  : baseline marginal utility and satiation for T1
###   theta_i, beta_i  : baseline marginal utility and satiation for Ti
###   phi_1, eta_1 = 1 : eta_1 = 1 fixes the utility scale (R1)
###   phi_i, eta_i     : baseline marginal utility and satiation for Ei
### Cholesky factor of the 3x3 covariance, ordered (Tw, T1, E1).
apollo_beta <- c(
  theta_w = 0,
  beta_w = 1,
  theta_1 = 0,
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

### Prediction returns all five quantities of the daily cycle
pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
head(pred)

cat("\n--- Fit of the daily cycle (means) ---\n")
for (v in c("Tw", "T1", "Ti", "E1", "Ei"))
  cat(sprintf("  %-3s predicted %9.3f   observed %9.3f\n", v, mean(pred[, v]), mean(database[[v]])))
cat("  time budget residual :", max(abs(tau - database$Tc - pred[, "Tw"] - pred[, "T1"] - pred[, "Ti"])), "\n")
cat("  money budget residual:", max(abs(database$w * pred[, "Tw"] - database$EcI - pred[, "E1"] - pred[, "Ei"])), "\n")

### -------------------------------------------------------------------------
### Shadow prices at the chosen anchor
### -------------------------------------------------------------------------
### Moving the anchor from 0 to A_w is invariance I2 with c = beta_w * A_w: every
### theta rises by c, mu rises by c, lambda and the allocations do not move.
b   <- model$estimate
ess <- get_additive_quadratic_essentials(
  list(b["theta_1"], b["theta_i"]), list(b["beta_1"], b["beta_i"]),
  list(b["phi_1"],   b["phi_i"]),   list(b["eta_1"],  b["eta_i"]), list(1, 1))

Tw_opt <- pred[, "Tw"]
lambda <- as.numeric((ess$Sph - (database$w * Tw_opt - database$EcI)) / ess$H)
mu0    <- as.numeric((ess$Sth - (tau - database$Tc - Tw_opt)) / ess$B)   # at A_w = 0

shift  <- as.numeric(b["beta_w"] * A_w)
mu     <- mu0 + shift
VoL    <- mu / lambda
VTAW   <- VoL - database$w

cat("\n--- Structural parameters at the reporting anchor A_w =", A_w, "h ---\n")
cat(sprintf("  theta_w %10.3f   beta_w %8.3f\n", b["beta_w"] * A_w, b["beta_w"]))
cat(sprintf("  theta_1 %10.3f   beta_1 %8.3f   bliss B_1 %7.3f h\n",
            b["theta_1"] + shift, b["beta_1"], (b["theta_1"] + shift) / b["beta_1"]))
cat(sprintf("  theta_i %10.3f   beta_i %8.3f   bliss B_i %7.3f h\n",
            b["theta_i"] + shift, b["beta_i"], (b["theta_i"] + shift) / b["beta_i"]))
cat(sprintf("  identified contrasts: theta_w - theta_1 %10.3f   theta_w - theta_i %10.3f\n",
            b["beta_w"] * A_w - (b["theta_1"] + shift), b["beta_w"] * A_w - (b["theta_i"] + shift)))

cat("\n--- Values of time at A_w =", A_w, "h ---\n")
cat(sprintf("  mean wage %8.3f $/h\n", mean(database$w)))
cat(sprintf("  mean VoL  %8.3f $/h   (%.3f x mean wage)\n", mean(VoL), mean(VoL) / mean(database$w)))
cat(sprintf("  mean VTAW %8.3f $/h   share VTAW < 0 : %.3f\n", mean(VTAW), mean(VTAW < 0)))
cat("  VoL - VTAW - w (invariant, should be 0):", max(abs(VoL - VTAW - database$w)), "\n")

### -------------------------------------------------------------------------
### Regularity checks
### -------------------------------------------------------------------------
cat("\n--- Regularity checks ---\n")
cat("  concavity, all beta and eta > 0 :",
    all(c(b["beta_w"], b["beta_1"], b["beta_i"], b["eta_1"], b["eta_i"]) > 0), "\n")
cat("  work curvature, beta_w + 1/B + w^2/H > 0 :",
    all(b["beta_w"] + 1 / ess$B + database$w^2 / ess$H > 0), "\n")
cat(sprintf("  share mu > 0 (time is scarce, set by A_w) : %.4f\n", mean(mu > 0)))
cat(sprintf("  share lambda > 0 (goods below bliss)      : %.4f\n", mean(lambda > 0)))
cat(sprintf("  share Tw* outside [0, tau - Tc]           : %.4f\n",
            mean(Tw_opt <= 0 | Tw_opt >= tau - database$Tc)))
for (v in c("T1", "Ti", "E1", "Ei"))
  cat(sprintf("  share %-3s < 0 : %.4f\n", v, mean(pred[, v] < 0)))

### -------------------------------------------------------------------------
### Anchor sensitivity: the level of the values of time is a choice, not a result
### -------------------------------------------------------------------------
### No re-estimation is involved: mu(A_w) = mu(0) + beta_w * A_w and nothing else
### moves, so the log-likelihood and every allocation are identical down the table.
cat("\n--- Anchor sensitivity (same fit, same allocations) ---\n")
cat("     A_w   share mu>0   share VTAW<0   mean VoL   VoL/mean w   mean VTAW\n")
for (A in c(0, 10, 20, 33.3, 41.1, 58.1)) {
  muA   <- mu0 + as.numeric(b["beta_w"]) * A
  VoLA  <- muA / lambda
  VTAWA <- VoLA - database$w
  cat(sprintf("  %6.1f   %10.3f   %12.3f   %8.3f   %10.3f   %9.3f\n",
              A, mean(muA > 0), mean(VTAWA < 0), mean(VoLA),
              mean(VoLA) / mean(database$w), mean(VTAWA)))
}

### -------------------------------------------------------------------------
### Labour supply curve for a representative worker
### -------------------------------------------------------------------------
w_seq  <- seq(4.5, 35, length.out = 300)
tw_sim <- as.numeric(get_tw_additive_quadratic(
  list(b["beta_w"] * A_w), list(b["beta_w"]), ess$Sth, ess$Sph, ess$B, ess$H,
  tau, mean(database$Tc), mean(database$EcI), w_seq))
df_supply <- data.frame(w = w_seq, Tw = tw_sim)

p_supply <- ggplot(database, aes(x = Tw, y = w)) +
  geom_point(colour = "grey65", alpha = 0.4, size = 1.6) +
  geom_path(data = df_supply, colour = "#006d2c", linewidth = 1.1) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 45)) +
  labs(title = "Labour supply, additive quadratic MTUEM",
       subtitle = sprintf("representative worker at mean Tc and Ec, N = %d", nrow(database)),
       x = "Work time (hours/week)", y = "Wage rate ($/hour)") +
  theme_bw(base_size = 11)

if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
ggsave(file.path("plots", "labour_supply_aq.png"), p_supply,
       width = 7, height = 4.5, dpi = 200, bg = "white")
print(p_supply)
