# =============================================================================
# Comparative Analysis: Cobb-Douglas MTUEM vs Additive Quadratic MTUEM
# Dataset: MAED (Workers, EcI > 0)
# =============================================================================

rm(list = ls())
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))
set.seed(42)

database <- mtuem::maed
database <- database[database$EcI > 0, ]
database$T1 <- database$Tf1
database$Ti <- database$Tf2
database$E1 <- database$Ef1
database$Ei <- database$Ef2 + database$Ef3

# =============================================================================
# 1. ESTIMATE MODEL 1: COBB-DOUGLAS MTUEM (Jara-Diaz et al., 2008)
# =============================================================================
apollo_initialise()
apollo_control <- list(
  modelName  = "comp_cd_3eq",
  modelDescr = "Cobb-Douglas MTUEM (Jara-Diaz 2008)",
  indivID    = "PeID",
  outputDirectory = "output"
)

apollo_beta_cd <- c(
  Theta = 1,
  Phi = 1,
  thw = 0,
  theta_1 = 0.5,
  phi_1 = 0.5,
  chk11 = 100,
  chk21 = 0, chk22 = 100,
  chk31 = 0, chk32 = 0, chk33 = 100
)
apollo_fixed_cd <- c("Theta")
apollo_inputs_cd <- apollo_validateInputs(apollo_beta = apollo_beta_cd, apollo_fixed = apollo_fixed_cd)

apollo_probabilities_cd <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))

  P <- list()
  work_elasticities  <- list(Theta = Theta, Phi = Phi, thw = thw)
  times_elasticities <- list(theta_1 = theta_1)
  goods_elasticities <- list(phi_1 = phi_1)

  mtuem_settings <- list(
    work_times = c("Tw"),
    free_times = c("Tf1"),
    free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1),
    work_elasticities  = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    cholesky = list(chk11 = chk11, chk21 = chk21, chk22 = chk22,
                    chk31 = chk31, chk32 = chk32, chk33 = chk33),
    Tc = Tc, Ec = EcI, w = w, tau = 168
  )

  P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

cat("\n>>> Estimating Cobb-Douglas MTUEM...\n")
model_cd <- apollo_estimate(apollo_beta_cd, apollo_fixed_cd,
                            apollo_probabilities_cd, apollo_inputs_cd,
                            list(estimate_settings = list(estimationRoutine = "bgw")))

# =============================================================================
# 2. ESTIMATE MODEL 2: ADDITIVE QUADRATIC MTUEM (theta_w = 0 Anchor)
# =============================================================================
apollo_initialise()
apollo_control <- list(
  modelName  = "comp_aq_3eq",
  modelDescr = "Additive Quadratic MTUEM (theta_w = 0)",
  indivID    = "PeID",
  outputDirectory = "output"
)

apollo_beta_aq <- c(
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
apollo_fixed_aq <- c("theta_w", "eta_1")
apollo_inputs_aq <- apollo_validateInputs(apollo_beta = apollo_beta_aq, apollo_fixed = apollo_fixed_aq)

apollo_probabilities_aq <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
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

cat("\n>>> Estimating Additive Quadratic MTUEM...\n")
model_aq <- apollo_estimate(apollo_beta_aq, apollo_fixed_aq,
                            apollo_probabilities_aq, apollo_inputs_aq,
                            list(estimate_settings = list(estimationRoutine = "bgw")))

# =============================================================================
# 3. PREDICTIONS & DETAILED COMPARISONS
# =============================================================================
pred_cd <- apollo_prediction(model_cd, apollo_probabilities_cd, apollo_inputs_cd)
pred_aq <- apollo_prediction(model_aq, apollo_probabilities_aq, apollo_inputs_aq)

Tw_obs <- database$Tw
Tw_cd  <- pred_cd$Tw
Tw_aq  <- pred_aq[, "Tw"]

w  <- database$w
Tc <- database$Tc
Ec <- database$EcI

VoL_cd  <- pred_cd$VoL
VTAW_cd <- pred_cd$VTAW
VoL_aq  <- pred_aq[, "VoL"]
VTAW_aq <- pred_aq[, "VTAW"]

cat("\n====================================================================\n")
cat(" 1. ESTIMATION & STATISTICAL FIT COMPARISON\n")
cat("====================================================================\n")
fit_df <- data.frame(
  Metric = c("Log-Likelihood (LL)", "Free Parameters", "AIC", "BIC"),
  Cobb_Douglas = c(
    as.character(round(model_cd$maximum, 2)),
    as.character(length(model_cd$estimate) - length(apollo_fixed_cd)),
    as.character(round(2*(length(model_cd$estimate) - length(apollo_fixed_cd)) - 2*model_cd$maximum, 2)),
    as.character(round(-2*model_cd$maximum + (length(model_cd$estimate) - length(apollo_fixed_cd))*log(nrow(database)), 2))
  ),
  Additive_Quadratic = c(
    as.character(round(model_aq$maximum, 2)),
    as.character(length(model_aq$estimate) - length(apollo_fixed_aq)),
    as.character(round(2*(length(model_aq$estimate) - length(apollo_fixed_aq)) - 2*model_aq$maximum, 2)),
    as.character(round(-2*model_aq$maximum + (length(model_aq$estimate) - length(apollo_fixed_aq))*log(nrow(database)), 2))
  )
)
print(fit_df, row.names = FALSE)

cat("\n====================================================================\n")
cat(" 2. PREDICTED WORK TIME (Tw) SUMMARY & ACCURACY\n")
cat("====================================================================\n")
rmse_cd <- sqrt(mean((Tw_obs - Tw_cd)^2))
mae_cd  <- mean(abs(Tw_obs - Tw_cd))
cor_cd  <- cor(Tw_obs, Tw_cd)

rmse_aq <- sqrt(mean((Tw_obs - Tw_aq)^2))
mae_aq  <- mean(abs(Tw_obs - Tw_aq))
cor_aq  <- cor(Tw_obs, Tw_aq)

tw_summary <- data.frame(
  Quantity = c("Mean Tw (h/week)", "Std Dev Tw", "Min Tw", "Median Tw", "Max Tw",
               "RMSE (vs Observed)", "MAE (vs Observed)", "Correlation (with Observed)"),
  Observed = c(round(mean(Tw_obs), 3), round(sd(Tw_obs), 3), round(min(Tw_obs), 3),
               round(median(Tw_obs), 3), round(max(Tw_obs), 3), 0, 0, 1.000),
  Cobb_Douglas = c(round(mean(Tw_cd), 3), round(sd(Tw_cd), 3), round(min(Tw_cd), 3),
                   round(median(Tw_cd), 3), round(max(Tw_cd), 3),
                   round(rmse_cd, 3), round(mae_cd, 3), round(cor_cd, 3)),
  Additive_Quadratic = c(round(mean(Tw_aq), 3), round(sd(Tw_aq), 3), round(min(Tw_aq), 3),
                        round(median(Tw_aq), 3), round(max(Tw_aq), 3),
                        round(rmse_aq, 3), round(mae_aq, 3), round(cor_aq, 3))
)
print(tw_summary, row.names = FALSE)

cat("\n====================================================================\n")
cat(" 3. VALUES OF TIME COMPARISON\n")
cat("====================================================================\n")
vot_summary <- data.frame(
  Quantity = c("Mean Wage ($/h)", "Mean VoL ($/h)", "Mean VTAW ($/h)", "Mean VoL / Mean Wage", "Share VoL <= 0"),
  Cobb_Douglas = c(round(mean(w), 3), round(mean(VoL_cd), 3), round(mean(VTAW_cd), 3),
                   round(mean(VoL_cd) / mean(w), 3), round(mean(VoL_cd <= 0), 4)),
  Additive_Quadratic = c(round(mean(w), 3), round(mean(VoL_aq), 3), round(mean(VTAW_aq), 3),
                        round(mean(VoL_aq) / mean(w), 3), round(mean(VoL_aq <= 0), 4))
)
print(vot_summary, row.names = FALSE)

cat("\n====================================================================\n")
cat(" 4. DETAILED LABOUR SUPPLY (Tw*) EQUATION FORMULATION\n")
cat("====================================================================\n")
cat("A) COBB-DOUGLAS MTUEM (Jara-Diaz et al., 2008):\n")
cat("   - First-Order Condition:\n")
cat("     thw / Tw + w * Phi / (w*Tw - Ec) - Theta / (tau - Tw - Tc) = 0\n")
cat("   - Labour Supply Formula (Closed-form root of a quadratic):\n")
cat("     Tw* = 0.5 * (tau - Tc) * [ thetaphiec + sqrt( thetaphiec^2 - 4*thw*ec*(Theta + Phi + thw) ) ] / (Theta + Phi + thw)\n")
cat("     where ec = Ec / [ w * (tau - Tc) ] and thetaphiec = (Phi + thw) + (Theta + thw)*ec\n")
cat("   - Key Property: Non-linear, captures proportional expenditure rate allocation.\n\n")

cat("B) ADDITIVE QUADRATIC MTUEM:\n")
cat("   - First-Order Condition:\n")
cat("     theta_w - beta_w*Tw + w * [ (Sph - (w*Tw - Ec))/H ] - [ (Sth - (tau - Tc - Tw))/B ] = 0\n")
cat("   - Labour Supply Formula (Linear-rational):\n")
cat("     Tw* = [ theta_w + w * (Sph + Ec)/H - (Sth - (tau - Tc))/B ] / [ beta_w + w^2/H + 1/B ]\n")
cat("     where Sth = Sum(theta_i/beta_i), B = Sum(1/beta_i), Sph = Sum(P_j*phi_j/eta_j), H = Sum(P_j^2/eta_j)\n")
cat("   - Key Property: Linear in committed time and committed expenses, balance of marginal utilities and bliss points.\n")
cat("====================================================================\n")
