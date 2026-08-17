# =============================================================================
# Script to Generate and Save Publication-Quality Plots:
# 1. Predictions against Real Working Time (Observed vs Predicted)
# 2. Prediction Cobb-Douglas vs Prediction Additive Quadratic
# 3. Labour Supply Curves (w on Y-axis, Tw on X-axis) for representative worker
# Output Folder: plots/
# =============================================================================

rm(list = ls())
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))
library(ggplot2)
set.seed(42)

# Load and prepare database
database <- mtuem::maed
database <- database[database$EcI > 0, ]
database$T1 <- database$Tf1
database$Ti <- database$Tf2
database$E1 <- database$Ef1
database$Ei <- database$Ef2 + database$Ef3

# Ensure output directories exist
plots_dir <- file.path(getwd(), "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

art_dir <- "C:/Users/pablo/.gemini/antigravity-cli/brain/a86158b8-36d4-4963-b612-60ac0d54c83e"
if (!dir.exists(art_dir)) dir.create(art_dir, recursive = TRUE)

# =============================================================================
# 1. ESTIMATE COBB-DOUGLAS MTUEM (Jara-Diaz et al., 2008)
# =============================================================================
apollo_initialise()
apollo_control <- list(
  modelName  = "plot_cd_3eq",
  modelDescr = "Cobb-Douglas MTUEM",
  indivID    = "PeID",
  outputDirectory = "output"
)

apollo_beta_cd <- c(
  Theta = 1, Phi = 1, thw = 0, theta_1 = 0.5, phi_1 = 0.5,
  chk11 = 100, chk21 = 0, chk22 = 100, chk31 = 0, chk32 = 0, chk33 = 100
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
    work_times = c("Tw"), free_times = c("Tf1"), free_goods = c("Ef1"),
    goods_cost = list(Ef1 = 1), work_elasticities = work_elasticities,
    times_elasticities = times_elasticities, goods_elasticities = goods_elasticities,
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
# 2. ESTIMATE ADDITIVE QUADRATIC MTUEM (theta_w = 0, Linear units)
# =============================================================================
apollo_initialise()
apollo_control <- list(
  modelName  = "plot_aq_3eq",
  modelDescr = "Additive Quadratic MTUEM",
  indivID    = "PeID",
  outputDirectory = "output"
)

apollo_beta_aq <- c(
  theta_w = 0, beta_w = 1, theta_1 = 1, beta_1 = 1, theta_i = 1, beta_i = 1,
  phi_1 = 1, eta_1 = 1, eta_i = 1, phi_i = 1,
  chk11 = 10, chk21 = 0, chk22 = 10, chk31 = 0, chk32 = 0, chk33 = 40
)
apollo_fixed_aq <- c("theta_w", "eta_1")
apollo_inputs_aq <- apollo_validateInputs(apollo_beta = apollo_beta_aq, apollo_fixed = apollo_fixed_aq)

apollo_probabilities_aq <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  mtuem_settings <- list(
    work_times = c("Tw"), free_times = c("T1", "Ti"), free_goods = c("E1", "Ei"),
    omitted_times = c("Ti"), omitted_goods = c("Ei"), goods_cost = list(E1 = 1, Ei = 1),
    work_elasticity = list(theta_w = theta_w), work_satiation = list(beta_w = beta_w),
    times_elasticities = list(theta_1 = theta_1, theta_i = theta_i),
    times_satiations = list(beta_1 = beta_1, beta_i = beta_i),
    goods_elasticities = list(phi_1 = phi_1, phi_i = phi_i),
    goods_satiations = list(eta_1 = eta_1, eta_i = eta_i),
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
# 3. COMPUTE PREDICTIONS
# =============================================================================
pred_cd <- apollo_prediction(model_cd, apollo_probabilities_cd, apollo_inputs_cd)
pred_aq <- apollo_prediction(model_aq, apollo_probabilities_aq, apollo_inputs_aq)

df_plot <- data.frame(
  Tw_obs = database$Tw,
  Tw_cd  = pred_cd$Tw,
  Tw_aq  = pred_aq[, "Tw"],
  w      = database$w,
  Tc     = database$Tc,
  Ec     = database$EcI
)

# =============================================================================
# 4. PLOT 1: PREDICTIONS AGAINST REAL WORKING TIME (Observed vs Predicted)
# =============================================================================
p1_cd <- ggplot(df_plot, aes(x = Tw_obs, y = Tw_cd)) +
  geom_point(color = "#2b5c8f", alpha = 0.45, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#d95f02", linewidth = 0.9) +
  geom_smooth(method = "lm", color = "#2b5c8f", se = FALSE, linewidth = 0.9) +
  labs(
    title = "Cobb-Douglas MTUEM",
    subtitle = sprintf("RMSE: %.2f h | MAE: %.2f h | r = %.3f",
                       sqrt(mean((df_plot$Tw_obs - df_plot$Tw_cd)^2)),
                       mean(abs(df_plot$Tw_obs - df_plot$Tw_cd)),
                       cor(df_plot$Tw_obs, df_plot$Tw_cd)),
    x = "Observed Work Time (hours/week)",
    y = "Predicted Work Time (hours/week)"
  ) +
  xlim(0, 80) + ylim(0, 80) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 13, color = "#111111"),
    plot.subtitle = element_text(size = 10, color = "#444444")
  )

p1_aq <- ggplot(df_plot, aes(x = Tw_obs, y = Tw_aq)) +
  geom_point(color = "#1b9e77", alpha = 0.45, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#d95f02", linewidth = 0.9) +
  geom_smooth(method = "lm", color = "#1b9e77", se = FALSE, linewidth = 0.9) +
  labs(
    title = "Additive Quadratic MTUEM",
    subtitle = sprintf("RMSE: %.2f h | MAE: %.2f h | r = %.3f",
                       sqrt(mean((df_plot$Tw_obs - df_plot$Tw_aq)^2)),
                       mean(abs(df_plot$Tw_obs - df_plot$Tw_aq)),
                       cor(df_plot$Tw_obs, df_plot$Tw_aq)),
    x = "Observed Work Time (hours/week)",
    y = "Predicted Work Time (hours/week)"
  ) +
  xlim(0, 80) + ylim(0, 80) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 13, color = "#111111"),
    plot.subtitle = element_text(size = 10, color = "#444444")
  )

p1_combined <- gridExtra::arrangeGrob(
  p1_cd, p1_aq, ncol = 2,
  top = grid::textGrob("Predicted vs Observed Working Time (Tw)",
                       gp = grid::gpar(fontsize = 15, fontface = "bold"))
)

cat("Saving Plot 1: Predictions vs Observed (300 DPI)...\n")
ggsave(file.path(plots_dir, "plot1_pred_vs_obs.png"), plot = p1_combined, width = 10, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(art_dir, "plot1_pred_vs_obs.png"), plot = p1_combined, width = 10, height = 5.2, dpi = 300, bg = "white")

# =============================================================================
# 5. PLOT 2: PREDICTION COBB-DOUGLAS VS PREDICTION ADDITIVE QUADRATIC
# =============================================================================
p2 <- ggplot(df_plot, aes(x = Tw_cd, y = Tw_aq)) +
  geom_point(color = "#7570b3", alpha = 0.55, size = 2.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#e7298a", linewidth = 0.9) +
  geom_smooth(method = "lm", color = "#386cb0", se = TRUE, alpha = 0.2, linewidth = 0.9) +
  labs(
    title = "Cobb-Douglas vs Additive Quadratic Predictions",
    subtitle = sprintf("Correlation r = %.3f | Mean CD: %.2f h | Mean AQ: %.2f h",
                       cor(df_plot$Tw_cd, df_plot$Tw_aq), mean(df_plot$Tw_cd), mean(df_plot$Tw_aq)),
    x = "Cobb-Douglas Predicted Work Time Tw (hours/week)",
    y = "Additive Quadratic Predicted Work Time Tw (hours/week)"
  ) +
  xlim(10, 75) + ylim(10, 75) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 13.5, color = "#111111"),
    plot.subtitle = element_text(size = 10.5, color = "#444444")
  )

cat("Saving Plot 2: Cobb-Douglas vs Additive Quadratic Predictions (300 DPI)...\n")
ggsave(file.path(plots_dir, "plot2_cd_vs_aq.png"), plot = p2, width = 7.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(art_dir, "plot2_cd_vs_aq.png"), plot = p2, width = 7.5, height = 5.5, dpi = 300, bg = "white")

# =============================================================================
# 6. PLOT 3: LABOUR SUPPLY CURVES FOR REPRESENTATIVE INDIVIDUAL
# =============================================================================
mean_Tc <- mean(database$Tc)
mean_Ec <- mean(database$EcI)
tau     <- 168

w_seq <- seq(4.5, 35, length.out = 300)

# Cobb-Douglas simulation
b_cd <- model_cd$estimate
Theta_cd <- 1
Phi_cd   <- b_cd["Phi"]
thw_cd   <- b_cd["thw"]

tw_sim_cd <- sapply(w_seq, function(w_val) {
  ec_val <- mean_Ec / (w_val * (tau - mean_Tc))
  base   <- Theta_cd + Phi_cd + thw_cd
  thetaphiec <- (Phi_cd + thw_cd) + (Theta_cd + thw_cd) * ec_val
  discriminant <- thetaphiec^2 - 4 * thw_cd * ec_val * base
  if (discriminant < 0) return(NA)
  tw <- 0.5 * (tau - mean_Tc) * (thetaphiec + sqrt(discriminant)) / base
  return(tw)
})

# Additive Quadratic simulation
b_aq <- model_aq$estimate
theta_w_aq <- 0
beta_w_aq  <- b_aq["beta_w"]
theta_1_aq <- b_aq["theta_1"]
beta_1_aq  <- b_aq["beta_1"]
theta_i_aq <- b_aq["theta_i"]
beta_i_aq  <- b_aq["beta_i"]
phi_1_aq   <- b_aq["phi_1"]
eta_1_aq   <- 1
phi_i_aq   <- b_aq["phi_i"]
eta_i_aq   <- b_aq["eta_i"]

Sth_aq <- theta_1_aq / beta_1_aq + theta_i_aq / beta_i_aq
B_aq   <- 1 / beta_1_aq + 1 / beta_i_aq
Sph_aq <- phi_1_aq / eta_1_aq + phi_i_aq / eta_i_aq
H_aq   <- 1 / eta_1_aq + 1 / eta_i_aq

tw_sim_aq <- sapply(w_seq, function(w_val) {
  num_1 <- theta_w_aq
  num_2 <- w_val * (Sph_aq + mean_Ec) / H_aq
  num_3 <- (Sth_aq - (tau - mean_Tc)) / B_aq
  den   <- beta_w_aq + w_val^2 / H_aq + 1 / B_aq
  tw    <- (num_1 + num_2 - num_3) / den
  return(tw)
})

df_supply <- rbind(
  data.frame(w = w_seq, Tw = tw_sim_cd, Model = "Cobb-Douglas MTUEM (Jara-Diaz 2008)"),
  data.frame(w = w_seq, Tw = tw_sim_aq, Model = "Additive Quadratic MTUEM (theta_w = 0)")
)
df_supply <- df_supply[!is.na(df_supply$Tw), ]

p3 <- ggplot() +
  # Real observed data points
  geom_point(data = database, aes(x = Tw, y = w),
             color = "#8c96c6", alpha = 0.35, size = 2.2, shape = 16) +
  # Simulated curves
  geom_path(data = df_supply, aes(x = Tw, y = w, color = Model, linetype = Model), linewidth = 1.3) +
  # Sample mean point
  geom_point(data = data.frame(Tw = mean(database$Tw), w = mean(database$w)),
             aes(x = Tw, y = w), inherit.aes = FALSE,
             color = "#d95f02", size = 4.2, shape = 18) +
  annotate("text", x = mean(database$Tw) + 1.2, y = mean(database$w) + 1.4,
           label = sprintf("Sample Mean\n(Tw = %.1f h, w = $%.2f)", mean(database$Tw), mean(database$w)),
           hjust = 0, size = 3.6, fontface = "bold", color = "#111111") +
  scale_color_manual(values = c("Cobb-Douglas MTUEM (Jara-Diaz 2008)" = "#08519c",
                                "Additive Quadratic MTUEM (theta_w = 0)" = "#006d2c")) +
  scale_linetype_manual(values = c("Cobb-Douglas MTUEM (Jara-Diaz 2008)" = "solid",
                                  "Additive Quadratic MTUEM (theta_w = 0)" = "dashed")) +
  xlim(0, 80) + ylim(0, 45) +
  labs(
    title = "Labour Supply Curves for the Representative Individual",
    subtitle = sprintf("Simulated for representative worker (Tc = %.2f h, Ec = $%.2f) overlaid on observed sample data (N = %d)",
                       mean_Tc, mean_Ec, nrow(database)),
    x = "Work Time Tw (hours/week)",
    y = "Hourly Wage Rate w ($/hour)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey75"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 13.5, color = "#111111"),
    plot.subtitle = element_text(size = 10, color = "#444444"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

cat("Saving Plot 3: Labour Supply Curves (300 DPI)...\n")
ggsave(file.path(plots_dir, "plot3_labour_supply.png"), plot = p3, width = 8, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(art_dir, "plot3_labour_supply.png"), plot = p3, width = 8, height = 5.8, dpi = 300, bg = "white")

cat("\n>>> All plots successfully updated with high DPI, white backgrounds, and real data scatter!\n")
