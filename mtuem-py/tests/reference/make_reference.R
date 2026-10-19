suppressMessages({library(apollo); library(mtuem); library(jsonlite); library(dplyr)})
args <- commandArgs(TRUE); case <- args[1]; outdir <- args[2]; setwd(outdir)
set.seed(42)
apollo_initialise()
lc <- grepl("^lc", case)
apollo_control <- list(modelName = case, modelDescr = case, indivID = "PeID", outputDirectory = ".")
if (lc) { apollo_control$noValidation <- TRUE; apollo_control$noDiagnostics <- TRUE }
est_set <- list(writeIter = FALSE, maxIterations = 500, estimationRoutine = "bgw", hessianRoutine = "numDeriv")
if (lc) { est_set$scaleHessian <- FALSE; est_set$scaleAfterConvergence <- FALSE; est_set$silent <- FALSE }
if (case == "lc_cholesky") est_set$estimationRoutine <- NULL
eq1like <- case %in% c("eq1", "eq1_bounded")
if (case == "eq1_bounded") {
  # A lower bound on thw that binds: the unconstrained eq1 estimate plus 10% of its absolute value,
  # imposed with maxLik's adaptive barrier (constrOptim) through Apollo's BFGS routine.
  thw_free <- fromJSON("eq1_ref.json")$estimate$thw
  thw_lower <- as.numeric(sprintf("%.15g", thw_free + 0.1 * abs(thw_free)))
  est_set$estimationRoutine <- "bfgs"
  est_set$constraints <- sprintf("thw > %.15g", thw_lower)
}
database <- mtuem::maed; database <- database[database$EcI > 0, ]

if (case == "eq1") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0)
} else if (case == "eq1_bounded") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0.5)
} else if (case == "eq3_profiled") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0, theta_1 = 0.5, phi_1 = 0.5)
} else if (case == "eq3_sig") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0, theta_1 = 0.5, phi_1 = 0.5, sig_1 = 100, sig_2 = 100, sig_3 = 100)
} else if (case == "eq3_full") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0, theta_1 = 0.5, phi_1 = 0.5, sig_1 = 100, sig_2 = 100, sig_3 = 100, rho_12 = 0, rho_13 = 0, rho_23 = 0)
} else if (case == "eq3_cholesky") {
  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0, theta_1 = 0.5, phi_1 = 0.5, chk11 = 100, chk21 = 0, chk22 = 100, chk31 = 0, chk32 = 0, chk33 = 100)
}
apollo_fixed <- "Theta"
if (!lc) {
  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    s <- list(work_times = "Tw", work_elasticities = list(Theta = Theta, Phi = Phi, thw = thw), Tc = Tc, Ec = EcI, w = w, tau = 168)
    if (!eq1like) {
      s$free_times <- "Tf1"; s$free_goods <- "Ef1"; s$goods_cost <- list(Ef1 = 1)
      s$times_elasticities <- list(theta_1 = theta_1); s$goods_elasticities <- list(phi_1 = phi_1)
    }
    if (case %in% c("eq3_sig", "eq3_full")) s$sig <- list(sig_1 = sig_1, sig_2 = sig_2, sig_3 = sig_3)
    if (case == "eq3_full") s$rho <- list(rho_12 = rho_12, rho_13 = rho_13, rho_23 = rho_23)
    if (case == "eq3_cholesky") s$cholesky <- list(chk11 = chk11, chk21 = chk21, chk22 = chk22, chk31 = chk31, chk32 = chk32, chk33 = chk33)
    P[["model"]] <- mtuem_likelihood(s, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
} else {
  chk_n <- c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33")
  base <- c(Theta = 1, Phi = 1, thw = 0, th1 = 0.3, ph1 = 0.3)
  memb <- c(asc = 0, x_female = 0, x_older45 = 0)
  chk_start <- c(chk11 = 10, chk21 = 0, chk22 = 20, chk31 = 0, chk32 = 0, chk33 = 20)
  mk <- function(v, s) setNames(v, paste0(names(v), "_", s))
  if (case == "lc_cholesky") {
    apollo_beta <- c(mk(base, 1), mk(chk_start, 1), mk(memb, 1), mk(base, 2), mk(chk_start, 2), mk(memb, 2))
  } else if (case == "lc_profiled") {
    apollo_beta <- c(mk(base, 1), mk(memb, 1), mk(base, 2), mk(memb, 2))
  } else {
    apollo_beta <- c(mk(base, 1), mk(memb, 1), mk(base, 2), mk(memb, 2), c(chk11 = 6, chk21 = 0, chk22 = 7, chk31 = 0, chk32 = 0, chk33 = 30))
  }
  apollo_fixed <- c("Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1")
  apollo_lcPars <- function(apollo_beta, apollo_inputs) {
    lcpars <- list()
    lcpars[["Theta"]] <- list(Theta_1, Theta_2)
    lcpars[["Phi"]] <- list(Phi_1, Phi_2)
    lcpars[["thw"]] <- list(thw_1, thw_2)
    lcpars[["th1"]] <- list(th1_1, th1_2)
    lcpars[["ph1"]] <- list(ph1_1, ph1_2)
    if (case == "lc_cholesky") {
      lcpars[["chk11"]] <- list(chk11_1, chk11_2)
      lcpars[["chk21"]] <- list(chk21_1, chk21_2)
      lcpars[["chk22"]] <- list(chk22_1, chk22_2)
      lcpars[["chk31"]] <- list(chk31_1, chk31_2)
      lcpars[["chk32"]] <- list(chk32_1, chk32_2)
      lcpars[["chk33"]] <- list(chk33_1, chk33_2)
    }
    V <- list()
    V[["class_1"]] <- asc_1 + female * x_female_1 + older45 * x_older45_1
    V[["class_2"]] <- asc_2 + female * x_female_2 + older45 * x_older45_2
    lcpars[["pi_values"]] <- apollo_classAlloc(list(classes = c(class_1 = 1, class_2 = 2), utilities = V))
    return(lcpars)
  }
  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    cs <- list()
    for (s in 1:length(pi_values)) {
      st <- list(work_times = "Tw", free_times = "Tf1", free_goods = "Ef1", goods_cost = list(Ef1 = 1),
                 work_elasticities = list(Theta = Theta[[s]], Phi = Phi[[s]], thw = thw[[s]]),
                 times_elasticities = list(th1 = th1[[s]]), goods_elasticities = list(ph1 = ph1[[s]]),
                 Tc = Tc, Ec = EcI, w = w, tau = 168, componentName = paste0("class_", s))
      if (case == "lc_cholesky") st$cholesky <- list(chk11[[s]], chk21[[s]], chk22[[s]], chk31[[s]], chk32[[s]], chk33[[s]])
      if (case == "lc_shared") st$cholesky <- list(chk11, chk21, chk22, chk31, chk32, chk33)
      if (case == "lc_profiled") {
        cs[[s]] <- st
      } else {
        P[[paste0("class_", s)]] <- mtuem_likelihood(st, functionality)
      }
    }
    if (case == "lc_profiled") {
      P[["model"]] <- mtuem_lc_profiled_likelihood(list(class_settings = cs, classProb = pi_values), functionality)
    } else {
      lc_settings <- list(inClassProb = P, classProb = pi_values)
      P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
    }
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
}
apollo_inputs <- apollo_validateInputs()

llf <- function(b) { p <- apollo_probabilities(b, apollo_inputs, "estimate"); if (is.list(p)) p <- p[["model"]]; sum(log(p)) }
b_test <- apollo_beta
free <- setdiff(names(b_test), apollo_fixed)
b_test[free] <- b_test[free] + seq(0.05, by = 0.01, length.out = length(free))
out <- list(case = case, beta_start = as.list(apollo_beta), fixed = apollo_fixed,
            ll_start = llf(apollo_beta), beta_test = as.list(b_test), ll_test = llf(b_test))

t0 <- Sys.time()
if (case == "lc_cholesky") {
  model <- apollo_lcEM(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                       lcEM_settings = list(EMmaxIterations = 10), estimate_settings = est_set)
} else {
  model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings = est_set)
}
if (case == "eq1_bounded") {
  # maxLik's barrier (Apollo's bfgs with constraints) stops short here: constrOptim2 leaves after one
  # outer iteration when maximising, and optim's BFGS stalls against the barrier before Phi has
  # converged (restarting from its estimate does not move Phi). Its result is kept for reference.
  # The bound binds, so the constrained maximum is the maximum with thw fixed at the bound (bgw);
  # the covariance matrix is then that of the unpenalised log-likelihood over all free parameters
  # at that point (numDeriv Hessian and BHHH sandwich), as maxLik and Apollo report it after constrOptim.
  out$barrier <- list(estimate = as.list(model$estimate), se = as.list(model$se), robse = as.list(model$robse),
                      maximum = model$maximum, ll_at_estimate = llf(model$estimate))
  est_kkt <- est_set; est_kkt$estimationRoutine <- "bgw"; est_kkt$constraints <- NULL
  fixed_free <- apollo_fixed
  apollo_beta["thw"] <- thw_lower; apollo_fixed <- c(fixed_free, "thw")
  apollo_inputs <- apollo_validateInputs()
  kkt <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings = est_kkt)
  model <- kkt
  free <- setdiff(names(model$estimate), fixed_free)
  lli <- function(bf) {
    b <- model$estimate; b[free] <- bf
    p <- apollo_probabilities(b, apollo_inputs, "estimate"); if (is.list(p)) p <- p[["model"]]; log(p)
  }
  H <- numDeriv::hessian(function(bf) sum(lli(bf)), model$estimate[free])
  S <- numDeriv::jacobian(lli, model$estimate[free])
  vc <- solve(-H); rvc <- vc %*% crossprod(S) %*% vc
  dimnames(H) <- dimnames(vc) <- dimnames(rvc) <- list(free, free)
}
out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
apollo_modelOutput(model)
apollo_saveOutput(model, saveOutput_settings = list(printCovar = TRUE, printCorr = TRUE))
if (case == "eq1_bounded") {
  apollo_beta <- model$estimate; apollo_fixed <- fixed_free
  apollo_inputs <- apollo_validateInputs()
  model$apollo_fixed <- fixed_free; model$hessian <- H; model$eigValue <- max(eigen(H)$values)
  model$varcov <- vc; model$robvarcov <- rvc
  model$se <- model$robse <- setNames(rep(NA_real_, length(model$estimate)), names(model$estimate))
  model$se[free] <- sqrt(diag(vc)); model$robse[free] <- sqrt(diag(rvc))
}
out$estimate <- as.list(model$estimate); out$se <- as.list(model$se); out$robse <- as.list(model$robse)
out$maximum <- model$maximum; out$eigValue <- model$eigValue; out$code <- model$code
out$ll_at_estimate <- llf(model$estimate)
if (case == "eq1_bounded") out$bounds <- list(thw = list(lower = thw_lower))

if (lc) {
  vot <- mtuem_values_of_time(model, apollo_inputs, apollo_probabilities, nClass = 2, Ec = "EcI")
} else {
  vot <- mtuem_values_of_time(model, apollo_inputs, Ec = "EcI")
}
out$vot <- vot$summary
if (!lc) {
  pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs, prediction_settings = list(silent = TRUE))
  out$pred_head <- head(pred, 10)
  if (!eq1like) {
    cr <- mtuem_get_corr(model, apollo_probabilities, apollo_inputs)
    out$corr <- list(covar = cr$covar, corr = cr$corr, sigma = as.list(cr$sigma))
  }
}
cn <- c("chk11", "chk21", "chk22", "chk31", "chk32", "chk33")
eqn <- c("Tw", "Tf1", "Ef1")
if (case %in% c("eq3_cholesky", "lc_shared")) out$chol_dm <- mtuem_cholesky_deltaMethod(model, cn, eqn)
if (case == "lc_cholesky") out$chol_dm <- list(
  mtuem_cholesky_deltaMethod(model, paste0(cn, "_1"), eqn),
  mtuem_cholesky_deltaMethod(model, paste0(cn, "_2"), eqn))
if (case == "lc_profiled") {
  invisible(apollo_probabilities(model$estimate, apollo_inputs, functionality = "estimate"))
  pr <- mtuem_get_profiled_sigma()
  out$profiled <- list(Sigma = pr$Sigma, shares = pr$shares)
  out$ic <- mtuem_profile_ic(model, nClass = 2, nEq = 3)
}
write_json(out, paste0(case, "_ref.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE)
