library(testthat)
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("mtuem_residuals refactor preserves mtuem_likelihood output (closed-form check)", {
  # Arrange: build a small multi-equation fixture and independently derive the
  # residuals and likelihood from the closed-form Theta-Phi equations, without
  # calling mtuem_residuals itself.
  n <- 8
  tau <- 168
  Tc <- rep(60, n)
  Ec <- rep(100, n)
  w  <- rep(15, n)

  work_elasticities <- list(Theta = 1, Phi = 1, thw = 0)
  times_elasticities <- list(th1 = 0.3, th2 = 0.2)

  tw_opt <- get_tw_thph(work_elasticities, tau, Tc, Ec, w)
  ti_opt <- get_ti_thph(times_elasticities, work_elasticities$Theta, tw_opt, tau, Tc)

  set.seed(2)
  eps <- matrix(rnorm(n * 3, sd = 2), ncol = 3)
  Tw  <- as.numeric(tw_opt) + eps[, 1]
  Tf1 <- ti_opt[, 1] + eps[, 2]
  Tf2 <- ti_opt[, 2] + eps[, 3]

  database <- data.frame(PeID = 1:n, Tw = Tw, Tf1 = Tf1, Tf2 = Tf2, Tc = Tc, EcI = Ec, w = w)
  apollo_inputs <- list(database = database)

  sig <- list(sig1 = 2.1, sig2 = 1.9, sig3 = 2.3)

  mtuem_settings <- list(
    work_times = "Tw",
    free_times = c("Tf1", "Tf2"),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    sig = sig,
    Tc = Tc, Ec = Ec, w = w, tau = tau
  )

  # Act
  L_actual <- mtuem_likelihood(mtuem_settings, functionality = "estimate")

  # With an explicit sig and no rho, the correlation matrix defaults to the
  # identity, so the multi-equation likelihood reduces exactly to a product
  # of independent normal densities. This is the reference used to check that
  # the refactor into mtuem_residuals did not change the returned likelihood.
  opt <- cbind(as.numeric(tw_opt), ti_opt)
  obs <- as.matrix(database[, c("Tw", "Tf1", "Tf2")])
  err <- obs - opt
  sig_vec <- unlist(sig)
  mu <- t(t(err) / sig_vec)
  ll_expected <- -0.5 * rowSums(mu^2) - sum(log(sig_vec)) - 0.5 * 3 * log(2 * base::pi)
  L_expected <- exp(ll_expected)

  # Assert
  expect_equal(as.numeric(L_actual), as.numeric(L_expected), tolerance = 1e-10)
})

test_that("Scaled Cholesky (lambda*L) preserves correlation and scales sigma (Option 2)", {
  # Arrange
  L <- matrix(c(
    2.0, 0,   0,
    0.5, 1.5, 0,
    0.3, 0.7, 1.2
  ), nrow = 3, byrow = TRUE)
  lambda <- 1.7

  # Act
  Sigma     <- L %*% t(L)
  Sigma_lam <- (lambda * L) %*% t(lambda * L)

  sig     <- sqrt(diag(Sigma))
  sig_lam <- sqrt(diag(Sigma_lam))

  rho     <- stats::cov2cor(Sigma)
  rho_lam <- stats::cov2cor(Sigma_lam)

  # Assert
  expect_equal(sig_lam, lambda * sig, tolerance = 1e-12)
  expect_equal(rho_lam, rho, tolerance = 1e-12)
})

test_that("Row-scaled Cholesky (diag(d) L) preserves correlation structure (Option 2, per-equation variant)", {
  # Arrange
  L <- matrix(c(
    2.0, 0,   0,
    0.5, 1.5, 0,
    0.3, 0.7, 1.2
  ), nrow = 3, byrow = TRUE)
  d <- c(1.3, 0.8, 2.1)
  D <- diag(d)

  # Act
  Sigma   <- L %*% t(L)
  Sigma_d <- D %*% Sigma %*% D

  rho   <- stats::cov2cor(Sigma)
  rho_d <- stats::cov2cor(Sigma_d)

  # Assert
  expect_equal(rho_d, rho, tolerance = 1e-12)
})

test_that("Row-major Cholesky extraction/reconstruction round-trip recovers Sigma", {
  # Arrange
  M <- 3
  a    <- c(1, -0.4, 0.6)
  sigw <- 2.5
  psi  <- c(1.1, 0.9, 1.3)

  Sigma <- sigw^2 * (a %*% t(a)) + diag(psi^2)
  L <- t(chol(Sigma))

  # Act: extract the lower triangle in row-major order. This pins the
  # L11, L21, L22, L31, L32, L33 ordering that mtuem_likelihood expects.
  rc <- do.call(rbind, lapply(1:M, function(r) cbind(r, 1:r)))
  cholesky <- L[rc]

  # Reconstruct exactly as mtuem_likelihood does (R/mtuem_likelihood.R lines 150-158)
  L_mat <- matrix(0, nrow = M, ncol = M)
  idx <- 1
  for (r in 1:M) {
    for (c in 1:r) {
      L_mat[r, c] <- cholesky[idx]
      idx <- idx + 1
    }
  }
  Sigma_rebuilt <- L_mat %*% t(L_mat)

  # Assert
  expect_equal(Sigma_rebuilt, Sigma, tolerance = 1e-10)
})

test_that("Profiled likelihood recovers class-specific Sigma and matches the ordinary mixture likelihood", {
  # Arrange: simulate a 2-class cross-section with known class-specific Sigma
  # and known (constant) class membership probabilities.
  set.seed(123)
  N <- 5000

  tau <- 168
  Tc <- stats::runif(N, 40, 90)
  w  <- stats::runif(N, 8, 25)
  Ec <- stats::runif(N, 50, 200)

  work_elasticities_1  <- list(Theta = 1, Phi = 1,   thw = 0)
  work_elasticities_2  <- list(Theta = 1, Phi = 2.5, thw = 0.5)
  times_elasticities_1 <- list(th1 = 0.4)
  times_elasticities_2 <- list(th1 = 0.15)

  tw_opt_1 <- get_tw_thph(work_elasticities_1, tau, Tc, Ec, w)
  ti_opt_1 <- get_ti_thph(times_elasticities_1, work_elasticities_1$Theta, tw_opt_1, tau, Tc)
  tw_opt_2 <- get_tw_thph(work_elasticities_2, tau, Tc, Ec, w)
  ti_opt_2 <- get_ti_thph(times_elasticities_2, work_elasticities_2$Theta, tw_opt_2, tau, Tc)

  Sigma_true_1 <- matrix(c(4.0, 1.0, 1.0, 2.25), nrow = 2)
  Sigma_true_2 <- matrix(c(9.0, -1.5, -1.5, 3.24), nrow = 2)

  pi1_true <- 0.6
  z <- rbinom(N, 1, 1 - pi1_true) + 1  # 1 or 2

  chol_1 <- chol(Sigma_true_1)  # upper triangular
  chol_2 <- chol(Sigma_true_2)
  e1 <- matrix(rnorm(N * 2), ncol = 2) %*% chol_1
  e2 <- matrix(rnorm(N * 2), ncol = 2) %*% chol_2

  opt_1 <- cbind(as.numeric(tw_opt_1), as.numeric(ti_opt_1))
  opt_2 <- cbind(as.numeric(tw_opt_2), as.numeric(ti_opt_2))

  obs <- matrix(NA_real_, N, 2)
  obs[z == 1, ] <- (opt_1 + e1)[z == 1, ]
  obs[z == 2, ] <- (opt_2 + e2)[z == 2, ]

  database <- data.frame(PeID = 1:N, Tw = obs[, 1], Tf1 = obs[, 2], Tc = Tc, EcI = Ec, w = w)
  apollo_inputs <- list(database = database)

  class_settings <- list(
    list(work_times = "Tw", free_times = "Tf1",
         work_elasticities = work_elasticities_1, times_elasticities = times_elasticities_1,
         Tc = Tc, Ec = Ec, w = w, tau = tau),
    list(work_times = "Tw", free_times = "Tf1",
         work_elasticities = work_elasticities_2, times_elasticities = times_elasticities_2,
         Tc = Tc, Ec = Ec, w = w, tau = tau)
  )
  classProb <- list(rep(pi1_true, N), rep(1 - pi1_true, N))
  lc_settings <- list(class_settings = class_settings, classProb = classProb)

  # Act
  P_profiled <- mtuem_lc_profiled_likelihood(lc_settings, functionality = "estimate")

  # Assert: (i) recovered Sigma is close to the true Sigma
  expect_true(all(is.finite(P_profiled)))
  expect_true(all(P_profiled > 0))

  profiled <- mtuem_get_profiled_sigma()
  Sigma_hat_1 <- profiled$Sigma[[1]]
  Sigma_hat_2 <- profiled$Sigma[[2]]

  rel_err_1 <- max(abs(Sigma_hat_1 - Sigma_true_1)) / max(abs(Sigma_true_1))
  rel_err_2 <- max(abs(Sigma_hat_2 - Sigma_true_2)) / max(abs(Sigma_true_2))
  expect_true(rel_err_1 < 0.2)
  expect_true(rel_err_2 < 0.2)

  # Assert: (ii) the profiled likelihood equals the ordinary mixture
  # likelihood evaluated at the recovered Sigma, computed independently here
  # via solve()/determinant() rather than the internal Cholesky-based density.
  dmvnorm_manual <- function(E, Sigma) {
    M <- ncol(E)
    Sinv <- solve(Sigma)
    logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
    q <- rowSums((E %*% Sinv) * E)
    exp(-0.5 * (q + logdet + M * log(2 * base::pi)))
  }

  E1 <- mtuem_residuals(class_settings[[1]], apollo_inputs)
  E2 <- mtuem_residuals(class_settings[[2]], apollo_inputs)

  dens1 <- dmvnorm_manual(E1, Sigma_hat_1)
  dens2 <- dmvnorm_manual(E2, Sigma_hat_2)
  P_manual <- pi1_true * dens1 + (1 - pi1_true) * dens2

  expect_equal(P_profiled, P_manual, tolerance = 1e-8)
})

test_that("Profiled likelihood returns zeros instead of erroring on non-finite residuals", {
  # Arrange: w = 0 forces a division by zero inside get_tw_thph, producing
  # non-finite residuals.
  N <- 10
  tau <- 168
  Tc <- rep(50, N)
  Ec <- rep(100, N)
  w  <- rep(0, N)

  work_elasticities <- list(Theta = 1, Phi = 1, thw = 0)
  database <- data.frame(PeID = 1:N, Tw = rep(50, N), Tc = Tc, EcI = Ec, w = w)
  apollo_inputs <- list(database = database)

  cs <- list(work_times = "Tw", work_elasticities = work_elasticities, Tc = Tc, Ec = Ec, w = w, tau = tau)
  class_settings <- list(cs, cs)
  classProb <- list(rep(0.5, N), rep(0.5, N))
  lc_settings <- list(class_settings = class_settings, classProb = classProb)

  # Act
  P <- mtuem_lc_profiled_likelihood(lc_settings, functionality = "estimate")

  # Assert
  expect_equal(P, rep(0, N))
})

test_that("Profiled likelihood returns zeros instead of erroring on a non-PD starting Sigma", {
  # Arrange: identical noise added to both equations makes the residuals
  # perfectly collinear, so the unweighted starting covariance is singular.
  N <- 30
  tau <- 168
  set.seed(7)
  Tc <- stats::runif(N, 40, 80)
  Ec <- stats::runif(N, 50, 150)
  w  <- stats::runif(N, 10, 20)

  work_elasticities <- list(Theta = 1, Phi = 1, thw = 0)
  times_elasticities <- list(th1 = 0.3)

  tw_opt <- get_tw_thph(work_elasticities, tau, Tc, Ec, w)
  ti_opt <- get_ti_thph(times_elasticities, work_elasticities$Theta, tw_opt, tau, Tc)

  noise <- rnorm(N, sd = 3)
  Tw  <- as.numeric(tw_opt) + noise
  Tf1 <- as.numeric(ti_opt) + noise

  database <- data.frame(PeID = 1:N, Tw = Tw, Tf1 = Tf1, Tc = Tc, EcI = Ec, w = w)
  apollo_inputs <- list(database = database)

  cs <- list(work_times = "Tw", free_times = "Tf1",
             work_elasticities = work_elasticities, times_elasticities = times_elasticities,
             Tc = Tc, Ec = Ec, w = w, tau = tau)
  class_settings <- list(cs, cs)
  classProb <- list(rep(0.5, N), rep(0.5, N))
  lc_settings <- list(class_settings = class_settings, classProb = classProb)

  # Act
  P <- mtuem_lc_profiled_likelihood(lc_settings, functionality = "estimate")

  # Assert
  expect_equal(P, rep(0, N))
})

test_that("Profiled likelihood returns zeros instead of erroring on class collapse", {
  # Arrange: class 2 has a vanishingly small prior everywhere, so its
  # posterior weight cannot reach the minimum weight guard.
  N <- 40
  tau <- 168
  set.seed(11)
  Tc <- stats::runif(N, 40, 80)
  Ec <- stats::runif(N, 50, 150)
  w  <- stats::runif(N, 10, 20)

  work_elasticities <- list(Theta = 1, Phi = 1, thw = 0)
  tw_opt <- get_tw_thph(work_elasticities, tau, Tc, Ec, w)
  Tw <- as.numeric(tw_opt) + rnorm(N, sd = 2)

  database <- data.frame(PeID = 1:N, Tw = Tw, Tc = Tc, EcI = Ec, w = w)
  apollo_inputs <- list(database = database)

  cs <- list(work_times = "Tw", work_elasticities = work_elasticities, Tc = Tc, Ec = Ec, w = w, tau = tau)
  class_settings <- list(cs, cs)
  classProb <- list(rep(1 - 1e-8, N), rep(1e-8, N))
  lc_settings <- list(class_settings = class_settings, classProb = classProb,
                       em_settings = list(minWeight = 5))

  # Act
  P <- mtuem_lc_profiled_likelihood(lc_settings, functionality = "estimate")

  # Assert
  expect_equal(P, rep(0, N))
})

test_that("customMultiStart normalizes shared unsuffixed covariance parameters without erroring", {
  # Arrange: a minimal 2-class, 1-equation latent class model with a single
  # covariance parameter (chk11) SHARED across both classes, i.e. not
  # suffixed with "_1"/"_2" in apollo_beta. Before the R/customMultiStart.R
  # fix, normalize() errored out while building the per-class suffixed pass
  # ("chk11_1", which does not exist in apollo_beta).
  apollo_initialise()
  apollo_control <- list(
    modelName = "test-shared-cov-normalize",
    indivID = "PeID",
    outputDirectory = tempdir(),
    noValidation = TRUE,
    noDiagnostics = TRUE
  )

  database <- mtuem::maed
  database <- database[database$EcI > 0, ][1:100, ]

  apollo_beta <- c(
    Theta_1 = 1, Phi_1 = 1, thw_1 = 0,
    asc_1 = 0, x_female_1 = 0, x_older45_1 = 0,
    Theta_2 = 1, Phi_2 = 1, thw_2 = 0,
    asc_2 = 0, x_female_2 = 0, x_older45_2 = 0,
    chk11 = 8  # single SHARED covariance parameter (M = 1 equation)
  )
  apollo_fixed <- c("Theta_1", "Theta_2", "asc_1", "x_female_1", "x_older45_1")

  apollo_lcPars <- function(apollo_beta, apollo_inputs) {
    lcpars <- list()
    lcpars[["Theta"]] <- list(Theta_1, Theta_2)
    lcpars[["Phi"]] <- list(Phi_1, Phi_2)
    lcpars[["thw"]] <- list(thw_1, thw_2)
    V <- list()
    V[["class_1"]] <- asc_1 + (female * x_female_1) + (older45 * x_older45_1)
    V[["class_2"]] <- asc_2 + (female * x_female_2) + (older45 * x_older45_2)
    classAlloc_settings <- list(classes = c(class_1 = 1, class_2 = 2), utilities = V)
    lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
    return(lcpars)
  }

  # apollo_validateInputs() looks up apollo_control, database, apollo_beta and
  # apollo_lcPars in the global environment, so they are mirrored there for
  # the duration of this test and removed again afterwards.
  globals <- list(apollo_control = apollo_control, database = database,
                   apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                   apollo_lcPars = apollo_lcPars)
  for (nm in names(globals)) assign(nm, globals[[nm]], envir = .GlobalEnv)
  on.exit(rm(list = names(globals), envir = .GlobalEnv), add = TRUE)

  apollo_inputs <- apollo_validateInputs()

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    for (s in 1:length(pi_values)) {
      work_elasticities <- list(Theta = Theta[[s]], Phi = Phi[[s]], thw = thw[[s]])
      mtuem_settings <- list(
        work_times = c("Tw"),
        work_elasticities = work_elasticities,
        cholesky = list(chk11 = chk11),
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

  # Act
  result <- tryCatch({
    customMultiStart(
      apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
      customMultistart_settings = list(nCandidates = 2),
      estimate_settings = list(maxIterations = 1, hessianRoutine = "none", estimationRoutine = "bfgs"),
      first_em = FALSE,
      verbose = FALSE,
      normalization = list(
        normalization = "theta_phi",
        work_elasticities = c("Theta", "Phi", "thw"),
        cholesky = c("chk11"),
        covs = c("asc", "x_female", "x_older45")
      ),
      nClass = 2
    )
    "ok"
  }, error = function(e) e)

  # Assert
  expect_false(inherits(result, "error"),
               info = if (inherits(result, "error")) conditionMessage(result) else NULL)
})
