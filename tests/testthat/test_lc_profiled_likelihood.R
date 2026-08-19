library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

# Two-class fixture with a single work equation plus one free time equation.
lc_fixture <- function(n = 400, seed = 1) {
  set.seed(seed)
  tau <- 168
  Tc <- stats::runif(n, 40, 90)
  Ec <- stats::runif(n, 50, 200)
  w  <- stats::runif(n, 8, 25)

  we1 <- list(Theta = 1, Phi = 1,   thw = 0)
  we2 <- list(Theta = 1, Phi = 2.5, thw = 0.5)
  te1 <- list(th1 = 0.4)
  te2 <- list(th1 = 0.15)

  tw1 <- get_tw_thph(we1, tau, Tc, Ec, w)
  tw2 <- get_tw_thph(we2, tau, Tc, Ec, w)
  opt1 <- cbind(as.numeric(tw1), as.numeric(get_ti_thph(te1, we1$Theta, tw1, tau, Tc)))
  opt2 <- cbind(as.numeric(tw2), as.numeric(get_ti_thph(te2, we2$Theta, tw2, tau, Tc)))

  S1 <- matrix(c(4.0, 1.0, 1.0, 2.25), 2)
  S2 <- matrix(c(9.0, -1.5, -1.5, 3.24), 2)
  pi1 <- 0.6
  z <- stats::rbinom(n, 1, 1 - pi1) + 1

  e1 <- matrix(stats::rnorm(n * 2), ncol = 2) %*% chol(S1)
  e2 <- matrix(stats::rnorm(n * 2), ncol = 2) %*% chol(S2)
  obs <- matrix(NA_real_, n, 2)
  obs[z == 1, ] <- (opt1 + e1)[z == 1, ]
  obs[z == 2, ] <- (opt2 + e2)[z == 2, ]

  database <- data.frame(PeID = 1:n, Tw = obs[, 1], Tf1 = obs[, 2],
                         Tc = Tc, EcI = Ec, w = w)
  class_settings <- list(
    list(work_times = "Tw", free_times = "Tf1", work_elasticities = we1,
         times_elasticities = te1, Tc = Tc, Ec = Ec, w = w, tau = tau),
    list(work_times = "Tw", free_times = "Tf1", work_elasticities = we2,
         times_elasticities = te2, Tc = Tc, Ec = Ec, w = w, tau = tau))

  list(apollo_inputs = list(database = database), database = database, n = n,
       lc_settings = list(class_settings = class_settings,
                          classProb = list(rep(pi1, n), rep(1 - pi1, n))),
       class_settings = class_settings, pi1 = pi1,
       Sigma_true = list(S1, S2))
}

test_that("the profiled likelihood handles the Apollo bookkeeping functionalities", {
  # Arrange
  fx <- lc_fixture(n = 60, seed = 71)
  apollo_inputs <- fx$apollo_inputs
  named <- fx$lc_settings; named$componentName <- "lc_component"

  # Act / Assert
  expect_equal(mtuem_lc_profiled_likelihood(fx$lc_settings, "validate"), rep(1, 60))
  expect_equal(mtuem_lc_profiled_likelihood(fx$lc_settings, "zero_LL"), rep(NA, 60))
  expect_equal(mtuem_lc_profiled_likelihood(named, "preprocess")$componentName,
               "lc_component")
  expect_equal(mtuem_lc_profiled_likelihood(fx$lc_settings, "preprocess")$componentName, "..")
  expect_error(mtuem_lc_profiled_likelihood(fx$lc_settings, "nonsense"),
               "Invalid value of argument")
})

test_that("the profiled likelihood is a proper mixture of the class densities", {
  # Arrange
  fx <- lc_fixture(n = 300, seed = 72)
  apollo_inputs <- fx$apollo_inputs

  # Act
  P <- mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate")
  prof <- mtuem_get_profiled_sigma()

  # Assert
  E1 <- mtuem_residuals(fx$class_settings[[1]], fx$apollo_inputs)
  E2 <- mtuem_residuals(fx$class_settings[[2]], fx$apollo_inputs)
  P_ref <- fx$pi1 * mtuem_dmvnorm_ref(E1, prof$Sigma[[1]]) +
    (1 - fx$pi1) * mtuem_dmvnorm_ref(E2, prof$Sigma[[2]])

  expect_equal(P, P_ref, tolerance = 1e-8)
  expect_true(all(P > 0))
})

test_that("the profiled likelihood is path independent across repeated calls", {
  # Arrange: the inner fixed point must restart from the deterministic
  # unweighted covariance every time, otherwise Apollo's numerical gradients
  # would be evaluated on a moving objective.
  fx <- lc_fixture(n = 200, seed = 73)
  apollo_inputs <- fx$apollo_inputs

  # Act
  P_first  <- mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate")
  invisible(mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate"))
  P_third  <- mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate")

  # Assert
  expect_identical(P_first, P_third)
})

test_that("the fixed point converges on a stationary covariance", {
  # Arrange
  fx <- lc_fixture(n = 300, seed = 74)
  apollo_inputs <- fx$apollo_inputs

  # Act
  invisible(mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate"))
  prof <- mtuem_get_profiled_sigma()

  # Assert: one more weighted second-moment update leaves Sigma unchanged.
  E <- lapply(fx$class_settings, mtuem_residuals, apollo_inputs = fx$apollo_inputs)
  W <- prof$posterior
  for (s in 1:2) {
    ws <- W[, s]
    Sigma_next <- crossprod(E[[s]] * sqrt(ws)) / sum(ws)
    expect_equal(unname(Sigma_next), unname(prof$Sigma[[s]]), tolerance = 1e-7)
  }
})

test_that("mtuem_get_profiled_sigma reports sigmas, correlations, posteriors and shares", {
  # Arrange
  fx <- lc_fixture(n = 300, seed = 75)
  apollo_inputs <- fx$apollo_inputs

  # Act
  invisible(mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate"))
  prof <- mtuem_get_profiled_sigma()

  # Assert
  for (s in 1:2) {
    expect_equal(prof$sig[[s]], sqrt(diag(prof$Sigma[[s]])), tolerance = 1e-12)
    expect_equal(unname(prof$rho[[s]]), unname(stats::cov2cor(prof$Sigma[[s]])),
                 tolerance = 1e-10)
    expect_equal(diag(prof$rho[[s]]), rep(1, 2), tolerance = 1e-12)
  }
  expect_equal(dim(prof$posterior), c(300L, 2L))
  expect_equal(rowSums(prof$posterior), rep(1, 300), tolerance = 1e-12)
  expect_equal(prof$shares, c(fx$pi1, 1 - fx$pi1), tolerance = 1e-12)
})

test_that("a single-class specification reduces to the plain multivariate normal", {
  # Arrange
  fx <- lc_fixture(n = 200, seed = 76)
  apollo_inputs <- fx$apollo_inputs
  one_class <- list(class_settings = fx$class_settings[1],
                    classProb = list(rep(1, 200)))

  # Act
  P <- mtuem_lc_profiled_likelihood(one_class, "estimate")

  # Assert: with a single class the posterior weights are all one, so the
  # profiled covariance is just the (uncentred) residual second moment.
  E <- mtuem_residuals(fx$class_settings[[1]], fx$apollo_inputs)
  Sigma <- crossprod(E) / nrow(E)
  expect_equal(P, mtuem_dmvnorm_ref(E, Sigma), tolerance = 1e-8)
})

test_that("prediction returns per-class allocations and posterior weights", {
  # Arrange
  fx <- lc_fixture(n = 100, seed = 77)
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- mtuem_lc_profiled_likelihood(fx$lc_settings, "prediction")

  # Assert
  expect_equal(nrow(pred), 100L)
  expect_true(all(c("post_class1", "post_class2") %in% colnames(pred)))
  expect_true(all(grepl("_class1$", colnames(pred)[1:6])))
  expect_equal(rowSums(pred[, c("post_class1", "post_class2")]), rep(1, 100),
               tolerance = 1e-12)
})

test_that("output matches estimate and report prints per-class values of time", {
  # Arrange
  fx <- lc_fixture(n = 80, seed = 78)
  apollo_inputs <- fx$apollo_inputs

  # Act / Assert
  expect_equal(mtuem_lc_profiled_likelihood(fx$lc_settings, "output"),
               mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate"))
  expect_output(mtuem_lc_profiled_likelihood(fx$lc_settings, "report"), "Class 2 VoL")
})

test_that("em_settings control the fixed point and change the answer when truncated", {
  # Arrange
  fx <- lc_fixture(n = 200, seed = 79)
  apollo_inputs <- fx$apollo_inputs
  truncated <- fx$lc_settings
  truncated$em_settings <- list(maxIterations = 1, tol = 0)

  # Act
  P_conv  <- mtuem_lc_profiled_likelihood(fx$lc_settings, "estimate")
  P_trunc <- mtuem_lc_profiled_likelihood(truncated, "estimate")

  # Assert
  expect_true(all(is.finite(P_trunc)))
  expect_false(isTRUE(all.equal(P_conv, P_trunc)))
})

test_that("mtuem_profile_ic counts the profiled covariance parameters", {
  # Arrange: 10 estimated entries of which 3 are fixed, 2 classes, 2 equations
  # per class, so k = 7 + 2 * 3 = 13.
  model <- list(estimate = stats::setNames(rep(0, 10), paste0("b", 1:10)),
                apollo_fixed = c("b1", "b2", "b3"),
                maximum = -1234.5,
                nObs = 500)

  # Act
  ic <- mtuem_profile_ic(model, nClass = 2, nEq = 2)

  # Assert
  expect_equal(ic$k, 13)
  expect_equal(ic$AIC, -2 * -1234.5 + 2 * 13)
  expect_equal(ic$BIC, -2 * -1234.5 + 13 * log(500))
})

test_that("mtuem_profile_ic penalises a larger covariance more heavily", {
  # Arrange
  model <- list(estimate = stats::setNames(rep(0, 6), paste0("b", 1:6)),
                apollo_fixed = character(0), maximum = -100, nObs = 200)

  # Act
  small <- mtuem_profile_ic(model, nClass = 2, nEq = 2)
  large <- mtuem_profile_ic(model, nClass = 2, nEq = 3)

  # Assert
  expect_equal(large$k - small$k, 2 * (6 - 3))
  expect_true(large$BIC > small$BIC)
})
