library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("preprocess returns the component name, falling back to componentName2", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 10)
  apollo_inputs <- fx$apollo_inputs
  named <- fx$mtuem_settings; named$componentName <- "work_component"
  alt   <- fx$mtuem_settings; alt$componentName2 <- "fallback_component"

  # Act / Assert
  expect_equal(mtuem_likelihood(named, "preprocess")$componentName, "work_component")
  expect_equal(mtuem_likelihood(alt, "preprocess")$componentName, "fallback_component")
  expect_equal(mtuem_likelihood(fx$mtuem_settings, "preprocess")$componentName, "..")
  expect_false(mtuem_likelihood(fx$mtuem_settings, "preprocess")$gradient)
})

test_that("validate and zero_LL return the shapes Apollo expects", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 17)
  apollo_inputs <- fx$apollo_inputs

  # Act
  v <- mtuem_likelihood(fx$mtuem_settings, "validate")
  z <- mtuem_likelihood(fx$mtuem_settings, "zero_LL")

  # Assert
  expect_equal(v, rep(1, 17))
  expect_equal(z, rep(NA, 17))
})

test_that("an unknown functionality is rejected", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 5)
  apollo_inputs <- fx$apollo_inputs

  # Act / Assert
  expect_error(mtuem_likelihood(fx$mtuem_settings, "not_a_functionality"),
               "Invalid value of argument")
})

test_that("single-equation likelihood equals the univariate normal density", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 40, times_elasticities = list(), seed = 31)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(sig1 = 2.3)

  # Act
  L <- mtuem_likelihood(settings, "estimate")

  # Assert
  err <- as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"])
  expect_equal(as.numeric(L), stats::dnorm(err, sd = 2.3), tolerance = 1e-12)
})

test_that("with no sig supplied the single-equation sigma is inferred from the residuals", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 60, times_elasticities = list(), seed = 32)
  apollo_inputs <- fx$apollo_inputs

  # Act
  L <- mtuem_likelihood(fx$mtuem_settings, "estimate")

  # Assert
  err <- as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"])
  expect_equal(as.numeric(L), stats::dnorm(err, sd = stats::sd(err)), tolerance = 1e-12)
})

test_that("sig without rho gives the product of independent normal densities", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 50, seed = 33)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  sig <- c(2.1, 1.9, 2.3)
  settings$sig <- as.list(sig)

  # Act
  L <- mtuem_likelihood(settings, "estimate")

  # Assert
  err <- fx$obs - fx$opt
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(err, diag(sig^2)), tolerance = 1e-10)
})

test_that("sig and rho reproduce the full multivariate normal density", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 45, seed = 34)
  apollo_inputs <- fx$apollo_inputs
  sig <- c(2.0, 1.5, 2.5)
  # rho is supplied as the upper triangle read horizontally: r12, r13, r23.
  r12 <- 0.35; r13 <- -0.2; r23 <- 0.15
  R <- matrix(c(1, r12, r13, r12, 1, r23, r13, r23, 1), 3)
  settings <- fx$mtuem_settings
  settings$sig <- as.list(sig)
  settings$rho <- list(r12, r13, r23)

  # Act
  L <- mtuem_likelihood(settings, "estimate")

  # Assert
  Sigma <- diag(sig) %*% R %*% diag(sig)
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(fx$obs - fx$opt, Sigma), tolerance = 1e-10)
})

test_that("a Cholesky specification reproduces the multivariate normal density", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 45, seed = 35)
  apollo_inputs <- fx$apollo_inputs
  Sigma <- matrix(c(4.0, 1.2, -0.6,
                    1.2, 2.25, 0.4,
                   -0.6, 0.4, 6.25), 3, byrow = TRUE)
  settings <- fx$mtuem_settings
  settings$cholesky <- as.list(mtuem_chol_vec(Sigma))

  # Act
  L <- mtuem_likelihood(settings, "estimate")

  # Assert
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(fx$obs - fx$opt, Sigma), tolerance = 1e-10)
})

test_that("a Cholesky specification agrees with the equivalent sig/rho specification", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 45, seed = 36)
  apollo_inputs <- fx$apollo_inputs
  Sigma <- matrix(c(4.0, 1.2, -0.6,
                    1.2, 2.25, 0.4,
                   -0.6, 0.4, 6.25), 3, byrow = TRUE)
  R <- stats::cov2cor(Sigma)

  chol_settings <- fx$mtuem_settings
  chol_settings$cholesky <- as.list(mtuem_chol_vec(Sigma))

  sigrho_settings <- fx$mtuem_settings
  sigrho_settings$sig <- as.list(sqrt(diag(Sigma)))
  sigrho_settings$rho <- list(R[1, 2], R[1, 3], R[2, 3])

  # Act
  L_chol <- mtuem_likelihood(chol_settings, "estimate")
  L_sig  <- mtuem_likelihood(sigrho_settings, "estimate")

  # Assert
  expect_equal(L_chol, L_sig, tolerance = 1e-10)
})

test_that("a Cholesky vector of the wrong length is rejected", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 10, seed = 37)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$cholesky <- as.list(c(1, 0.2, 1))  # 3 entries for a 3-equation model

  # Act / Assert
  expect_error(mtuem_likelihood(settings, "estimate"), "must match M")
})

test_that("with neither sig nor rho the covariance is inferred from the residuals", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 300, seed = 38)
  apollo_inputs <- fx$apollo_inputs

  # Act
  L <- mtuem_likelihood(fx$mtuem_settings, "estimate")

  # Assert: sigma comes from the residual covariance and the correlation from
  # the residual correlation, so the result is the normal density at that
  # empirical covariance matrix.
  err <- fx$obs - fx$opt
  sig <- sqrt(diag(stats::cov(err)))
  R <- as.matrix(Matrix::nearPD(stats::cor(err), corr = TRUE, keepDiag = TRUE)$mat)
  Sigma <- diag(sig) %*% R %*% diag(sig)
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(err, Sigma), tolerance = 1e-8)
})

test_that("a degenerate sigma yields a zero likelihood instead of an error", {
  # Arrange
  fx1 <- mtuem_fixture_thph(n = 12, times_elasticities = list(), seed = 39)
  apollo_inputs <- fx1$apollo_inputs
  one_eq <- fx1$mtuem_settings; one_eq$sig <- list(0)

  # Act
  L1 <- mtuem_likelihood(one_eq, "estimate")

  # Assert
  expect_equal(L1, rep(0, 12))
})

test_that("a non-positive-definite correlation matrix yields a zero likelihood", {
  # Arrange: r12 = r13 = 0.99 with r23 = -0.99 is not a valid correlation
  # matrix, so the eigenvalue guard must fire.
  fx <- mtuem_fixture_thph(n = 20, seed = 40)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(2, 2, 2)
  settings$rho <- list(0.99, 0.99, -0.99)

  # Act
  L <- mtuem_likelihood(settings, "estimate")

  # Assert
  expect_equal(L, rep(0, 20))
})

test_that("conditionals, raw and output return the same likelihood as estimate", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 25, seed = 41)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(2, 2, 2)

  # Act
  L_est  <- mtuem_likelihood(settings, "estimate")
  L_cond <- mtuem_likelihood(settings, "conditionals")
  L_raw  <- mtuem_likelihood(settings, "raw")
  L_out  <- mtuem_likelihood(settings, "output")

  # Assert
  expect_equal(L_cond, L_est)
  expect_equal(L_raw, L_est)
  expect_equal(L_out, L_est)
})

test_that("report prints the average values of time and still returns the likelihood", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 15, seed = 42)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(2, 2, 2)

  # Act / Assert
  expect_output(L <- mtuem_likelihood(settings, "report"), "VoL")
  expect_equal(L, mtuem_likelihood(settings, "estimate"))
})

test_that("prediction returns the optimal allocations and closes both budgets", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 30,
                           times_elasticities = list(th1 = 0.3, th2 = 0.2),
                           goods_elasticities = list(ph1 = 0.4),
                           seed = 43)
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- mtuem_likelihood(fx$mtuem_settings, "prediction")

  # Assert
  expect_equal(colnames(pred), c("Tw", "Tf1", "Tf2", "Tfi", "Ef1", "Xfj", "VoL", "VTAW"))
  expect_equal(pred[, "Tw"] + pred[, "Tf1"] + pred[, "Tf2"] + pred[, "Tfi"],
               fx$tau - fx$Tc, tolerance = 1e-9)
  expect_equal(pred[, "Ef1"] + pred[, "Xfj"], fx$w * pred[, "Tw"] - fx$Ec,
               tolerance = 1e-9)
  expect_equal(pred[, "VoL"], fx$w + pred[, "VTAW"], tolerance = 1e-9)
})

test_that("prediction reports aggregate T and X when no free categories are given", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 20, times_elasticities = list(), seed = 44)
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- mtuem_likelihood(fx$mtuem_settings, "prediction")

  # Assert
  expect_equal(colnames(pred), c("Tw", "T", "X", "VoL", "VTAW"))
  expect_equal(pred[, "Tw"] + pred[, "T"], fx$tau - fx$Tc, tolerance = 1e-9)
  expect_equal(pred[, "X"], fx$w * pred[, "Tw"] - fx$Ec, tolerance = 1e-9)
})

test_that("prediction agrees with the residual computation used for estimation", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 20, seed = 45)
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- mtuem_likelihood(fx$mtuem_settings, "prediction")
  err  <- mtuem_residuals(fx$mtuem_settings, fx$apollo_inputs)

  # Assert
  eqs <- c("Tw", "Tf1", "Tf2")
  expect_equal(unname(pred[, eqs]), unname(as.matrix(fx$database[, eqs]) - err),
               tolerance = 1e-10)
})

test_that("the gradient functionality is a documented stub", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 10, seed = 46)
  apollo_inputs <- fx$apollo_inputs

  # Act
  g <- mtuem_likelihood(fx$mtuem_settings, "gradient")

  # Assert: gradients are not implemented, so Apollo must fall back to
  # numeric derivatives.
  expect_equal(g, list(like = NA, grad = NA))
})
