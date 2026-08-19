library(testthat)
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("a Cholesky parameter count that is not M*(M+1)/2 is rejected", {
  # Arrange
  model <- mtuem_fake_apollo_model(c(a = 1, b = 2), diag(2) * 0.01)

  # Act / Assert
  expect_error(mtuem_cholesky_deltaMethod(model, c("a", "b")),
               "M \\* \\(M \\+ 1\\) / 2")
})

test_that("delta method recovers the sigmas and correlation implied by the Cholesky factor", {
  # Arrange
  est <- c(chk11 = 2, chk21 = 0.5, chk22 = 1.5)
  model <- mtuem_fake_apollo_model(est, diag(c(0.01, 0.02, 0.03)))

  # Act
  res <- suppressMessages(mtuem_cholesky_deltaMethod(model, names(est)))

  # Assert: Sigma = L L', so sig1 = L11, sig2 = sqrt(L21^2 + L22^2) and
  # rho = L21 * L11 / (sig1 * sig2).
  L <- matrix(c(2, 0, 0.5, 1.5), 2, byrow = TRUE)
  Sigma <- L %*% t(L)
  sig <- sqrt(diag(Sigma))
  rho <- Sigma[2, 1] / (sig[1] * sig[2])

  expect_equal(as.character(res$Expression), c("sig_1", "sig_2", "rho_1_2"))
  expect_equal(res$Value, c(sig, rho), tolerance = 1e-4)
  expect_true(all(res[["s.e."]] > 0))
})

test_that("equation names label the delta method output when supplied", {
  # Arrange
  est <- c(chk11 = 2, chk21 = 0.5, chk22 = 1.5, chk31 = 0.3, chk32 = 0.2, chk33 = 1.1)
  model <- mtuem_fake_apollo_model(est, diag(rep(0.01, 6)))

  # Act
  res <- suppressMessages(
    mtuem_cholesky_deltaMethod(model, names(est), eq_names = c("Tw", "Tf1", "Ef1")))

  # Assert
  expect_equal(as.character(res$Expression),
               c("sig_Tw", "sig_Tf1", "sig_Ef1",
                 "rho_Tw_Tf1", "rho_Tw_Ef1", "rho_Tf1_Ef1"))
})

test_that("delta method values match the covariance implied by a 3-equation factor", {
  # Arrange
  est <- c(chk11 = 2, chk21 = 0.5, chk22 = 1.5, chk31 = 0.3, chk32 = 0.2, chk33 = 1.1)
  model <- mtuem_fake_apollo_model(est, diag(rep(0.01, 6)))
  L <- matrix(c(2, 0, 0,
                0.5, 1.5, 0,
                0.3, 0.2, 1.1), 3, byrow = TRUE)
  Sigma <- L %*% t(L)
  sig <- sqrt(diag(Sigma))
  R <- stats::cov2cor(Sigma)

  # Act
  res <- suppressMessages(mtuem_cholesky_deltaMethod(model, names(est)))

  # Assert: sigmas first, then correlations in (1,2), (1,3), (2,3) order.
  expect_equal(res$Value,
               c(sig, R[2, 1], R[3, 1], R[3, 2]),
               tolerance = 1e-4)
})
