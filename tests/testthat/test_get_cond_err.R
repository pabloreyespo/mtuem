library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

make_rho <- function(M, seed = 1) {
  set.seed(seed)
  A <- matrix(stats::rnorm(M * M), M)
  stats::cov2cor(A %*% t(A) + 2 * diag(M))
}

test_that("get_cond_err reproduces the multivariate quadratic form and log determinant", {
  # The sequential conditioning is what turns the multivariate normal density
  # into a product of univariate ones inside the likelihoods, so it must
  # reproduce both mu' rho^-1 mu and log|rho| for any number of equations.
  for (M in 2:6) {
    # Arrange
    rho <- make_rho(M, seed = M)
    mu  <- matrix(stats::rnorm(10 * M), ncol = M)

    # Act
    cond <- get_cond_err(mu, rho)

    # Assert
    expect_equal(rowSums(cond$cond_mu^2),
                 rowSums((mu %*% solve(rho)) * mu),
                 tolerance = 1e-8, info = paste("M =", M))
    expect_equal(sum(log(cond$cond_sd)),
                 as.numeric(determinant(rho, logarithm = TRUE)$modulus),
                 tolerance = 1e-8, info = paste("M =", M))
  }
})

test_that("get_cond_err leaves standardised errors untouched under independence", {
  # Arrange
  M <- 4
  rho <- diag(M)
  mu <- matrix(stats::rnorm(7 * M), ncol = M)

  # Act
  cond <- get_cond_err(mu, rho)

  # Assert
  expect_equal(unname(cond$cond_mu), unname(mu), tolerance = 1e-12)
  expect_equal(cond$cond_sd, rep(1, M), tolerance = 1e-12)
})

test_that("get_cond_err matches the textbook bivariate conditioning", {
  # Arrange
  r <- 0.6
  rho <- matrix(c(1, r, r, 1), 2)
  mu <- matrix(c(1.2, -0.4, 0.9, 2.1), ncol = 2)

  # Act
  cond <- get_cond_err(mu, rho)

  # Assert
  expect_equal(cond$cond_sd, c(1, 1 - r^2), tolerance = 1e-12)
  expect_equal(unname(cond$cond_mu[, 1]), mu[, 1], tolerance = 1e-12)
  expect_equal(unname(cond$cond_mu[, 2]), (mu[, 2] - r * mu[, 1]) / sqrt(1 - r^2),
               tolerance = 1e-12)
})

test_that("get_cond_err preserves the column names of the standardised errors", {
  # Arrange
  rho <- make_rho(3, seed = 21)
  mu <- matrix(stats::rnorm(15), ncol = 3, dimnames = list(NULL, c("Tw", "Tf1", "Ef1")))

  # Act
  cond <- get_cond_err(mu, rho)

  # Assert
  expect_equal(colnames(cond$cond_mu), c("Tw", "Tf1", "Ef1"))
  expect_equal(dim(cond$cond_mu), dim(mu))
})

test_that("get_cond_err is invariant to the equation ordering, in the density it implies", {
  # Arrange: conditioning proceeds in column order, so a permutation changes
  # the intermediate conditional means but must leave the implied joint
  # density unchanged.
  M <- 3
  rho <- make_rho(M, seed = 5)
  mu  <- matrix(stats::rnorm(9 * M), ncol = M)
  perm <- c(3, 1, 2)

  # Act
  a <- get_cond_err(mu, rho)
  b <- get_cond_err(mu[, perm], rho[perm, perm])

  # Assert
  ll_a <- -0.5 * rowSums(a$cond_mu^2) - 0.5 * sum(log(a$cond_sd))
  ll_b <- -0.5 * rowSums(b$cond_mu^2) - 0.5 * sum(log(b$cond_sd))
  expect_equal(ll_a, ll_b, tolerance = 1e-8)
})
