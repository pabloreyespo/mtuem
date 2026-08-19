library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

aq_params <- function() {
  list(
    times_elasticities = list(theta_1 = 800, theta_i = 3000),
    times_satiations   = list(beta_1 = 20,  beta_i = 300),
    goods_elasticities = list(phi_1 = 100, phi_i = 80),
    goods_satiations   = list(eta_1 = 1,   eta_i = 0.5),
    goods_cost         = list(P1 = 1, Pi = 2)
  )
}

test_that("get_additive_quadratic_essentials computes the four sufficient statistics", {
  # Arrange
  p <- aq_params()
  theta <- unlist(p$times_elasticities); beta <- unlist(p$times_satiations)
  phi <- unlist(p$goods_elasticities);   eta  <- unlist(p$goods_satiations)
  P <- unlist(p$goods_cost)

  # Act
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)

  # Assert
  expect_equal(ess$Sth, sum(theta / beta), tolerance = 1e-12)
  expect_equal(ess$B,   sum(1 / beta),     tolerance = 1e-12)
  expect_equal(ess$Sph, sum(P * phi / eta), tolerance = 1e-12)
  expect_equal(ess$H,   sum(P^2 / eta),     tolerance = 1e-12)
})

test_that("get_additive_quadratic_essentials recycles a scalar price across all goods", {
  # Arrange
  p <- aq_params()

  # Act
  scalar <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                              p$goods_elasticities, p$goods_satiations,
                                              goods_cost = 1)
  vector <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                              p$goods_elasticities, p$goods_satiations,
                                              goods_cost = c(1, 1))

  # Assert
  expect_equal(scalar, vector)
})

test_that("get_tw_additive_quadratic matches the analytical work-time expression", {
  # Arrange
  p <- aq_params()
  tau <- 168; Tc <- 80; Ec <- 200; w <- 15
  theta_w <- 0; beta_w <- 10
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)

  # Act
  tw <- get_tw_additive_quadratic(theta_w, beta_w, ess$Sth, ess$Sph, ess$B, ess$H,
                                  tau, Tc, Ec, w)

  # Assert
  expected <- (theta_w + w * (ess$Sph + Ec) / ess$H - (ess$Sth - (tau - Tc)) / ess$B) /
    (beta_w + w^2 / ess$H + 1 / ess$B)
  expect_equal(dim(tw), c(1L, 1L))
  expect_equal(as.numeric(tw), expected, tolerance = 1e-9)
})

test_that("optimal work time zeroes the work first-order condition", {
  # Arrange
  p <- aq_params()
  tau <- 168; Tc <- c(70, 80, 60); Ec <- c(150, 200, 120); w <- c(12, 15, 20)
  theta_w <- 5; beta_w <- 10
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)

  # Act
  tw <- as.numeric(get_tw_additive_quadratic(theta_w, beta_w, ess$Sth, ess$Sph,
                                             ess$B, ess$H, tau, Tc, Ec, w))

  # Assert: theta_w - beta_w*Tw - mu + w*lambda = 0, with mu and lambda the
  # time and money shadow prices implied by that Tw.
  mu     <- (ess$Sth - (tau - tw - Tc)) / ess$B
  lambda <- (ess$Sph - (w * tw - Ec)) / ess$H
  expect_equal(theta_w - beta_w * tw - mu + w * lambda, rep(0, 3), tolerance = 1e-9)
})

test_that("optimal times and goods satisfy the time and money budgets exactly", {
  # Arrange
  p <- aq_params()
  tau <- 168; Tc <- c(70, 65); Ec <- c(150, 180); w <- c(12, 16)
  theta_w <- 0; beta_w <- 8
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)
  tw <- get_tw_additive_quadratic(theta_w, beta_w, ess$Sth, ess$Sph, ess$B, ess$H,
                                  tau, Tc, Ec, w)

  # Act
  ti <- get_ti_additive_quadratic(p$times_elasticities, p$times_satiations,
                                  ess$Sth, ess$B, tau, tw, Tc)
  xi <- get_xi_additive_quadratic(p$goods_elasticities, p$goods_satiations,
                                  p$goods_cost, ess$Sph, ess$H, tw, Ec, w)

  # Assert: times exhaust tau - Tc, and price-weighted consumption exhausts
  # disposable income (the goods equations return quantities, not expenses).
  expect_equal(rowSums(ti) + as.numeric(tw), tau - Tc, tolerance = 1e-9)
  expect_equal(as.numeric(xi %*% unlist(p$goods_cost)), w * as.numeric(tw) - Ec,
               tolerance = 1e-9)
})

test_that("additive quadratic allocations decrease linearly in the shadow prices", {
  # Arrange
  p <- aq_params()
  tau <- 168; Tc <- 70; Ec <- 150; w <- 12
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)
  Tw <- matrix(c(40, 45), ncol = 1)

  # Act
  ti <- get_ti_additive_quadratic(p$times_elasticities, p$times_satiations,
                                  ess$Sth, ess$B, tau, Tw, Tc)
  xi <- get_xi_additive_quadratic(p$goods_elasticities, p$goods_satiations,
                                  p$goods_cost, ess$Sph, ess$H, Tw, Ec, w)

  # Assert
  theta <- unlist(p$times_elasticities); beta <- unlist(p$times_satiations)
  phi <- unlist(p$goods_elasticities);   eta  <- unlist(p$goods_satiations)
  P <- unlist(p$goods_cost)
  mu     <- (ess$Sth - (tau - as.numeric(Tw) - Tc)) / ess$B
  lambda <- (ess$Sph - (w * as.numeric(Tw) - Ec)) / ess$H

  expect_equal(dim(ti), c(2L, 2L))
  expect_equal(ti, outer(rep(1, 2), theta / beta) - outer(mu, 1 / beta), tolerance = 1e-12)
  expect_equal(xi, outer(rep(1, 2), phi / eta) - outer(lambda, P / eta), tolerance = 1e-12)
})

test_that("values of time equal the shadow price ratio and satisfy VoL = w + VTAW", {
  # Arrange
  p <- aq_params()
  tau <- 168; Tc <- c(75, 70); Ec <- c(180, 150); w <- c(14, 18)
  theta_w <- 0; beta_w <- 11
  ess <- get_additive_quadratic_essentials(p$times_elasticities, p$times_satiations,
                                           p$goods_elasticities, p$goods_satiations,
                                           p$goods_cost)
  tw <- get_tw_additive_quadratic(theta_w, beta_w, ess$Sth, ess$Sph, ess$B, ess$H,
                                  tau, Tc, Ec, w)

  # Act
  vot <- get_values_of_time_additive_quadratic(ess$Sth, ess$Sph, ess$B, ess$H,
                                               tw, tau, Tc, Ec, w)

  # Assert
  mu     <- (ess$Sth - (tau - as.numeric(tw) - Tc)) / ess$B
  lambda <- (ess$Sph - (w * as.numeric(tw) - Ec)) / ess$H
  expect_true(all(is.finite(vot)))
  expect_equal(as.numeric(vot[, 1]), mu / lambda, tolerance = 1e-12)
  expect_equal(as.numeric(vot[, 1]), w + as.numeric(vot[, 2]), tolerance = 1e-12)
})
