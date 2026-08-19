library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("get_tw_thph returns the root of the work first-order condition", {
  # Arrange
  tau <- 168
  Tc <- c(50, 65, 70)
  Ec <- c(100, 120, 90)
  w  <- c(12, 18, 22)
  we <- list(Theta = 1.2, Phi = 0.8, thw = 0.15)

  # Act
  tw <- get_tw_thph(we, tau, Tc, Ec, w)

  # Assert: u = Tw / (tau - Tc) solves base*u^2 - A*u + thw*ec = 0, where
  # base = Theta + Phi + thw and A = Phi + thw + (Theta + thw)*ec.
  ec   <- Ec / (w * (tau - Tc))
  base <- we$Theta + we$Phi + we$thw
  A    <- we$Phi + we$thw + (we$Theta + we$thw) * ec
  u <- as.numeric(tw) / (tau - Tc)

  expect_equal(base * u^2 - A * u + we$thw * ec, rep(0, 3), tolerance = 1e-12)
})

test_that("get_tw_thph collapses to the closed form when thw = 0", {
  # Arrange
  tau <- 168; Tc <- 60; Ec <- 100; w <- 15
  we <- list(Theta = 1, Phi = 2, thw = 0)
  ec <- Ec / (w * (tau - Tc))

  # Act
  tw <- get_tw_thph(we, tau, Tc, Ec, w)

  # Assert
  expected <- (tau - Tc) * (we$Phi + we$Theta * ec) / (we$Theta + we$Phi)
  expect_equal(as.numeric(tw), expected, tolerance = 1e-12)
})

test_that("get_tw_thph returns an N x 1 matrix regardless of input length", {
  # Arrange
  we <- list(Theta = 1, Phi = 1, thw = 0)

  # Act
  tw1 <- get_tw_thph(we, 168, 60, 100, 15)
  twN <- get_tw_thph(we, 168, rep(60, 7), rep(100, 7), rep(15, 7))

  # Assert
  expect_equal(dim(tw1), c(1L, 1L))
  expect_equal(dim(twN), c(7L, 1L))
  expect_equal(as.numeric(twN), rep(as.numeric(tw1), 7), tolerance = 1e-12)
})

test_that("get_tw_thph leaves a strictly positive free-time residual", {
  # Arrange
  tau <- 168
  set.seed(4)
  Tc <- stats::runif(50, 40, 90)
  Ec <- stats::runif(50, 50, 200)
  w  <- stats::runif(50, 8, 30)
  we <- list(Theta = 1, Phi = 1, thw = 0)

  # Act
  tw <- as.numeric(get_tw_thph(we, tau, Tc, Ec, w))

  # Assert
  expect_true(all(tw > 0))
  expect_true(all(tw < tau - Tc))
  expect_true(all(w * tw - Ec > 0))
})

test_that("get_ti_thph splits residual time in proportion to the time elasticities", {
  # Arrange
  tau <- 168
  Tc <- c(60, 70)
  Theta <- 1.5
  te <- list(th1 = 0.3, th2 = 0.45)
  Tw <- matrix(c(40, 35), ncol = 1)

  # Act
  ti <- get_ti_thph(te, Theta, Tw, tau, Tc)

  # Assert
  expect_equal(dim(ti), c(2L, 2L))
  expect_equal(colnames(ti), c("th1", "th2"))
  expect_equal(ti[, "th1"], (0.3 / Theta) * (tau - as.numeric(Tw) - Tc), tolerance = 1e-12)
  expect_equal(ti[, "th2"] / ti[, "th1"], rep(0.45 / 0.3, 2), tolerance = 1e-12)
})

test_that("get_ti_thph exhausts the time budget when the elasticities sum to Theta", {
  # Arrange: the leave-one-out constraint is exactly binding here, so the
  # parameterised activities absorb the whole residual time.
  tau <- 168; Tc <- 60; Theta <- 1
  te <- list(th1 = 0.4, th2 = 0.6)
  Tw <- matrix(40, ncol = 1)

  # Act
  ti <- get_ti_thph(te, Theta, Tw, tau, Tc)

  # Assert
  expect_equal(sum(ti), tau - as.numeric(Tw) - Tc, tolerance = 1e-12)
})

test_that("get_xi_thph splits residual income by elasticity over price", {
  # Arrange
  Phi <- 2
  ge <- list(ph1 = 0.5, ph2 = 1.0)
  gc <- list(p1 = 1, p2 = 2)
  Tw <- matrix(c(40, 45), ncol = 1)
  Ec <- c(100, 120)
  w  <- c(15, 20)

  # Act
  xj <- get_xi_thph(ge, gc, Phi, Tw, Ec, w)

  # Assert
  income <- w * as.numeric(Tw) - Ec
  expect_equal(dim(xj), c(2L, 2L))
  expect_equal(xj[, 1], (0.5 / 1 / Phi) * income, tolerance = 1e-12)
  expect_equal(xj[, 2], (1.0 / 2 / Phi) * income, tolerance = 1e-12)
})

test_that("get_values_of_time_thph satisfies VoL = w + VTAW", {
  # Arrange
  tau <- 168
  Tc <- c(60, 70, 55)
  Ec <- c(100, 150, 80)
  w  <- c(15, 12, 25)
  we <- list(Theta = 1.4, Phi = 0.9, thw = 0.1)
  tw <- get_tw_thph(we, tau, Tc, Ec, w)

  # Act
  vot <- get_values_of_time_thph(tw, we, Tc, Ec, w)

  # Assert
  # Columns are returned positionally as (VoL, VTAW); the callers in
  # mtuem_likelihood attach the names themselves.
  expect_equal(ncol(vot), 2L)
  expect_equal(as.numeric(vot[, 1]), w + as.numeric(vot[, 2]), tolerance = 1e-12)
})

test_that("get_values_of_time_thph equals the marginal rate of substitution formula", {
  # Arrange
  tau <- 168; Tc <- 60; Ec <- 100; w <- 15
  we <- list(Theta = 1.4, Phi = 0.9, thw = 0.1)
  tw <- get_tw_thph(we, tau, Tc, Ec, w)

  # Act
  vot <- get_values_of_time_thph(tw, we, Tc, Ec, w)

  # Assert: VoL = (Theta/Phi) * (w*Tw - Ec) / (tau - Tw - Tc). NOTE: the
  # implementation hardcodes 168 in place of tau, so this reference uses
  # tau = 168 to match; it would diverge for any other time budget.
  expected <- (we$Theta / we$Phi) * (w * as.numeric(tw) - Ec) / (168 - as.numeric(tw) - Tc)
  expect_equal(as.numeric(vot[, 1]), expected, tolerance = 1e-12)
})

test_that("VoL is increasing in Theta/Phi, holding the allocation fixed", {
  # Arrange
  tau <- 168; Tc <- 60; Ec <- 100; w <- 15
  tw <- matrix(45, ncol = 1)

  # Act
  low  <- get_values_of_time_thph(tw, list(Theta = 1, Phi = 2, thw = 0), Tc, Ec, w)
  high <- get_values_of_time_thph(tw, list(Theta = 2, Phi = 1, thw = 0), Tc, Ec, w)

  # Assert
  expect_true(high[, 1] > low[, 1])
})
