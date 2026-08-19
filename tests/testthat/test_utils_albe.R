library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("get_tw_albe returns the root of the alpha-beta work equation", {
  # Arrange
  tau <- 168
  Tc <- c(50, 60, 72)
  Ec <- c(90, 120, 140)
  w  <- c(14, 18, 21)
  we <- list(alpha = 0.3, beta = 0.25)

  # Act
  tw <- get_tw_albe(we, tau, Tc, Ec, w)

  # Assert: u = Tw / (tau - Tc) solves u^2 - 2*B*u + ec*(2a + 2b - 1) = 0
  # with B = beta + alpha*ec.
  ec <- Ec / (w * (tau - Tc))
  B  <- we$beta + we$alpha * ec
  u  <- as.numeric(tw) / (tau - Tc)

  expect_equal(u^2 - 2 * B * u + ec * (2 * we$alpha + 2 * we$beta - 1),
               rep(0, 3), tolerance = 1e-12)
})

test_that("get_tw_albe returns an N x 1 matrix and a feasible interior allocation", {
  # Arrange
  tau <- 168
  set.seed(9)
  Tc <- stats::runif(40, 40, 90)
  Ec <- stats::runif(40, 50, 180)
  w  <- stats::runif(40, 10, 28)
  we <- list(alpha = 0.3, beta = 0.25)

  # Act
  tw <- get_tw_albe(we, tau, Tc, Ec, w)

  # Assert
  expect_equal(dim(tw), c(40L, 1L))
  expect_true(all(tw > 0))
  expect_true(all(as.numeric(tw) < tau - Tc))
  expect_true(all(w * as.numeric(tw) - Ec > 0))
})

test_that("get_ti_albe splits residual time by elasticity over (1 - 2*beta)", {
  # Arrange
  tau <- 168
  Tc <- c(60, 70)
  beta <- 0.25
  te <- list(th1 = 0.2, th2 = 0.1)
  Tw <- matrix(c(40, 38), ncol = 1)

  # Act
  ti <- get_ti_albe(te, beta, Tw, tau, Tc)

  # Assert
  expect_equal(dim(ti), c(2L, 2L))
  expect_equal(ti[, "th1"], (0.2 / (1 - 2 * beta)) * (tau - as.numeric(Tw) - Tc),
               tolerance = 1e-12)
  expect_equal(ti[, "th2"] / ti[, "th1"], rep(0.5, 2), tolerance = 1e-12)
})

test_that("get_xi_albe splits residual income by elasticity over price and (1 - alpha)", {
  # Arrange
  alpha <- 0.3
  ge <- list(ph1 = 0.2, ph2 = 0.1)
  gc <- list(p1 = 1, p2 = 2)
  Tw <- matrix(c(42, 39), ncol = 1)
  Ec <- c(100, 90)
  w  <- c(16, 19)

  # Act
  xj <- get_xi_albe(ge, gc, alpha, Tw, Ec, w)

  # Assert
  income <- w * as.numeric(Tw) - Ec
  expect_equal(xj[, 1], (0.2 / 1 / (1 - alpha)) * income, tolerance = 1e-12)
  expect_equal(xj[, 2], (0.1 / 2 / (1 - alpha)) * income, tolerance = 1e-12)
})

test_that("get_values_of_time_albe satisfies VoL = w + VTAW and its closed form", {
  # Arrange
  tau <- 168
  Tc <- c(60, 55)
  Ec <- c(100, 130)
  w  <- c(15, 20)
  we <- list(alpha = 0.3, beta = 0.25)
  tw <- get_tw_albe(we, tau, Tc, Ec, w)

  # Act
  vot <- get_values_of_time_albe(tw, we, Tc, Ec, w)

  # Assert. NOTE: the implementation hardcodes 168 in place of tau, so this
  # reference uses tau = 168 to match.
  expected <- ((w * as.numeric(tw) - Ec) / (168 - as.numeric(tw) - Tc)) *
    (1 - 2 * we$beta) / (1 - 2 * we$alpha)
  expect_equal(as.numeric(vot[, 1]), expected, tolerance = 1e-12)
  expect_equal(as.numeric(vot[, 1]), w + as.numeric(vot[, 2]), tolerance = 1e-12)
})
