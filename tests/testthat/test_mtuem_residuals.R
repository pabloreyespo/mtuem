library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("mtuem_residuals returns observed minus optimal for every equation", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 25, seed = 12)

  # Act
  err <- mtuem_residuals(fx$mtuem_settings, fx$apollo_inputs)

  # Assert
  expect_equal(dim(err), dim(fx$opt))
  expect_equal(colnames(err), c("Tw", "Tf1", "Tf2"))
  expect_equal(unname(err), unname(fx$obs - fx$opt), tolerance = 1e-10)
})

test_that("mtuem_residuals fills in optional settings with their defaults", {
  # Arrange: a work-only specification omitting free_times, free_goods,
  # goods_cost, the elasticity lists and optimal_tw entirely.
  fx <- mtuem_fixture_thph(n = 15, times_elasticities = list(), seed = 3)
  minimal <- list(work_times = "Tw",
                  work_elasticities = fx$mtuem_settings$work_elasticities,
                  Tc = fx$Tc, Ec = fx$Ec, w = fx$w, tau = fx$tau)

  # Act
  err <- mtuem_residuals(minimal, fx$apollo_inputs)

  # Assert
  expect_equal(dim(err), c(15L, 1L))
  expect_equal(colnames(err), "Tw")
  expect_equal(as.numeric(err), as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"]),
               tolerance = 1e-10)
})

test_that("mtuem_residuals covers free goods as well as free times", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 30,
                           times_elasticities = list(th1 = 0.3),
                           goods_elasticities = list(ph1 = 0.4),
                           seed = 6)

  # Act
  err <- mtuem_residuals(fx$mtuem_settings, fx$apollo_inputs)

  # Assert
  expect_equal(colnames(err), c("Tw", "Tf1", "Ef1"))
  expect_equal(unname(err), unname(fx$obs - fx$opt), tolerance = 1e-10)
})

test_that("mtuem_residuals uses obs_matrix in preference to the database", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 20, seed = 8)
  shifted <- fx$obs + 5
  inputs_with_matrix <- list(database = fx$database, obs_matrix = shifted)

  # Act
  err_db  <- mtuem_residuals(fx$mtuem_settings, fx$apollo_inputs)
  err_mat <- mtuem_residuals(fx$mtuem_settings, inputs_with_matrix)

  # Assert
  expect_equal(unname(err_mat - err_db), matrix(5, nrow = 20, ncol = 3),
               tolerance = 1e-10)
})

test_that("optimal_tw = FALSE derives the other equations from observed work time", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 20, times_elasticities = list(th1 = 0.3), seed = 15)
  settings <- fx$mtuem_settings
  settings$optimal_tw <- FALSE

  # Act
  err <- mtuem_residuals(settings, fx$apollo_inputs)

  # Assert: the free-time equation is now built from the observed Tw, while
  # the work equation still compares against the optimal Tw.
  ti_from_obs <- get_ti_thph(settings$times_elasticities,
                             settings$work_elasticities$Theta,
                             fx$database$Tw, fx$tau, fx$Tc)
  expect_equal(as.numeric(err[, "Tf1"]),
               as.numeric(fx$database$Tf1 - ti_from_obs), tolerance = 1e-10)
  expect_equal(as.numeric(err[, "Tw"]),
               as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"]), tolerance = 1e-10)
})

test_that("mtuem_residuals is exactly zero when observations equal the optimum", {
  # Arrange
  fx <- mtuem_fixture_thph(n = 12, seed = 2, sd = 0)

  # Act
  err <- mtuem_residuals(fx$mtuem_settings, fx$apollo_inputs)

  # Assert
  expect_equal(unname(err), matrix(0, nrow = 12, ncol = 3), tolerance = 1e-10)
})
