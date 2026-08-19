library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

# Additive quadratic fixture: observed allocations are the optimal ones plus noise.
aq_fixture <- function(n = 60, seed = 1, sd = 1, goods_cost = list(P1 = 1, P2 = 1)) {
  set.seed(seed)
  tau <- 168
  Tc <- stats::runif(n, 40, 80)
  Ec <- stats::runif(n, 50, 150)
  w  <- stats::runif(n, 10, 25)

  times_elasticities <- list(theta_1 = 800, theta_2 = 3000)
  times_satiations   <- list(beta_1 = 20,  beta_2 = 300)
  goods_elasticities <- list(phi_1 = 100, phi_2 = 80)
  goods_satiations   <- list(eta_1 = 1,   eta_2 = 0.5)
  work_elasticity <- 0
  work_satiation  <- 10

  ess <- get_additive_quadratic_essentials(times_elasticities, times_satiations,
                                           goods_elasticities, goods_satiations, goods_cost)
  tw <- get_tw_additive_quadratic(work_elasticity, work_satiation, ess$Sth, ess$Sph,
                                  ess$B, ess$H, tau, Tc, Ec, w)
  ti <- get_ti_additive_quadratic(times_elasticities, times_satiations, ess$Sth, ess$B, tau, tw, Tc)
  xj <- get_xi_additive_quadratic(goods_elasticities, goods_satiations, goods_cost,
                                  ess$Sph, ess$H, tw, Ec, w)

  opt <- cbind(tw, ti, xj)
  colnames(opt) <- c("Tw", "Tf1", "Tf2", "Ef1", "Ef2")
  obs <- opt + matrix(stats::rnorm(n * 5, sd = sd), ncol = 5)
  colnames(obs) <- colnames(opt)
  database <- cbind(data.frame(PeID = 1:n, Tc = Tc, EcI = Ec, w = w), as.data.frame(obs))

  list(apollo_inputs = list(database = database), database = database,
       opt = opt, obs = obs, ess = ess, n = n, tau = tau, Tc = Tc, Ec = Ec, w = w,
       mtuem_settings = list(
         work_times = "Tw",
         free_times = c("Tf1", "Tf2"),
         free_goods = c("Ef1", "Ef2"),
         goods_cost = goods_cost,
         work_elasticity = work_elasticity,
         work_satiation = work_satiation,
         times_elasticities = times_elasticities,
         times_satiations = times_satiations,
         goods_elasticities = goods_elasticities,
         goods_satiations = goods_satiations,
         Tc = Tc, Ec = Ec, w = w, tau = tau))
}

test_that("additive quadratic likelihood reproduces the multivariate normal density", {
  # Arrange: with all five equations kept, the residual covariance is singular
  # (the budgets make them collinear), so a Cholesky factor is supplied for a
  # reduced set of equations instead.
  fx <- aq_fixture(n = 50, seed = 61)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- "Tf2"
  settings$omitted_goods <- "Ef2"

  Sigma <- matrix(c(1.5, 0.3, -0.2,
                    0.3, 1.0,  0.1,
                   -0.2, 0.1,  2.0), 3, byrow = TRUE)
  settings$cholesky <- as.list(mtuem_chol_vec(Sigma))

  # Act
  L <- additive_quadratic_mtuem_likelihood(settings, "estimate")

  # Assert
  eqs <- c("Tw", "Tf1", "Ef1")
  expect_equal(as.numeric(L),
               mtuem_dmvnorm_ref(fx$obs[, eqs] - fx$opt[, eqs], Sigma),
               tolerance = 1e-10)
})

test_that("omitting all but one equation falls back to the univariate normal", {
  # Arrange
  fx <- aq_fixture(n = 30, seed = 62)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- c("Tf1", "Tf2")
  settings$omitted_goods <- c("Ef1", "Ef2")
  settings$sig <- list(1.3)

  # Act
  L <- additive_quadratic_mtuem_likelihood(settings, "estimate")

  # Assert
  err <- as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"])
  expect_equal(as.numeric(L), stats::dnorm(err, sd = 1.3), tolerance = 1e-12)
})

test_that("omitting an equation keeps its parameters in the aggregates", {
  # Arrange: dropping Tf2 from the likelihood must not change the predicted
  # Tw, because theta_2 and beta_2 still enter S_theta and B.
  fx <- aq_fixture(n = 30, seed = 63)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- "Tf2"
  settings$omitted_goods <- "Ef2"

  # Act
  pred_full <- additive_quadratic_mtuem_likelihood(fx$mtuem_settings, "prediction")
  pred_omit <- additive_quadratic_mtuem_likelihood(settings, "prediction")

  # Assert
  expect_equal(pred_omit, pred_full, tolerance = 1e-12)
})

test_that("omitted categories must name existing free categories", {
  # Arrange
  fx <- aq_fixture(n = 10, seed = 64)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- "Tf9"

  # Act / Assert
  expect_error(additive_quadratic_mtuem_likelihood(settings, "estimate"),
               "not in free_times / free_goods")
})

test_that("omitting every equation is rejected", {
  # Arrange
  fx <- aq_fixture(n = 10, seed = 65)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$work_times <- character(0)
  settings$omitted_times <- c("Tf1", "Tf2")
  settings$omitted_goods <- c("Ef1", "Ef2")

  # Act / Assert
  expect_error(additive_quadratic_mtuem_likelihood(settings, "estimate"),
               "nothing left to estimate")
})

test_that("the Cholesky length is checked against the surviving equations", {
  # Arrange: two equations are omitted, so three Cholesky entries are needed
  # for the three that remain, not fifteen for all five.
  fx <- aq_fixture(n = 10, seed = 66)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- c("Tf1", "Tf2")
  settings$omitted_goods <- "Ef2"
  settings$cholesky <- as.list(rep(1, 15))

  # Act / Assert
  expect_error(additive_quadratic_mtuem_likelihood(settings, "estimate"), "must match M")
})

test_that("additive quadratic prediction closes both budgets and reports values of time", {
  # Arrange
  fx <- aq_fixture(n = 25, seed = 67, goods_cost = list(P1 = 1, P2 = 2))
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- additive_quadratic_mtuem_likelihood(fx$mtuem_settings, "prediction")

  # Assert: no residual time or goods category is reported, because the
  # parameterised categories already exhaust both budgets.
  expect_equal(colnames(pred), c("Tw", "Tf1", "Tf2", "Ef1", "Ef2", "VoL", "VTAW"))
  expect_equal(pred[, "Tw"] + pred[, "Tf1"] + pred[, "Tf2"], fx$tau - fx$Tc, tolerance = 1e-9)
  expect_equal(1 * pred[, "Ef1"] + 2 * pred[, "Ef2"], fx$w * pred[, "Tw"] - fx$Ec,
               tolerance = 1e-9)
  expect_equal(pred[, "VoL"], fx$w + pred[, "VTAW"], tolerance = 1e-9)
})

test_that("additive quadratic handles the Apollo bookkeeping functionalities", {
  # Arrange
  fx <- aq_fixture(n = 11, seed = 68)
  apollo_inputs <- fx$apollo_inputs
  named <- fx$mtuem_settings; named$componentName <- "aq_component"

  # Act / Assert
  expect_equal(additive_quadratic_mtuem_likelihood(fx$mtuem_settings, "validate"), rep(1, 11))
  expect_equal(additive_quadratic_mtuem_likelihood(fx$mtuem_settings, "zero_LL"), rep(NA, 11))
  expect_equal(additive_quadratic_mtuem_likelihood(named, "preprocess")$componentName,
               "aq_component")
  expect_error(additive_quadratic_mtuem_likelihood(fx$mtuem_settings, "nonsense"),
               "Invalid value of argument")
})

test_that("additive quadratic output matches estimate and report prints the values of time", {
  # Arrange
  fx <- aq_fixture(n = 20, seed = 69)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- "Tf2"
  settings$omitted_goods <- "Ef2"
  settings$sig <- list(1, 1, 1)

  # Act / Assert
  expect_equal(additive_quadratic_mtuem_likelihood(settings, "output"),
               additive_quadratic_mtuem_likelihood(settings, "estimate"))
  expect_output(additive_quadratic_mtuem_likelihood(settings, "report"), "VoL")
})

test_that("optimal_tw = FALSE builds the other equations from observed work time", {
  # Arrange
  fx <- aq_fixture(n = 25, seed = 70)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$omitted_times <- "Tf2"
  settings$omitted_goods <- c("Ef1", "Ef2")
  settings$optimal_tw <- FALSE
  settings$sig <- list(1, 1)

  # Act
  L <- additive_quadratic_mtuem_likelihood(settings, "estimate")

  # Assert
  ti_from_obs <- get_ti_additive_quadratic(settings$times_elasticities,
                                           settings$times_satiations,
                                           fx$ess$Sth, fx$ess$B, fx$tau,
                                           fx$database$Tw, fx$Tc)
  err <- cbind(fx$obs[, "Tw"] - fx$opt[, "Tw"],
               fx$obs[, "Tf1"] - ti_from_obs[, 1])
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(err, diag(2)), tolerance = 1e-10)
})
