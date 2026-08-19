library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

# Alpha-beta fixture: observed allocations are the optimal ones plus noise.
ab_fixture <- function(n = 60, times_elasticities = list(th1 = 0.1),
                       goods_elasticities = list(), seed = 1, sd = 2) {
  set.seed(seed)
  tau <- 168
  Tc <- stats::runif(n, 40, 80)
  Ec <- stats::runif(n, 50, 150)
  w  <- stats::runif(n, 10, 25)
  we <- list(alpha = 0.3, beta = 0.25)

  tw <- get_tw_albe(we, tau, Tc, Ec, w)
  opt <- as.matrix(tw); nms <- "Tw"
  if (length(times_elasticities) > 0) {
    opt <- cbind(opt, get_ti_albe(times_elasticities, we$beta, tw, tau, Tc))
    nms <- c(nms, paste0("Tf", seq_along(times_elasticities)))
  }
  if (length(goods_elasticities) > 0) {
    opt <- cbind(opt, get_xi_albe(goods_elasticities, 1, we$alpha, tw, Ec, w))
    nms <- c(nms, paste0("Ef", seq_along(goods_elasticities)))
  }
  colnames(opt) <- nms

  obs <- opt + matrix(stats::rnorm(n * ncol(opt), sd = sd), ncol = ncol(opt))
  colnames(obs) <- nms
  database <- cbind(data.frame(PeID = 1:n, Tc = Tc, EcI = Ec, w = w), as.data.frame(obs))

  list(apollo_inputs = list(database = database), database = database,
       opt = opt, obs = obs, n = n, tau = tau, Tc = Tc, Ec = Ec, w = w,
       mtuem_settings = list(
         work_times = "Tw",
         free_times = if (length(times_elasticities) > 0) nms[2:(1 + length(times_elasticities))] else list(),
         free_goods = if (length(goods_elasticities) > 0) utils::tail(nms, length(goods_elasticities)) else list(),
         work_elasticities = we,
         times_elasticities = times_elasticities,
         goods_elasticities = goods_elasticities,
         Tc = Tc, Ec = Ec, w = w, tau = tau))
}

test_that("alpha-beta single-equation likelihood equals the univariate normal density", {
  # Arrange
  fx <- ab_fixture(n = 40, times_elasticities = list(), seed = 51)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(2.4)

  # Act
  L <- mtuem_ab_likelihood(settings, "estimate")

  # Assert
  err <- as.numeric(fx$obs[, "Tw"] - fx$opt[, "Tw"])
  expect_equal(as.numeric(L), stats::dnorm(err, sd = 2.4), tolerance = 1e-12)
})

test_that("alpha-beta Cholesky specification reproduces the multivariate normal density", {
  # Arrange
  fx <- ab_fixture(n = 50, times_elasticities = list(th1 = 0.1), seed = 52)
  apollo_inputs <- fx$apollo_inputs
  Sigma <- matrix(c(4.0, 1.1, 1.1, 2.25), 2)
  settings <- fx$mtuem_settings
  settings$cholesky <- as.list(mtuem_chol_vec(Sigma))

  # Act
  L <- mtuem_ab_likelihood(settings, "estimate")

  # Assert
  expect_equal(as.numeric(L), mtuem_dmvnorm_ref(fx$obs - fx$opt, Sigma), tolerance = 1e-10)
})

test_that("alpha-beta rejects a Cholesky vector of the wrong length", {
  # Arrange
  fx <- ab_fixture(n = 10, seed = 53)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$cholesky <- as.list(c(1, 0.2))  # 2 entries for a 2-equation model

  # Act / Assert
  expect_error(mtuem_ab_likelihood(settings, "estimate"), "must match M")
})

test_that("alpha-beta guards against a degenerate sigma", {
  # Arrange
  fx <- ab_fixture(n = 14, times_elasticities = list(), seed = 54)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(0)

  # Act / Assert
  expect_equal(mtuem_ab_likelihood(settings, "estimate"), rep(0, 14))
})

test_that("alpha-beta get_covar returns the residual covariance, correlation and sigmas", {
  # Arrange
  fx <- ab_fixture(n = 80, times_elasticities = list(th1 = 0.1), seed = 55)
  apollo_inputs <- fx$apollo_inputs

  # Act
  cv <- mtuem_ab_likelihood(fx$mtuem_settings, "get_covar")

  # Assert
  err <- fx$obs - fx$opt
  expect_equal(unname(cv$covar), unname(stats::cov(err)), tolerance = 1e-10)
  expect_equal(unname(cv$corr),  unname(stats::cor(err)), tolerance = 1e-10)
  expect_equal(unname(cv$sigma), unname(sqrt(diag(stats::cov(err)))), tolerance = 1e-10)
})

test_that("alpha-beta handles the Apollo bookkeeping functionalities", {
  # Arrange
  fx <- ab_fixture(n = 13, seed = 56)
  apollo_inputs <- fx$apollo_inputs
  named <- fx$mtuem_settings; named$componentName <- "ab_component"

  # Act / Assert
  expect_equal(mtuem_ab_likelihood(fx$mtuem_settings, "validate"), rep(1, 13))
  expect_equal(mtuem_ab_likelihood(fx$mtuem_settings, "zero_LL"), rep(NA, 13))
  expect_equal(mtuem_ab_likelihood(named, "preprocess")$componentName, "ab_component")
  expect_equal(mtuem_ab_likelihood(fx$mtuem_settings, "gradient"), list(like = NA, grad = NA))
  expect_error(mtuem_ab_likelihood(fx$mtuem_settings, "nonsense"), "Invalid value of argument")
})

test_that("alpha-beta prediction closes both budgets and reports values of time", {
  # Arrange
  fx <- ab_fixture(n = 25, times_elasticities = list(th1 = 0.1),
                   goods_elasticities = list(ph1 = 0.2), seed = 57)
  apollo_inputs <- fx$apollo_inputs

  # Act
  pred <- mtuem_ab_likelihood(fx$mtuem_settings, "prediction")

  # Assert
  expect_equal(colnames(pred), c("Tw", "Tf1", "Tfi", "Ef1", "Xfj", "VoL", "VTAW"))
  expect_equal(pred[, "Tw"] + pred[, "Tf1"] + pred[, "Tfi"], fx$tau - fx$Tc, tolerance = 1e-9)
  expect_equal(pred[, "Ef1"] + pred[, "Xfj"], fx$w * pred[, "Tw"] - fx$Ec, tolerance = 1e-9)
  expect_equal(pred[, "VoL"], fx$w + pred[, "VTAW"], tolerance = 1e-9)
})

test_that("alpha-beta output matches estimate and report prints the values of time", {
  # Arrange
  fx <- ab_fixture(n = 20, seed = 58)
  apollo_inputs <- fx$apollo_inputs
  settings <- fx$mtuem_settings
  settings$sig <- list(2, 2)

  # Act / Assert
  expect_equal(mtuem_ab_likelihood(settings, "output"),
               mtuem_ab_likelihood(settings, "estimate"))
  expect_output(mtuem_ab_likelihood(settings, "report"), "VTAW")
})
