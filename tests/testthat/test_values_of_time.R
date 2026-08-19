test_that("single class values of time match the closed form", {
  set.seed(7)
  n <- 50
  db <- data.frame(Tc = runif(n, 40, 80), EcI = runif(n, 50, 150), w = runif(n, 10, 25))
  we <- list(Theta = 1, Phi = 1.4, thw = 0.2)
  model <- list(estimate = c(Theta = we$Theta, Phi = we$Phi, thw = we$thw))

  res <- mtuem_values_of_time(model, list(database = db), Ec = "EcI",
                              standard_errors = FALSE)

  Tw <- as.numeric(get_tw_thph(we, 168, db$Tc, db$EcI, db$w))
  expected_VoL  <- mean((db$w * Tw - db$EcI) / (168 - Tw - db$Tc)) * we$Theta / we$Phi
  expected_VTAW <- mean((db$w * Tw - db$EcI) / Tw) * we$thw / we$Phi

  expect_equal(res$summary$mean[res$summary$quantity == "VoL"], expected_VoL)
  expect_equal(res$summary$mean[res$summary$quantity == "VTAW"], expected_VTAW)
  expect_equal(res$summary$mean[res$summary$quantity == "wage"], mean(db$w))
  expect_equal(nrow(res$individual), n)
})

test_that("VoL minus VTAW is the wage rate per observation", {
  set.seed(8)
  n <- 40
  db <- data.frame(Tc = runif(n, 40, 80), EcI = runif(n, 50, 150), w = runif(n, 10, 25))
  model <- list(estimate = c(Theta = 1, Phi = 1.1, thw = 0.35))

  res <- mtuem_values_of_time(model, list(database = db), Ec = "EcI",
                              standard_errors = FALSE)

  expect_equal(res$individual$VoL - res$individual$VTAW, db$w)
})

test_that("class weights change the reported means", {
  set.seed(9)
  n <- 60
  db <- data.frame(Tc = runif(n, 40, 80), EcI = runif(n, 50, 150), w = runif(n, 10, 25))
  model <- list(estimate = c(Theta_1 = 1, Phi_1 = 1.4, thw_1 = 0.2,
                             Theta_2 = 1, Phi_2 = 1.4, thw_2 = 0.2))

  # Identical parameters in both classes, so only the weights can differ:
  # class 1 keeps the first half of the sample, class 2 the second half.
  W <- cbind(c(rep(1, n / 2), rep(0, n / 2)), c(rep(0, n / 2), rep(1, n / 2)))
  res <- mtuem_values_of_time(model, list(database = db), nClass = 2, weights = W,
                              Ec = "EcI", standard_errors = FALSE)

  expect_equal(res$summary$mean[res$summary$class == 1 & res$summary$quantity == "wage"],
               mean(db$w[1:(n / 2)]))
  expect_equal(res$summary$mean[res$summary$class == 2 & res$summary$quantity == "wage"],
               mean(db$w[(n / 2 + 1):n]))
  expect_false(isTRUE(all.equal(
    res$summary$mean[res$summary$class == 1 & res$summary$quantity == "VoL"],
    res$summary$mean[res$summary$class == 2 & res$summary$quantity == "VoL"])))
})

test_that("equal weights reproduce the unweighted mean", {
  set.seed(10)
  n <- 30
  db <- data.frame(Tc = runif(n, 40, 80), EcI = runif(n, 50, 150), w = runif(n, 10, 25))
  model <- list(estimate = c(Theta = 1, Phi = 1.2, thw = 0.1))

  weighted <- mtuem_values_of_time(model, list(database = db), nClass = 1,
                                   weights = rep(0.5, n), Ec = "EcI",
                                   standard_errors = FALSE)
  plain <- mtuem_values_of_time(model, list(database = db), nClass = 1, Ec = "EcI",
                                standard_errors = FALSE)

  expect_equal(weighted$summary$mean, plain$summary$mean)
})

test_that("a missing class parameter is reported by name", {
  db <- data.frame(Tc = 60, EcI = 100, w = 15)
  model <- list(estimate = c(Theta = 1, Phi = 1, thw = 0))

  expect_error(
    mtuem_values_of_time(model, list(database = db), nClass = 2, weights = matrix(1, 1, 2),
                         Ec = "EcI", standard_errors = FALSE),
    "Theta_1")
})

test_that("weights of the wrong shape are rejected", {
  db <- data.frame(Tc = c(60, 62), EcI = c(100, 110), w = c(15, 16))
  model <- list(estimate = c(Theta_1 = 1, Phi_1 = 1, thw_1 = 0,
                             Theta_2 = 1, Phi_2 = 1, thw_2 = 0))

  expect_error(
    mtuem_values_of_time(model, list(database = db), nClass = 2,
                         weights = matrix(1, nrow = 3, ncol = 2), Ec = "EcI",
                         standard_errors = FALSE),
    "must have 2 rows")
})
