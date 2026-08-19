library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("all documented datasets ship with the package", {
  # Arrange / Act
  items <- utils::data(package = "mtuem")$results[, "Item"]

  # Assert
  expect_true(all(c("maed", "enut.i", "enut.i.raw", "enut.ii", "enut.ii.raw") %in% items))
})

test_that("maed carries the columns the model examples rely on", {
  # Arrange
  d <- mtuem::maed

  # Act / Assert
  expect_true(is.data.frame(d))
  expect_gt(nrow(d), 0)
  expect_true(all(c("PeID", "Tw", "Tc", "Tf1", "Tf2", "EcI", "w",
                    "Ef1", "Ef2", "ta",
                    "female", "older45", "professional", "self_employed") %in% names(d)))
  expect_equal(anyDuplicated(d$PeID), 0L)
})

test_that("maed has no missing values", {
  # Arrange
  d <- mtuem::maed

  # Act / Assert
  expect_equal(sum(is.na(d)), 0L)
})

test_that("maed satisfies the sign restrictions the likelihood assumes", {
  # Arrange
  d <- mtuem::maed

  # Act / Assert
  expect_true(all(d$w > 0))
  expect_true(all(d$Tw > 0))
  expect_true(all(d$Tc > 0))
  expect_true(all(d$Tf1 >= 0))
  expect_true(all(d$Tf2 >= 0))
  expect_true(all(d$ta == 168))
  expect_true(all(d$Tw + d$Tc < d$ta))
})

test_that("the maed estimation subset has strictly positive committed expenditure", {
  # Arrange: the examples filter on EcI > 0 because a non-positive committed
  # expenditure breaks the work equation, so the subset must be non-empty.
  d <- mtuem::maed

  # Act
  keep <- d[d$EcI > 0, ]

  # Assert
  expect_gt(nrow(keep), 500)

  # Disposable income w*Tw - EcI is negative for a minority of the sample
  # (42 of 712 rows at the time of writing). Those rows give a negative
  # money shadow price and therefore a negative value of leisure, so this
  # pins the share rather than asserting it away.
  disposable <- keep$w * keep$Tw - keep$EcI
  expect_lt(mean(disposable <= 0), 0.1)
})

test_that("maed time use exhausts the time budget, so the 3-equation time model is unidentified", {
  # Arrange: Tw + Tf1 + Tf2 + Tc equals the time budget up to rounding, which
  # means those four categories leave nothing out. Estimating work plus both
  # free times together therefore violates the leave-one-out requirement and
  # the likelihood is unbounded. This test pins the data property that makes
  # that so, so the constraint is not silently broken by a data update.
  d <- mtuem::maed

  # Act
  slack <- d$ta - (d$Tw + d$Tf1 + d$Tf2 + d$Tc)

  # Assert
  expect_true(all(abs(slack) <= 1.05))
})

test_that("the maed binary covariates really are binary", {
  # Arrange
  d <- mtuem::maed

  # Act / Assert
  for (v in c("female", "older45", "professional", "self_employed")) {
    expect_true(all(as.numeric(d[[v]]) %in% c(0, 1)), info = v)
  }
})
