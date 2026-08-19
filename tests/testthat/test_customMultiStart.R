library(testthat)
library(apollo)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

# Minimal single-class, single-equation model on a small maed subset. Returns
# everything customMultiStart needs. apollo_validateInputs() reads
# apollo_control, database, apollo_beta and apollo_fixed from the global
# environment, so they are mirrored there and removed by the caller.
cms_setup <- function(env) {
  apollo_initialise()
  apollo_control <- list(modelName = "test-cms", indivID = "PeID",
                         outputDirectory = tempdir(),
                         noValidation = TRUE, noDiagnostics = TRUE)

  database <- mtuem::maed
  database <- database[database$EcI > 0, ][1:120, ]

  apollo_beta <- c(Theta = 1, Phi = 1, thw = 0, sig1 = 8)
  apollo_fixed <- c("Theta")

  globals <- list(apollo_control = apollo_control, database = database,
                  apollo_beta = apollo_beta, apollo_fixed = apollo_fixed)
  for (nm in names(globals)) assign(nm, globals[[nm]], envir = .GlobalEnv)
  assign("cms_globals", names(globals), envir = env)

  apollo_inputs <- suppressMessages(apollo_validateInputs())

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    mtuem_settings <- list(
      work_times = "Tw",
      work_elasticities = list(Theta = Theta, Phi = Phi, thw = thw),
      sig = list(sig1 = sig1),
      Tc = Tc, Ec = EcI, w = w, tau = 168)
    P[["model"]] <- mtuem_likelihood(mtuem_settings, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }

  list(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
       apollo_inputs = apollo_inputs, apollo_probabilities = apollo_probabilities)
}

cms_teardown <- function(env) {
  nms <- get("cms_globals", envir = env)
  suppressWarnings(rm(list = nms, envir = .GlobalEnv))
}

fast_estimate_settings <- list(maxIterations = 2, hessianRoutine = "none",
                               estimationRoutine = "bfgs")

test_that("customMultiStart returns the best model, its log-likelihood and every candidate", {
  # Arrange
  env <- environment()
  s <- cms_setup(env)
  on.exit(cms_teardown(env), add = TRUE)

  # Act
  out <- suppressWarnings(utils::capture.output(
    res <- customMultiStart(s$apollo_beta, s$apollo_fixed, s$apollo_probabilities,
                            s$apollo_inputs,
                            customMultistart_settings = list(nCandidates = 3),
                            estimate_settings = fast_estimate_settings,
                            verbose = FALSE)))

  # Assert
  expect_named(res, c("best_model", "best_ll", "models"))
  expect_length(res$models, 3)
  succeeded <- Filter(function(m) inherits(m, "apollo"), res$models)
  expect_gt(length(succeeded), 0)
  if (!identical(res$best_model, NA)) {
    expect_equal(res$best_ll, res$best_model$maximum)
    expect_true(res$best_ll >= max(vapply(succeeded, function(m) m$maximum, numeric(1))) - 1e-8)
  }
})

test_that("customMultiStart samples candidates inside the requested box", {
  # Arrange: pinning the box to a single point makes every candidate identical,
  # so all of them must reach the same log-likelihood.
  env <- environment()
  s <- cms_setup(env)
  on.exit(cms_teardown(env), add = TRUE)
  pinned <- s$apollo_beta

  # Act
  out <- suppressWarnings(utils::capture.output(
    res <- customMultiStart(s$apollo_beta, s$apollo_fixed, s$apollo_probabilities,
                            s$apollo_inputs,
                            customMultistart_settings = list(nCandidates = 2,
                                                             apolloBetaMin = pinned,
                                                             apolloBetaMax = pinned),
                            estimate_settings = fast_estimate_settings,
                            verbose = FALSE)))

  # Assert
  lls <- vapply(Filter(function(m) inherits(m, "apollo"), res$models),
                function(m) m$maximum, numeric(1))
  expect_gt(length(lls), 0)
  expect_equal(diff(range(lls)), 0, tolerance = 1e-8)
})

test_that("theta_phi normalization keeps every sampled candidate feasible", {
  # Arrange: a sampling box wide enough to produce negative Phi and negative
  # sigma. Without normalization those candidates give a non-finite
  # log-likelihood and get discarded; with it, all of them must survive the
  # initial screening.
  env <- environment()
  s <- cms_setup(env)
  on.exit(cms_teardown(env), add = TRUE)
  wide_min <- s$apollo_beta - 12
  wide_max <- s$apollo_beta + 12

  # Act
  printed <- suppressWarnings(utils::capture.output(
    res <- customMultiStart(s$apollo_beta, s$apollo_fixed, s$apollo_probabilities,
                            s$apollo_inputs,
                            customMultistart_settings = list(nCandidates = 5,
                                                             apolloBetaMin = wide_min,
                                                             apolloBetaMax = wide_max),
                            estimate_settings = fast_estimate_settings,
                            verbose = FALSE,
                            normalization = list(
                              normalization = "theta_phi",
                              work_elasticities = c("Theta", "Phi", "thw"),
                              sig = c("sig1")))))

  # Assert
  expect_true(any(grepl("5 candidates available", printed)))
})

test_that("the deprecated estimation_settings alias is still honoured", {
  # Arrange
  env <- environment()
  s <- cms_setup(env)
  on.exit(cms_teardown(env), add = TRUE)

  # Act
  result <- tryCatch({
    suppressWarnings(utils::capture.output(
      res <- customMultiStart(s$apollo_beta, s$apollo_fixed, s$apollo_probabilities,
                              s$apollo_inputs,
                              customMultistart_settings = list(nCandidates = 2),
                              estimation_settings = fast_estimate_settings,
                              verbose = FALSE)))
    "ok"
  }, error = function(e) e)

  # Assert
  expect_identical(result, "ok")
})
