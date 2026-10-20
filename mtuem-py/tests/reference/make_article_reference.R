# Reference fixtures for the activity participation article (latent class Cholesky hobbies model).
# Run from the root of the article repository (mtuem-lc-activity-participation), which is only read:
#   Rscript <mtuem-py>/tests/reference/make_article_reference.R <run> <penalization> [<output folder>]
# <run> is 2015, 2015-new, 2022, 2022-bigcities (conce_valpo location set). The output
# folder defaults to tests/reference/article next to this script. It writes <run>_p<k>.parquet (the
# model database restricted to the columns the likelihood uses plus the id, in R's row order) and
# <run>_p<k>.json (model definition, start values, R estimates read from the article's saved fits)
# there and nothing else. It also checks that the article's apollo_probabilities at the saved
# estimate, evaluated on the exported parquet, reproduces the saved maximum.
suppressMessages({library(tidyverse); library(apollo); library(mtuem); library(jsonlite); library(arrow)})
args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: make_article_reference.R <run> <penalization> [<output folder>]")
run <- args[1]
level <- as.numeric(args[2])
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
outdir <- if (length(args) >= 3) args[3] else file.path(dirname(normalizePath(script_file)), "article")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists("R/models/hobbies_cholesky.R")) stop("Run from the root of the article repository")

set.seed(42)
study_year <- run
# R/config/project_setup.R without ensure_project_directories(), which would create folders in the
# article repository.
for (expression in parse("R/config/project_setup.R")) {
  if (is.call(expression) && identical(expression[[1]], as.name("ensure_project_directories"))) next
  eval(expression, globalenv())
}
source("R/config/hobbies_config.R")
source("R/data/enut_data.R")
source("R/data/hobbies_data.R")
source("R/models/initial_values.R")
source("R/models/hobbies_cholesky.R")
source("R/models/penalization_grid.R")
apollo_initialise()

prepared_data <- prepare_hobbies_data(config = hobbies_config())
database <- prepared_data$database
nClass <- prepared_data$n_classes
socioecon <- prepared_data$socioeconomic_variables
especificas <- prepared_data$class_covariates
model_directory <- project_paths$penalization_models
if (identical(run, "2022-bigcities")) {
  # As scripts/2022-bigcities/experiments/fit_penalization_grid.R with the default --location.
  location <- c("conce", "valpo")
  socioecon <- c(setdiff(prepared_data$socioeconomic_variables, c("Metropolitana", "conce", "valpo")), location)
  especificas <- c("asc", paste0("x_", socioecon))
  model_directory <- file.path(dirname(project_paths$penalization_models), paste(location, collapse = "_"))
}
initialise_hobbies_cholesky()

model_name <- paste0("corner_hobbies_p", level)
path <- function(suffix) file.path(model_directory, paste0(model_name, suffix))
model <- readRDS(path("_model.rds"))
estimates_csv <- read.csv(path("_estimates.csv"), row.names = 1, check.names = FALSE)
convergence <- if (file.exists(path("_convergence.csv"))) read.csv(path("_convergence.csv")) else NULL
output_lines <- readLines(path("_output.txt"), warn = FALSE)
output_nobs <- as.integer(sub(".*:\\s*", "", grep("^Number of rows in database", output_lines, value = TRUE)[1]))

beta_start <- get_apollo_beta()
fixed <- get_apollo_fixed()
indiv_id <- model$apollo_control$indivID
if (!setequal(names(beta_start), names(model$estimate))) stop("get_apollo_beta() and the saved estimate differ in names")
if (!setequal(fixed, model$apollo_fixed)) stop("get_apollo_fixed() and the saved fixed parameters differ")
if (!isTRUE(all.equal(unname(estimates_csv[names(model$estimate), "Estimate"]), unname(model$estimate), tolerance = 1e-12))) {
  stop("estimates.csv and model.rds disagree")
}

# Floor and warm start of the kept attempt. fit_penalization_grid() records the floor in the
# convergence CSV; a warm start fits by apollo_estimate (no EM) and is labelled "warm start from p<k>"
# in its attempts column. The 2015 CSVs predate that column, so the Apollo routine string is the
# evidence there: EM fits start from get_apollo_beta(), warm starts never run EM.
class_likelihood_floor <- if (!is.null(convergence)) convergence$class_likelihood_floor[1] else 0
attempts <- if (!is.null(convergence) && !is.null(convergence$attempts)) convergence$attempts[1] else NA_character_
em_fit <- grepl("^EM algorithm", model$estimationRoutine)
warm_started <- !em_fit || (!is.na(attempts) && grepl("warm start", attempts))
warm_evidence <- paste0("estimationRoutine '", model$estimationRoutine, "'",
                        if (is.na(attempts)) "; convergence CSV has no attempts column" else paste0("; attempts '", attempts, "'"))
# An EM fit starts from get_apollo_beta() (fit_hobbies_cholesky() with start = NULL); Apollo saves
# the post EM values as apollo_beta, so only a warm start's saved start is its own.
beta_start_used <- if (warm_started) model$apollo_beta else beta_start

# Estimation settings: the saved ones when Apollo stored them (bgw fits), otherwise est_set of
# R/models/hobbies_cholesky.R with the routine of the saved fit. apollo_lcEM ran with its defaults.
est_used <- if (!is.null(model$estimate_settings)) model$estimate_settings[names(est_set)] else est_set
routine <- sub("^.*Maximum likelihood \\((.*)\\)$", "\\1", model$estimationRoutine)
if (!identical(routine, est_used$estimationRoutine)) stop("Routine mismatch: ", routine, " vs ", est_used$estimationRoutine)
lcEM_settings <- list(EMstoppingCriterion = 1e-5, EMmaxIterations = 100, postEM = 2, silent = FALSE)

# Columns of the likelihood: work and free times, Tc, Ec, w (attached database columns), the
# allocation covariates, and the id.
columns <- c(indiv_id, "Tw", times_internal, "Tc", "Ec", "w", socioecon)
exported <- database[, columns]
for (v in setdiff(columns, indiv_id)) exported[[v]] <- as.numeric(exported[[v]])
if (is.numeric(exported[[indiv_id]])) exported[[indiv_id]] <- as.numeric(exported[[indiv_id]])
if (nrow(exported) != model$nObs) stop("Row count ", nrow(exported), " differs from the saved nObs ", model$nObs)
tag <- paste0(run, "_p", level)
parquet_file <- file.path(outdir, paste0(tag, ".parquet"))
write_parquet(exported, parquet_file)

# Self-consistency: the article's apollo_probabilities at the saved estimate on the exported data.
penalization <- level
apollo_control <- list(modelName = model_name, modelDescr = model_name, indivID = indiv_id,
                       outputDirectory = tempdir(), noValidation = TRUE, noDiagnostics = TRUE)
apollo_beta <- model$estimate
apollo_fixed <- fixed
apollo_inputs <- apollo_validateInputs(apollo_control = apollo_control, apollo_beta = apollo_beta,
                                       apollo_fixed = apollo_fixed, database = as.data.frame(read_parquet(parquet_file)),
                                       silent = TRUE)
P <- apollo_probabilities(model$estimate, apollo_inputs, functionality = "estimate")
if (is.list(P)) P <- P[["model"]]
ll_check <- sum(log(P))
rel <- abs(ll_check - model$maximum) / abs(model$maximum)
cat(sprintf("%s: LL at r_estimate %.10f, saved maximum %.10f, rel %.3e\n", tag, ll_check, model$maximum, rel))
if (rel > 1e-8) stop("LL at the saved estimate does not reproduce the saved maximum")

classes <- seq_along(act_mix)
cholesky_params <- lapply(classes, function(s) {
  active <- c(TRUE, act_mix[[s]])
  unlist(lapply(seq_len(nequations), function(i) unlist(lapply(seq_len(i), function(j) {
    if (active[i] && active[j]) paste0("chk", i, j, "_", s)
  }))))
})
out <- list(
  run = run, penalization = level, model_name = model_name,
  likelihood_floor = class_likelihood_floor, class_likelihood_floor = class_likelihood_floor,
  warm_started = warm_started, warm_started_evidence = warm_evidence,
  indiv_id = indiv_id, n_classes = length(classes),
  work_times = "Tw", times_internal = times_internal, free_times = times_internal,
  thetas = thetas, Tc = "Tc", Ec = "Ec", w = "w", tau = 168, optimal_tw = TRUE, nequations = nequations,
  act_mix = lapply(act_mix, function(m) setNames(as.list(m), times_internal)),
  cholesky_params = cholesky_params,
  socioecon = socioecon, allocation = socioecon,
  allocation_params = c("asc", paste0("x_", tolower(socioecon))),
  penalty = "class likelihood times exp(-penalization * sum of the hours of times_internal the class does not practice), then pmax(., class_likelihood_floor) when the floor is positive",
  beta_start = as.list(beta_start), beta_start_used = as.list(beta_start_used), fixed = fixed,
  est_set = est_used, lcEM_settings = lcEM_settings, estimation_routine = model$estimationRoutine,
  constraints = list(),
  r_estimate = as.list(model$estimate), r_se = as.list(model$se[names(model$estimate)]),
  r_robse = as.list(model$robse[names(model$estimate)]),
  r_maximum = model$maximum, r_loglik = model$maximum, r_ll_check = ll_check,
  r_nobs = model$nObs, r_nobs_output = output_nobs,
  r_minutes = if (!is.null(convergence)) convergence$minutes[1] else as.numeric(model$timeTaken) / 60,
  r_converged = if (!is.null(convergence)) convergence$converged[1] else NA,
  r_message = model$message
)
write_json(out, file.path(outdir, paste0(tag, ".json")), digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("Wrote", parquet_file, "\n")
