# Polished R optimum of one article target, added to its reference fixture.
# Run from the root of the article repository (mtuem-lc-activity-participation), which is only read,
# after make_article_reference.R has written the fixture:
#   Rscript <mtuem-py>/tests/reference/make_article_polish.R <run> <penalization>
# The article's fits stop before the optimum: Apollo's bfgs is maxLik over optim BFGS, whose default
# reltol (about 1.5e-8) lets it stop once an iteration gains about 1e-3 on an LL near 6.6e4. This
# script rebuilds the target's model as make_article_reference.R does, on the fixture parquet, and
# runs apollo_estimate from r_estimate with the fixture est_set plus maxLik_settings reltol = 1e-15.
# Other settings: estimationRoutine is bfgs for every target (2015_p5 was fitted with bgw, which
# ignores maxLik_settings), and maxIterations and iterlim are 2000 so the iteration limit is not
# the stopping rule. It then checks that max |gradient| at the polished point (numDeriv::grad,
# Richardson, over the free parameters) is below 1e-2 and adds r_polished_estimate, r_polished_se,
# r_polished_robse, r_polished_loglik, r_polished_max_grad and r_polished_settings to
# tests/reference/article/<run>_p<k>.json. Apollo's output goes to a temporary folder.
suppressMessages({library(tidyverse); library(apollo); library(mtuem); library(jsonlite); library(arrow); library(numDeriv)})
args <- commandArgs(TRUE)
if (length(args) != 2) stop("usage: make_article_polish.R <run> <penalization>")
run <- args[1]
level <- as.numeric(args[2])
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
article_dir <- file.path(dirname(normalizePath(script_file)), "article")
if (!file.exists("R/models/hobbies_cholesky.R")) stop("Run from the root of the article repository")
max_grad_allowed <- 1e-2

set.seed(42)
study_year <- run
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
nClass <- prepared_data$n_classes
socioecon <- prepared_data$socioeconomic_variables
especificas <- prepared_data$class_covariates
if (identical(run, "2022-bigcities")) {
  location <- c("conce", "valpo")
  socioecon <- c(setdiff(prepared_data$socioeconomic_variables, c("Metropolitana", "conce", "valpo")), location)
  especificas <- c("asc", paste0("x_", socioecon))
}
initialise_hobbies_cholesky()

tag <- paste0(run, "_p", level)
json_file <- file.path(article_dir, paste0(tag, ".json"))
# These simplify options and the write_json options below rewrite the existing fields byte for byte.
fixture <- fromJSON(json_file, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
r_est <- unlist(fixture$r_estimate)
fixed <- unlist(fixture$fixed)
if (!setequal(fixed, get_apollo_fixed())) stop("fixed parameters differ from get_apollo_fixed()")
free <- setdiff(names(r_est), fixed)

penalization <- level
apollo_control <- list(modelName = paste0(tag, "_polish"), modelDescr = tag, indivID = fixture$indiv_id,
                       outputDirectory = tempdir(), noValidation = TRUE, noDiagnostics = TRUE)
apollo_beta <- r_est
apollo_fixed <- fixed
database <- as.data.frame(read_parquet(file.path(article_dir, paste0(tag, ".parquet"))))
apollo_inputs <- apollo_validateInputs(apollo_control = apollo_control, apollo_beta = apollo_beta,
                                       apollo_fixed = apollo_fixed, database = database, silent = TRUE)
ll_at <- function(b) {
  P <- apollo_probabilities(b, apollo_inputs, functionality = "estimate")
  if (is.list(P)) P <- P[["model"]]
  sum(log(P))
}
ll_r <- ll_at(r_est)
if (abs(ll_r - fixture$r_maximum) / abs(fixture$r_maximum) > 1e-8) stop("LL at r_estimate does not reproduce r_maximum")

settings <- lapply(fixture$est_set, identity)
settings$estimationRoutine <- "bfgs"
settings$maxIterations <- 2000
settings$maxLik_settings <- list(reltol = 1e-15, iterlim = 2000)
t0 <- Sys.time()
model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, estimate_settings = settings)
minutes <- as.numeric(Sys.time() - t0, units = "mins")
est <- model$estimate[names(r_est)]
ll_end <- ll_at(est)
if (abs(ll_end - model$maximum) > 1e-6) stop("LL at the polished estimate differs from Apollo's maximum")
ll_free <- function(x) { b <- est; b[free] <- x; ll_at(b) }
grad <- numDeriv::grad(ll_free, est[free], method = "Richardson")
max_grad <- max(abs(grad))
cat(sprintf("%s: LL %.8f -> %.8f, max|grad| %.3g, %d iterations, %.1f min, '%s'\n",
            tag, ll_r, ll_end, max_grad, as.integer(model$nIter), minutes, model$message))
if (max_grad > max_grad_allowed) stop(sprintf("max |gradient| %.3g at the polished point exceeds %g", max_grad, max_grad_allowed))
if (ll_end < ll_r) stop("polishing lowered the LL")

na_null <- function(x) lapply(x, function(v) if (is.na(v)) NULL else v)
fixture$r_polished_estimate <- as.list(est)
fixture$r_polished_se <- na_null(model$se[names(r_est)])
fixture$r_polished_robse <- na_null(model$robse[names(r_est)])
fixture$r_polished_loglik <- model$maximum
fixture$r_polished_max_grad <- max_grad
fixture$r_polished_settings <- list(
  script = "tests/reference/make_article_polish.R", start = "r_estimate",
  estimate_settings = settings[c("estimationRoutine", "maxIterations", "maxLik_settings", "hessianRoutine",
                                 "scaleHessian", "scaleAfterConvergence")],
  gradient = "numDeriv::grad, method Richardson, over the free parameters",
  max_grad_allowed = max_grad_allowed, message = model$message, iterations = as.integer(model$nIter),
  minutes = minutes, ll_start = ll_r, apollo_version = as.character(packageVersion("apollo")),
  r_version = R.version.string
)
write_json(fixture, json_file, digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
cat("Wrote", json_file, "\n")
