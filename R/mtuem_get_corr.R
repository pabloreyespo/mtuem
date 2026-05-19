#' Function for getting the correlation between errors in the model
#'
#' @param model Object. Model returned by apollo estimation routines.
#' @param apollo_inputs List grouping most common inputs. Created by function apollo_validateInputs.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated.
#' @param vars Character vector containing the name of the time-use-expenditure allocations in the database
#' @return array of likelihood for each individual
#' @export

# mtuem_get_corr <- function(model, apollo_probabilities, apollo_inputs, vars) {
#   pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs, prediction_settings = list(silent = T))
#   database <- apollo_inputs$database
#
#   return(list(
#     covar = stats::cov(database[, vars]-pred[,vars], use = "pairwise.complete.obs"),
#     corr  = stats::cor(database[, vars]-pred[,vars], use = "pairwise.complete.obs")
#   ))
# }

mtuem_get_corr <- function(model, apollo_probabilities, apollo_inputs, vars) {
  pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs,
                            prediction_settings = list(silent = TRUE))
  database <- apollo_inputs$database

  err <- database[, vars] - pred[, vars]

  # Extract estimated sig from model in the same order as vars
  sig_est <- model$estimate[paste0("sig_", seq_along(vars))]

  # Standardize residuals the same way the likelihood does
  err_std <- sweep(err, MARGIN = 2, sig_est, "/")

  return(list(
    covar     = stats::cov(err,     use = "pairwise.complete.obs"),
    corr      = stats::cor(err,     use = "pairwise.complete.obs"),
    corr_std  = stats::cor(err_std, use = "pairwise.complete.obs")  # comparable to rho
  ))
}