#' Function for getting the correlation between errors in the model
#'
#' @param model Object. Model returned by apollo estimation routines.
#' @param apollo_inputs List grouping most common inputs. Created by function apollo_validateInputs.
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated.
#' @param vars Character vector containing the name of the time-use-expenditure allocations in the database
#' @return array of likelihood for each individual
#' @export
mtuem_get_corr <- function(model, apollo_probabilities, apollo_inputs, vars) {
  pred <- apollo_prediction(model, apollo_probabilities, apollo_inputs)
  database <- apollo_inputs$database

  return(list(
    covar = stats::cov(database[, vars]-pred[,vars], use = "pairwise.complete.obs"),
    corr  = stats::cor(database[, vars]-pred[,vars], use = "pairwise.complete.obs")
  ))
}