#' Compute residuals from optimal allocations for a MTUEM specification (Theta-Phi parameterization)
#'
#' Computes the optimal work time, free time, and free goods allocations implied by
#' \code{mtuem_settings} and returns the difference between the observed allocations
#' and these optimal allocations (observed minus optimal). This is the same residual
#' computation used internally by \link{mtuem_likelihood} for the "estimate",
#' "conditionals", and "raw" functionalities, factored out so it can be reused by
#' other likelihood functions (e.g. latent class specifications that concentrate the
#' error covariance out of the likelihood).
#'
#' @param mtuem_settings List of arguments describing the MTUEM specification. See \link{mtuem_likelihood}
#'                      for the full list of fields. Only the fields relevant to computing optimal
#'                      allocations are used: \code{work_times}, \code{free_times}, \code{free_goods},
#'                      \code{goods_cost}, \code{work_elasticities}, \code{times_elasticities},
#'                      \code{goods_elasticities}, \code{optimal_tw}, \code{Tc}, \code{Ec}, \code{w}, and \code{tau}.
#'                      \code{sig}, \code{rho}, and \code{cholesky} are not used and may be omitted.
#' @param apollo_inputs List grouping most common inputs, as created by \link{apollo_validateInputs}. Must
#'                      contain either \code{obs_matrix} (a matrix of observed allocations) or a \code{database}
#'                      with columns matching \code{work_times}, \code{free_times}, and \code{free_goods}.
#' @return N x M numeric matrix of residuals (observed minus optimal), with columns
#'         \code{c(work_times, free_times, free_goods)}.
#' @export
mtuem_residuals <- function(mtuem_settings, apollo_inputs) {
  default = list(
    free_times = list(),
    free_goods = list(),
    goods_cost = 1,
    times_elasticities = list(),
    goods_elasticities = list(),
    optimal_tw = T
  )

  tmp <- names(default)[!(names(default) %in% names(mtuem_settings))]
  for (i in tmp) mtuem_settings[[i]] <- default[[i]]
  rm(tmp)

  work_times <- mtuem_settings$work_times
  free_times <- unlist(mtuem_settings$free_times)
  free_goods <- unlist(mtuem_settings$free_goods)
  goods_cost <- mtuem_settings$goods_cost
  work_elasticities <- mtuem_settings$work_elasticities
  times_elasticities <- mtuem_settings$times_elasticities
  goods_elasticities <- mtuem_settings$goods_elasticities
  optimal_tw <- mtuem_settings$optimal_tw
  Tc  <- mtuem_settings$Tc
  Ec  <- mtuem_settings$Ec
  w   <- mtuem_settings$w
  tau <- mtuem_settings$tau

  tw_opt <- get_tw_thph(work_elasticities, tau, Tc, Ec, w)
  if (optimal_tw) {
    ti_opt <- get_ti_thph(times_elasticities, work_elasticities$Theta, tw_opt, tau, Tc)
    xj_opt <- get_xi_thph(goods_elasticities, goods_cost, work_elasticities$Phi, tw_opt, Ec, w)
  } else {
    ti_opt <- get_ti_thph(times_elasticities, work_elasticities$Theta, apollo_inputs$database[, work_times], tau, Tc)
    xj_opt <- get_xi_thph(goods_elasticities, goods_cost, work_elasticities$Phi, apollo_inputs$database[, work_times], Ec, w)
  }

  opt <- cbind(tw_opt, unlist(ti_opt), unlist(xj_opt))
  colnames(opt) <- c(work_times, free_times, free_goods)

  if (is.null(apollo_inputs$obs_matrix)) {
    obs <- as.matrix(apollo_inputs$database[, colnames(opt)])
  } else {
    obs <- apollo_inputs$obs_matrix
  }

  return(obs - opt)
}
