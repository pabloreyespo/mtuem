#' Likelihood function for the parameter estimation of a Microeconomic Time-Use-Expenditure Model (Alpha-Beta parameterization).
#'
#' The implementation is based on Jara-Diaz et al. (2008) consisting of the maximization of a Cobb-Douglas utility function dependent on Time-use (T) and Goods Consumption (X).
#' The model is subject to time budget, monetary budget and technical time-use-expenditure constraints (committed time and expenditures).
#'
#' @details
#' This function uses the Alpha-Beta parameterization of the MTUEM, an alternative to the
#' Theta-Phi formulation in \code{\link{mtuem_likelihood}}. The work equation is specified
#' through \code{work_elasticities} with mandatory names \code{alpha} and \code{beta}.
#'
#' Free time and free goods allocations are derived from the work equation using the corresponding
#' elasticity parameters. For identification, at least one free time and one free good must be
#' left out of the model (the "leave one out" constraint). The omitted categories are recovered
#' as residuals from the time and monetary budgets.
#'
#' The error structure can be specified in three ways:
#' \itemize{
#'   \item No \code{sig} and no \code{rho}: The covariance matrix is inferred from residuals.
#'     This is the simplest approach and is only valid for single-equation models.
#'   \item \code{sig} specified but no \code{rho}: Standard deviations are estimated as parameters,
#'     and correlations are inferred from residuals. Valid for multi-equation models.
#'   \item Both \code{sig} and \code{rho} specified: Full covariance estimation where both
#'     standard deviations and correlations are estimated as parameters.
#' }
#'
#' When \code{functionality = "get_covar"}, the function returns a list with the covariance matrix,
#' correlation matrix, and standard deviations computed from the residuals.
#'
#' @param mtuem_settings List of arguments to the functions. It must contain the following. User input is required for all settings except those with a default or marked as optional.
#'                      \itemize{
#'                        \item \strong{\code{work_times}} Character. Name of the work time column in the database.
#'                        \item \strong{\code{free_times}} Character vector. Optional. Names of the freely allocated time columns. Leave at least one out for identification.
#'                        \item \strong{\code{free_goods}} Character vector. Optional. Names of the freely consumed goods columns. Leave at least one out for identification.
#'                        \item \strong{\code{goods_cost}} Named list. Optional. Cost associated with the consumption of each good. If not specified, assumed as 1 (expenses equal consumption of the good).
#'                        \item \strong{\code{work_elasticities}} Named list. Value of the elasticity parameters associated with the work equation. Mandatory names are \code{alpha} and \code{beta}.
#'                        \item \strong{\code{times_elasticities}} List. Optional. Value of the elasticity parameters associated with freely allocated times. Must be in the same order as \code{free_times}.
#'                        \item \strong{\code{goods_elasticities}} List. Optional. Value of the elasticity parameters associated with freely allocated goods. Must be in the same order as \code{free_goods}.
#'                        \item \strong{\code{sig}} List. Optional. Value of the sigma parameters (standard deviations) of the error covariance matrix. One element per equation. If not specified, inferred from the observed errors.
#'                        \item \strong{\code{rho}} List. Optional. Value of the correlation parameters of the upper diagonal of the error covariance matrix. Must be ordered horizontally. Only used if \code{sig} is specified. If \code{rho} is not specified, correlations are assumed to be zero.
#'                        \item \strong{\code{cholesky}} List. Optional. Value of the lower Cholesky decomposition triangular matrix elements (L11, L21, L22, L31, L32, L33, ...). If specified, overrides \code{sig} and \code{rho}.
#'                        \item \strong{\code{Tc}} Numeric vector. Committed time for each observation.
#'                        \item \strong{\code{Ec}} Numeric vector. Committed expenses for each observation.
#'                        \item \strong{\code{w}} Numeric vector. Wage rate for each observation.
#'                        \item \strong{\code{tau}} Numeric vector or scalar. Time budget (typically 168 hours per week).
#'                        \item \strong{\code{componentName}} Character. Name given to model component. If not provided by the user, Apollo will set the name automatically.
#'                       }
#'
#' @param functionality Character. Setting instructing Apollo what processing to apply to the likelihood function. This is in general controlled by the functions that call \code{apollo_probabilities}, though the user can also call \code{apollo_probabilities} manually with a given functionality for testing/debugging. Possible values are:
#'                      \itemize{
#'                        \item \strong{\code{"estimate"}}: For model estimation, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"conditionals"}}: For conditionals, produces likelihood of the full model, at the level of individual inter-individual draws.
#'                        \item \strong{\code{"output"}}: Prepares output for post-estimation reporting.
#'                        \item \strong{\code{"prediction"}}: For model prediction, produces probabilities for individual alternatives and individual model components (if multiple components are present) at the level of an observation, after averaging across draws.
#'                        \item \strong{\code{"preprocess"}}: Prepares likelihood functions for use in estimation.
#'                        \item \strong{\code{"raw"}}: For debugging, produces probabilities of all alternatives and individual model components at the level of an observation, at the level of individual draws.
#'                        \item \strong{\code{"report"}}: Prepares output summarising model and choiceset structure.
#'                        \item \strong{\code{"validate"}}: Validates model specification, produces likelihood of the full model, at the level of individual decision-makers, after averaging across draws.
#'                        \item \strong{\code{"get_covar"}}: Returns the covariance matrix, correlation matrix, and standard deviations computed from the residuals.
#'                      }
#' @return For estimation functionalities, an array of likelihood for each individual. For \code{"get_covar"}, a list with components \code{covar}, \code{corr}, and \code{sigma}. For \code{"prediction"}, a matrix with predicted values and values of time.
#' @examples
#' \dontrun{
#' # Single equation model (work only) using the alpha-beta parameterization
#' library(apollo)
#' library(mtuem)
#'
#' apollo_initialise()
#' apollo_control <- list(modelName = "maed-1eq-ab", indivID = "PeID")
#'
#' database <- arrow::read_parquet("data/maed.parquet.gzip")
#' database <- database[database$EcI > 0, ]
#'
#' apollo_beta <- c(alpha = 0.3, beta = 0.3)
#' apollo_fixed <- c()
#'
#' apollo_inputs <- apollo_validateInputs()
#'
#' apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
#'   apollo_attach(apollo_beta, apollo_inputs)
#'   on.exit(apollo_detach(apollo_beta, apollo_inputs))
#'
#'   P <- list()
#'   work_elasticities <- list(alpha = alpha, beta = beta)
#'
#'   mtuem_settings <- list(
#'     work_times = c("Tw"),
#'     work_elasticities = work_elasticities,
#'     Tc = Tc, Ec = EcI, w = w, tau = 168
#'   )
#'
#'   P[["model"]] <- mtuem_ab_likelihood(mtuem_settings, functionality)
#'   P <- apollo_prepareProb(P, apollo_inputs, functionality)
#'   return(P)
#' }
#'
#' model <- apollo_estimate(apollo_beta, apollo_fixed,
#'                           apollo_probabilities, apollo_inputs)
#' }
#' @export
#'
mtuem_ab_likelihood <- function(mtuem_settings, functionality="estimate"){
  # Rename input if necessary
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE))

  default = list(
    free_times = list(),
    free_goods = list(),
    goods_cost = 1,
    times_elasticities = list(),
    goods_elasticities = list(),
    sig = list(),
    rho = list(),
    cholesky = list()
  )

  tmp <- names(default)[!(names(default) %in% names(mtuem_settings))]
  for (i in tmp) mtuem_settings[[i]] <- default[[i]]
  rm(tmp)

  # Copy variables from list to environment
  for(i in 1:length(mtuem_settings)) assign(names(mtuem_settings)[i], mtuem_settings[[i]])
  N <- nrow(apollo_inputs$database)

  ##### PARAMETROS COMUNES

  free_times <- unlist(free_times)
  free_goods <- unlist(free_goods)
  sig <- unlist(sig)
  rho <- unlist(rho)
  cholesky <- unlist(cholesky)

  flag_times <- length(times_elasticities)>0
  flag_goods <- length(goods_elasticities)>0
  estimate_sig <- length(sig) > 0
  estimate_rho   <- length(rho) > 0
  estimate_cholesky <- length(cholesky) > 0

  if (estimate_cholesky) {
    M <- 1 + length(free_times) + length(free_goods)
    if (length(cholesky) != M * (M + 1) / 2) {
      stop(paste0("The length of cholesky (", length(cholesky), ") must match M * (M + 1) / 2 where M is the number of equations (", M, ")."))
    }
    L_mat <- matrix(0, nrow = M, ncol = M)
    idx <- 1
    for (r in 1:M) {
      for (c in 1:r) {
        L_mat[r, c] <- cholesky[idx]
        idx <- idx + 1
      }
    }
    Sigma <- L_mat %*% t(L_mat)
    sig <- sqrt(diag(Sigma))
    D_inv_sqrt <- diag(1 / sig, nrow = M)
    rho <- D_inv_sqrt %*% Sigma %*% D_inv_sqrt
    diag(rho) <- 1
    estimate_sig <- TRUE
    estimate_rho <- TRUE
  }

  if (estimate_rho && !estimate_cholesky) {
    tmp <- rho
    rho <- diag(rep(1, length(sig)))
    #rho[upper.tri(rho, diag = FALSE)] <- rho[lower.tri(rho, diag = FALSE)]  <- tmp
    rho[upper.tri(rho, diag = FALSE)] <- tmp
    rho[lower.tri(rho, diag = FALSE)] <- t(rho)[lower.tri(rho, diag = FALSE)]
    rm(tmp)
  }

  if(functionality=="preprocess"){
    preproc_settings <- list(componentName = "..", gradient = FALSE)
    if(!is.null(mtuem_settings$componentName)){
      preproc_settings$componentName <- mtuem_settings$componentName
    } else if(!is.null(mtuem_settings$componentName2)){
      preproc_settings$componentName <- mtuem_settings$componentName2
    }
    return(preproc_settings)
  }

  # -------------- #
  #### VALIDATE ####
  # -------------- #
  if(functionality %in% c("validate")){
    return(invisible( rep(1, N) ))
  }

  # ------------- #
  #### ZERO LL ####
  # ------------- #
  if(functionality=="zero_LL"){
    ans <- rep(NA, N)
    return(ans)
  }

  # ------------------------------------ #
  #### ESTIMATE, CONDITIONALS AND RAW ####
  # ------------------------------------ #
  if(functionality %in% c("estimate", "conditionals", "raw", "get_covar")){
    tw_opt <- get_tw_albe(work_elasticities, tau, Tc, Ec, w)
    ti_opt <- get_ti_albe(times_elasticities, work_elasticities$beta, tw_opt, tau, Tc)
    xj_opt <- get_xi_albe(goods_elasticities, goods_cost, work_elasticities$alpha, tw_opt, Ec, w)
    opt <- cbind(tw_opt, unlist(ti_opt), unlist(xj_opt))
    colnames(opt) <- c(work_times, free_times, free_goods)
    obs <- as.matrix(apollo_inputs$database[, colnames(opt)] )
    err <- obs - opt

    if (functionality == "get_covar") {
      return(list(
        covar = stats::cov(err, use = "complete.obs"),
        corr =  stats::cor(err, use = "complete.obs"),
        sigma = sqrt(diag(stats::cov(err, use = "complete.obs")))))
    }

    err_ll <- err
    #err_ll[is.na(err_ll)] <- 0

    if (!estimate_sig) {
      sig <- stats::cov(err, use = "complete.obs")
      sig <- sqrt(diag(sig))
    }

    if (!(flag_times | flag_goods)) {
      if (any(is.na(sig)) || any(sig <= 1e-6) || any(is.infinite(sig))) {
        return(rep(0, N))
      }
      mu = err_ll/sig
      ll = -0.5*mu^2 -log(sig) -0.5*log(2*base::pi)
    } else {

      if (!estimate_rho) {
        if (estimate_sig) {
          rho <- diag(rep(1, length(sig)))
        } else {
          rho <- stats::cor(err, use = "complete.obs")
          rho_pd <- Matrix::nearPD(rho, corr = TRUE, keepDiag = TRUE)
          rho    <- as.matrix(rho_pd$mat)
        }
      }

      if (any(is.na(sig)) || any(sig <= 1e-6) || any(is.infinite(sig)) || 
          any(is.na(rho)) || any(is.infinite(rho))) {
        return(rep(0, N))
      }

      eigs <- eigen(rho, symmetric = TRUE, only.values = TRUE)$values
      if (any(is.na(eigs)) || any(eigs <= 1e-6)) {
        return(rep(0, N))
      }

      mu = sweep(err_ll, MARGIN = 2, sig, "/")

      cond <- get_cond_err(mu, rho)
      cond_mu  <- cond$cond_mu
      cond_sd  <- cond$cond_sd

      term = -0.5*(cond_mu^2)
      ll = sweep(term ,  MARGIN = 2, log(sig * sqrt(cond_sd)), "-") - 0.5*log(2*base::pi)
    }

    if(is.matrix(ll)) ll <- rowSums(ll)
    L <- exp(ll)
    return(L)
  }

  # ------------ #
  #### OUTPUT ####
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- mtuem_ab_likelihood(mtuem_settings, functionality="estimate")
    if (functionality == "report") {
      tw_opt <- get_tw_albe(work_elasticities, tau, Tc, Ec, w)
      vot <- get_values_of_time_albe(tw_opt, work_elasticities, Tc, Ec, w)
      colnames(vot) <- c("VoL", "VTAW")
      vot <- colMeans(vot)
      cat("VoL:", vot[1], "\n")
      cat("VTAW:", vot[2], "\n")
    }

    return(ans)
  }

  # ---------------- #
  #### PREDICTION ####
  # ---------------- #
  if(functionality=="prediction"){
    tw_opt <- get_tw_albe(work_elasticities, tau, Tc, Ec, w)
    colnames(tw_opt) <- work_times
    opt <- tw_opt

    if (flag_times) {
      ti_opt <- get_ti_albe(times_elasticities, work_elasticities$beta, tw_opt, tau, Tc)
      colnames(ti_opt) <- free_times
      ti_other = matrix(tau - Tc - rowSums(ti_opt) - tw_opt, ncol=1)
      colnames(ti_other) <- c("Tfi")
      opt <- cbind(opt, ti_opt, ti_other)
    } else {
      ti = matrix(tau - Tc - tw_opt, ncol = 1)
      colnames(ti) <- c("T")
      opt <- cbind(opt, ti)
    }

    if (flag_goods) {
      xj_opt <- get_xi_albe(goods_elasticities, goods_cost, work_elasticities$alpha, tw_opt, Ec, w)
      colnames(xj_opt) <- free_goods
      xj_other = matrix(w*tw_opt - Ec - rowSums(xj_opt), ncol=1)
      colnames(xj_other) <- c("Xfj")
      opt <- cbind(opt, xj_opt, xj_other)
    } else {
      xj = matrix(w*tw_opt - Ec, ncol=1)
      colnames(xj) <- c("X")
      opt <- cbind(opt, xj)
    }

    vot <- get_values_of_time_albe(tw_opt, work_elasticities, Tc, Ec, w)
    colnames(vot) <- c("VoL", "VTAW")
    opt <- cbind(opt, vot)

    return(opt)
  }

  if(functionality=="gradient"){
    return(list(like = NA, grad = NA))
  }

  # End of function
  stop("Invalid value of argument 'functionality'")
}
