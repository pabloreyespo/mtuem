#' Likelihood function for the parameter estimation of a Microeconomic Time-Use-Expenditure Model (additive quadratic parameterization).
#'
#' The implementation maximizes an additively separable quadratic utility function dependent on
#' Time-use (T) and Goods Consumption (X), subject to time budget, monetary budget and technical
#' time-use-expenditure constraints (committed time and expenditures).
#'
#' @details
#' Utility is \eqn{\sum_k (\theta_k T_k - \tfrac{1}{2}\beta_k T_k^2) + \sum_j (\phi_j X_j - \tfrac{1}{2}\eta_j X_j^2)},
#' so marginal utilities are linear and the first order conditions solve in closed form. Writing
#' the aggregates over the freely allocated activities and goods as
#' \eqn{S_\theta=\sum \theta_i/\beta_i}, \eqn{B=\sum 1/\beta_i},
#' \eqn{S_\phi=\sum P_j\phi_j/\eta_j} and \eqn{H=\sum P_j^2/\eta_j}, the optimal work time is
#'
#' \deqn{T_w^* = \frac{\theta_w + w(S_\phi + E_c)/H - (S_\theta - (\tau - T_c))/B}{\beta_w + w^2/H + 1/B}}
#'
#' and the free allocations follow from the two shadow prices
#' \eqn{\mu = (S_\theta - (\tau - T_w^* - T_c))/B} and \eqn{\lambda = (S_\phi - (wT_w^* - E_c))/H} as
#' \eqn{T_i^* = \theta_i/\beta_i - \mu/\beta_i} and \eqn{X_j^* = \phi_j/\eta_j - P_j\lambda/\eta_j}.
#' The values of time are \eqn{VoL = \mu/\lambda} and \eqn{VTAW = VoL - w}.
#'
#' The goods aggregates are price weighted because the monetary budget adds expenditures
#' \eqn{P_jX_j}, not quantities. With the default \code{goods_cost} of 1 the two coincide.
#'
#' Unlike the Cobb-Douglas parameterizations, the aggregates \eqn{S_\theta, B, S_\phi, H} are built
#' from the parameters supplied in \code{times_elasticities}, \code{times_satiations},
#' \code{goods_elasticities} and \code{goods_satiations}, and these must span all freely allocated
#' activities and goods, including any category left out of the estimated equations. Note also that
#' \eqn{\beta_k > 0} and \eqn{\eta_j > 0} are required for concavity and for the aggregates to be
#' well defined.
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
#' correlation matrix, and standard deviations computed from the residuals.
#'
#' @param mtuem_settings List of arguments to the functions. It must contain the following. User input is required for all settings except those with a default or marked as optional.
#'                      \itemize{
#'                        \item \strong{\code{work_times}} Character. Name of the work time column in the database.
#'                        \item \strong{\code{free_times}} Character vector. Optional. Names of the freely allocated time columns. Leave at least one out for identification.
#'                        \item \strong{\code{free_goods}} Character vector. Optional. Names of the freely consumed goods columns. Leave at least one out for identification.
#'                        \item \strong{\code{goods_cost}} Named list. Optional. Cost associated with the consumption of each good. If not specified, assumed as 1 (expenses equal consumption of the good).
#'                        \item \strong{\code{work_elasticity}} Named list with a single element. Baseline marginal utility of work time (\eqn{\theta_w}).
#'                        \item \strong{\code{work_satiation}} Named list with a single element. Satiation parameter of work time (\eqn{\beta_w}), which must be strictly positive.
#'                        \item \strong{\code{times_elasticities}} List. Baseline marginal utilities (\eqn{\theta_i}) of the freely allocated times. Must be in the same order as \code{free_times}.
#'                        \item \strong{\code{times_satiations}} List. Satiation parameters (\eqn{\beta_i}) of the freely allocated times, strictly positive, in the same order as \code{times_elasticities}.
#'                        \item \strong{\code{goods_elasticities}} List. Baseline marginal utilities (\eqn{\phi_j}) of the freely consumed goods. Must be in the same order as \code{free_goods}.
#'                        \item \strong{\code{goods_satiations}} List. Satiation parameters (\eqn{\eta_j}) of the freely consumed goods, strictly positive, in the same order as \code{goods_elasticities}.
#'                        \item \strong{\code{omitted_times}} Character vector. Optional. Names from \code{free_times} whose equation is dropped from the likelihood. Their parameters stay in \eqn{S_\theta} and \eqn{B}, and they are still predicted. Use this when an equation is a budget identity and therefore collinear with the others.
#'                        \item \strong{\code{omitted_goods}} Character vector. Optional. Same for \code{free_goods}.
#'                        \item \strong{\code{sig}} List. Optional. Value of the sigma parameters (standard deviations) of the error covariance matrix. One element per equation. If not specified, inferred from the observed errors.
#'                        \item \strong{\code{rho}} List. Optional. Value of the correlation parameters of the upper diagonal of the error covariance matrix. Must be ordered horizontally. Only used if \code{sig} is specified. If \code{rho} is not specified, correlations are assumed to be zero.
#'                        \item \strong{\code{cholesky}} List. Optional. Value of the lower Cholesky decomposition triangular matrix elements (L11, L21, L22, L31, L32, L33, ...). If specified, overrides \code{sig} and \code{rho}.
#'                        \item \strong{\code{optimal_tw}} Logical. Optional. If TRUE (default), uses the theoretical optimal work time for computing free time and goods allocations. If FALSE, uses the observed work time from the database.
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
#'
#'                      }
#' @return For estimation functionalities, an array of likelihood for each individual.. For \code{"prediction"}, a matrix with predicted values and values of time.
#' @examples
#' \dontrun{
#' # Three equation model (work, one free time, one free good) using the maed dataset
#' library(apollo)
#' library(mtuem)
#'
#' apollo_initialise()
#' apollo_control <- list(modelName = "maed-3eq-aq", indivID = "PeID")
#'
#' database <- mtuem::maed
#' database <- database[database$EcI > 0, ]
#'
#' apollo_beta <- c(theta_w = 0, beta_w = 1, theta_1 = 1, beta_1 = 1,
#'                  phi_1 = 1, eta_1 = 1, sig_Tw = 10, sig_Tf1 = 10, sig_Ef1 = 10)
#' apollo_fixed <- c("eta_1")
#'
#' apollo_inputs <- apollo_validateInputs()
#'
#' apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
#'   apollo_attach(apollo_beta, apollo_inputs)
#'   on.exit(apollo_detach(apollo_beta, apollo_inputs))
#'
#'   P <- list()
#'
#'   mtuem_settings <- list(
#'     work_times = c("Tw"),
#'     free_times = c("Tf1"),
#'     free_goods = c("Ef1"),
#'     goods_cost = list(Ef1 = 1),
#'     work_elasticity = list(theta_w = theta_w),
#'     work_satiation = list(beta_w = beta_w),
#'     times_elasticities = list(theta_1 = theta_1),
#'     times_satiations = list(beta_1 = beta_1),
#'     goods_elasticities = list(phi_1 = phi_1),
#'     goods_satiations = list(eta_1 = eta_1),
#'     sig = list(sig_Tw = sig_Tw, sig_Tf1 = sig_Tf1, sig_Ef1 = sig_Ef1),
#'     Tc = Tc, Ec = EcI, w = w, tau = 168
#'   )
#'
#'   P[["model"]] <- additive_quadratic_mtuem_likelihood(mtuem_settings, functionality)
#'   P <- apollo_prepareProb(P, apollo_inputs, functionality)
#'   return(P)
#' }
#'
#' model <- apollo_estimate(apollo_beta, apollo_fixed,
#'                           apollo_probabilities, apollo_inputs)
#' }
#' @export
#'
additive_quadratic_mtuem_likelihood <- function(mtuem_settings, functionality="estimate"){
  # Rename input if necessary
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE))

  default = list(
    free_times = list(),
    free_goods = list(),
    goods_cost = 1,
    times_elasticities = list(),
    times_satiations = list(),
    goods_elasticities = list(),
    goods_satiations = list(),
    omitted_times = character(0),
    omitted_goods = character(0),
    sig = list(),
    rho = list(),
    cholesky = list(),
    optimal_tw = T
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
  omitted_times <- unlist(omitted_times)
  omitted_goods <- unlist(omitted_goods)
  sig <- unlist(sig)
  rho <- unlist(rho)
  cholesky <- unlist(cholesky)

  # Equations that enter the likelihood. Categories named in omitted_times /
  # omitted_goods keep their parameters in the aggregates S_theta, B, S_phi, H
  # and are still predicted, but their equation is dropped from the likelihood.
  # This is what allows a category whose equation is a budget identity, and
  # therefore collinear with the others, to be excluded without deleting it
  # from the model.
  bad <- setdiff(c(omitted_times, omitted_goods), c(free_times, free_goods))
  if (length(bad) > 0) {
    stop(paste0("omitted_times / omitted_goods name categories that are not in free_times / free_goods: ",
                paste(bad, collapse = ", "), "."))
  }
  eq_names <- setdiff(c(work_times, free_times, free_goods), c(omitted_times, omitted_goods))
  if (length(eq_names) == 0) stop("Every equation has been omitted; nothing left to estimate.")

  flag_times <- length(times_elasticities)>0
  flag_goods <- length(goods_elasticities)>0
  estimate_sig <- length(sig) > 0
  estimate_rho   <- length(rho) > 0
  estimate_cholesky <- length(cholesky) > 0

  if (estimate_cholesky) {
    M <- length(eq_names)
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

  temp <- get_additive_quadratic_essentials(times_elasticities, times_satiations, goods_elasticities, goods_satiations, goods_cost)
  Sth <- temp$Sth
  Sph <- temp$Sph
  B <- temp$B
  H <- temp$H

  if(functionality=="preprocess"){
    preproc_settings <- list(componentName = "..", gradient = FALSE) # TODO verificar que esto haga sentido
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
  if(functionality %in% c("validate")){ # (Engañito)
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
  #### ESTIMATE, CONDITIONALS AND RAW #### ----> Calculan la verosimilitud (La que tengo en likelihoods)
  # ------------------------------------ #
  if(functionality %in% c("estimate", "conditionals", "raw")){
    tw_opt <- get_tw_additive_quadratic(work_elasticity, work_satiation, Sth, Sph, B, H, tau, Tc, Ec, w)
    if (optimal_tw) {
      ti_opt <- get_ti_additive_quadratic(times_elasticities, times_satiations, Sth, B, tau, tw_opt, Tc)
      xj_opt <- get_xi_additive_quadratic(goods_elasticities, goods_satiations, goods_cost, Sph, H, tw_opt, Ec, w)
    } else {
      ti_opt <- get_ti_additive_quadratic(times_elasticities, times_satiations, Sth, B, tau, apollo_inputs$database[, work_times], Tc)
      xj_opt <- get_xi_additive_quadratic(goods_elasticities, goods_satiations, goods_cost, Sph, H, apollo_inputs$database[, work_times], Ec, w)
    }

    opt <- cbind(tw_opt, ti_opt, xj_opt)
    colnames(opt) <- c(work_times, free_times, free_goods)
    opt <- opt[, eq_names, drop = FALSE]

    if (is.null(apollo_inputs$obs_matrix)) {
      obs <- as.matrix(apollo_inputs$database[, colnames(opt)] )
    } else {
      obs <- apollo_inputs$obs_matrix
    }
    err <- obs - opt

    err_ll <- err
    #err_ll[is.na(err_ll)] <- 0

    if (!estimate_sig) {
      sig <- stats::cov(err, use = "complete.obs")
      sig <- sqrt(diag(sig))
    }

    # Single remaining equation -> univariate normal. Driven by how many
    # equations survive omitted_times / omitted_goods, not by whether free
    # categories were parameterised.
    if (ncol(opt) == 1) {
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

      mu = t(t(err_ll) / sig)

      ######################################################################
      # maybe delete
      #eigs <- eigen(rho, symmetric = TRUE, only.values = TRUE)$values
      #if (any(eigs <= 0)) {
        # Return zero likelihood for this parameter vector so the
        # optimizer steps away from the non-PD region:
      #  return(rep(0, N))
      #}
      ####################################################################
      cond <- get_cond_err(mu, rho)
      cond_mu  <- cond$cond_mu
      cond_sd  <- cond$cond_sd

      #term = -0.5*(cond_mu^2)
      #ll = sweep(term ,  MARGIN = 2, log(sig * sqrt(cond_sd)), "-") - 0.5*log(2*base::pi)
      const <- sum(log(sig)) + 0.5 * sum(log(cond_sd)) + 0.5 * ncol(cond_mu) * log(2 * base::pi)
      ll    <- -0.5 * rowSums(cond_mu^2) - const
    }

    if(is.matrix(ll)) ll <- rowSums(ll)
    L <- exp(ll)
    return(L)
  }

  # ------------ #
  #### OUTPUT #### ---> En base a los datos que entrego me hace un reporte, puedo hacerlo pero aun no es necesario
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- additive_quadratic_mtuem_likelihood(mtuem_settings, functionality="estimate")
    if (functionality == "report") {
      # Compute values of time and print
      tw_opt <- get_tw_additive_quadratic(work_elasticity, work_satiation, Sth, Sph, B, H, tau, Tc, Ec, w)
      vot <- get_values_of_time_additive_quadratic(Sth, Sph, B, H, tw_opt, tau, Tc, Ec, w)
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
    tw_opt <- get_tw_additive_quadratic(work_elasticity, work_satiation, Sth, Sph, B, H, tau, Tc, Ec, w)
    colnames(tw_opt) <- work_times
    opt <- tw_opt

    # No residual time category is reported. Unlike the Cobb-Douglas
    # parameterisations, where the aggregate Theta covers activities left out of
    # free_times, here S_theta and B are built only from the activities passed
    # in, so those activities ARE the whole free set and the leftover
    # tau - Tc - sum(Ti*) - Tw* is structurally zero.
    if (flag_times) {
      ti_opt <- get_ti_additive_quadratic(times_elasticities, times_satiations, Sth, B, tau, tw_opt, Tc)
      colnames(ti_opt) <- free_times
      opt <- cbind(opt, ti_opt)
    } else {
      ti = matrix(tau - Tc - tw_opt, ncol = 1)
      colnames(ti) <- c("T")
      opt <- cbind(opt, ti)
    }

    if (flag_goods) {
      xj_opt <- get_xi_additive_quadratic(goods_elasticities, goods_satiations, goods_cost, Sph, H, tw_opt, Ec, w)
      colnames(xj_opt) <- free_goods
      opt <- cbind(opt, xj_opt)
    } else {
      xj = matrix(w*tw_opt - Ec, ncol=1)
      colnames(xj) <- c("X")
      opt <- cbind(opt, xj)
    }

    vot <- get_values_of_time_additive_quadratic(Sth, Sph, B, H, tw_opt, tau, Tc, Ec, w)
    colnames(vot) <- c("VoL", "VTAW")
    opt <- cbind(opt, vot)

    return(opt)
  }

  if(functionality=="gradient"){
    ec = Ec / (w * (tau-Tc))
    if ("Tw" %in% names(mtuem_settings)) {
        # TODO algún dia arreglar los gradientes, aunque la verdad puede que no sea tan necesario
        x <- (PH + TH + thw)
        thetaphiec = PH + thw + (TH + thw)*ec
        aux_sqrt    = sqrt(thetaphiec^2 - 4*thw*ec*x)
        topt_work = (ta-Tc)*(thetaphiec + aux_sqrt) / (2*x)
        mu = (Tw - topt_work)
        ll = -0.5*(mu/sig)^2 -log(sig) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)

        MuparcPH   = (topt_work*(x-aux_sqrt) - thw*ec*(ta-Tc)) / (x *aux_sqrt)
        MuparcTHW  = (topt_work*(x*(1+ec)-aux_sqrt) - (PH + TH + 2*thw)*ec*(ta-Tc)) / (x *aux_sqrt)
        LLparcPH = (mu / sig^2) * MuparcPH
        LLparcTHW  = (mu / sig^2) * MuparcTHW
        LLparcSigma = (1/sig) * ((mu/sig)^2 - 1)

        G <- list()
        G[["PH"]]      =  LLparcPH * L
        G[["thw"]] =  LLparcTHW  * L
        G[["sig"]]   =  LLparcSigma * L
      # output = list(like = L, grad = G)
    } else {
      output = list(like = NA, grad = NA)
    }
    return ( list(like = NA, grad = NA))
  }

  # End of function
  stop("Invalid value of argument 'functionality'")
}
