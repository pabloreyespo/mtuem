#' Likelihood function for the parameter estimation of a Microeconomic Time-Use-Expenditure Model.
#'
#' The implementation is based on Jara-Diaz et.al. (2008) consisting of the maximization of a Cobb-Douglass utility function dependent on Time-use (T) and Goods Consumption (X).
#' The model is subject to time budget, monetary budget and technical time-use-expenditure constraints (committed time and expenditures).
#'
#' @param mtuem_settings List of arguments to the functions. It must contain the following. User input is required for all settings except those with a default or marked as optional.
#'                      \itemize{
#'                        \item \strong{\code{work_times}} Character. Name of the work time column
#'                        \item \strong{\code{free_times}} Character vector. Optional. Names of the freely allocated times. (Remember to left one out for identification)
#'                        \item \strong{\code{free_goods}} Character vector. Optional. Names of the freely consumed goods.(Remember to left one out for identification)
#'                        \item \strong{\code{goods_cost}} Named list. Optional. Cost associated to the consumption of each good. If not especified assumed as 1. Then Expeneses equal consumption of the good.
#'                        \item \strong{\code{work_elasticities}} Named list. Value of the elasticity parameters associated to work equation. Mandatory names are alpha and beta.
#'                        \item \strong{\code{times_elasticities}} List. Optional. Value of the elasticity parameters associated to freely allocated times. Must be in the same order of \code{free_times}.
#'                        \item \strong{\code{goods_elasticities}} List. Optional. Value of the elasticity parameters associated to freely allocated goods Must be in the same order of \code{free_goods}.
#'                        \item \strong{\code{sig}} List. Optional. Value of the sigma parameters (standard deviations) of the error covariance matrix. One element for equation. If not specified inferred by the observed errors.
#'                        \item \strong{\code{rho}} List. Optional. Value of the correlation parameters of the upper diagonal of the error covariance matrix. Must be ordered horizontally. Only used if \code{sigma} is specified. If rho is not correlations are assumed to be zero.
#'                        \item \strong{\code{Tc}} Numeric vector. Budget for each observation.
#'                        \item \strong{\code{Ec}} Numeric vector. Committed expenses for each observation.
#'                        \item \strong{\code{w}} Numeric vector. Wage rate for each observation.
#'                        \item \strong{\code{tau}} Numeric vector or scalar. Time Budget.
#'                        \item \strong{\code{componentName}}: Character. Name given to model component. If not provided by the user, Apollo will set the name automatically according to the element in \code{P} to which the function output is directed.
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
#'                      }
#' @return array of likelihood for each individual
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
    rho = list()
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

  flag_times <- length(times_elasticities)>0
  flag_goods <- length(goods_elasticities)>0
  estimate_sig <- length(sig) > 0
  estimate_rho   <- length(rho) > 0

  if (estimate_rho) {
    tmp <- rho
    rho <- diag(rep(1, length(sig)))
    rho[upper.tri(rho, diag = FALSE)] <- rho[lower.tri(rho, diag = FALSE)]  <- tmp
    rm(tmp)
  }

  get_tw <- function(work_elasticities, tau, Tc, Ec, w) {
    ec   <- Ec / (w * (tau-Tc))
    betaalphaec <- work_elasticities$beta + work_elasticities$alpha*ec
    sqrt_term <-  sqrt(betaalphaec^2 - ec*(2*work_elasticities$alpha+2*work_elasticities$beta-1))
    tw_opt <- (tau - Tc) * (betaalphaec + sqrt_term)
    tw_opt <- matrix(tw_opt, ncol = 1)
    return(tw_opt)
  }

  get_ti <- function(times_elasticities, beta, Tw, tau, Tc) {
    ti_opt <- sapply(unlist(times_elasticities), function(x) (x / (1-2*beta)) * (tau - Tw - Tc))
    return(ti_opt)
  }

  get_xi <- function(goods_elasticities, goods_cost, alpha, Tw, Ec, w) {
    xj_opt <- sapply(unlist(goods_elasticities) / unlist(goods_cost), function(x) (x / (1-alpha)) * (w*Tw - Ec))
    return(xj_opt)
  }

  get_cond_err <- function(mu, rho) {
    neq <- ncol(mu)
    conditional_mu        <-  matrix(0, nrow = N, ncol = neq)
    colnames(conditional_mu) <- colnames(mu)
    conditional_sd        <-  rep(1, neq)

    for (j in 2:neq) {
      if (j == 2) {
        conditional_mu[,2] <- rho[2,1] * mu[,1]
        conditional_sd[2]  <- 1 - rho[2,1]^2
      } else if (j ==3) {
        conditional_mu[,3] <- ((rho[2,3]-rho[1,3]*rho[1,2])* mu[,2] + (rho[1,3]-rho[2,3]*rho[1,2])*mu[,1]) / conditional_sd[2]
        conditional_sd[3]  <- 1- (rho[2,3]^2 -2*rho[2,1]*rho[2,3]*rho[1,3] +rho[1,3]^2) / conditional_sd[2]
      } else {
        i <- j-1
        conditional_mu[,j] <- rho[j,i:1] %*% solve(rho[i:1,i:1]) %*% t(mu[,i:1]) #(3|2,1)
        conditional_sd[j]  <- rho[j,j] - rho[j,i:1] %*% solve(rho[i:1,i:1]) %*% rho[i:1,j]
      }
    }

    conditional_mu = sweep(mu - conditional_mu, MARGIN = 2, sqrt(conditional_sd), "/")
    return(list(cond_mu = conditional_mu, cond_sd = conditional_sd))
  }

  get_values_of_time <- function(tw, work_elasticities, Tc, Ec, w) {
      cteVoLi  <- (w*tw - Ec) / (168-tw-Tc)
      coefVoL <-  (1-2*work_elasticities$beta) / (1-2*work_elasticities$alpha)

      #cteVTAWi <- (w*tw - Ec) / (tw)
      #coefVTAW <- (2*work_elasticities$beta + 2*work_elasticities$alpha - 1)/ (1-2*work_elasticities$alpha)

      VoL  <- cteVoLi * coefVoL
      VTAW <- VoL - w

      return(cbind(VoL, VTAW))
    }

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
    tw_opt <- get_tw(work_elasticities, tau, Tc, Ec, w)
    colnames(tw_opt) <- work_times
    opt <- tw_opt

    if (flag_times) {
      ti_opt <- get_ti(times_elasticities, work_elasticities$beta, tw_opt, tau, Tc)
      colnames(ti_opt) <- free_times
      opt <- cbind(opt, ti_opt)
    }

    if (flag_goods) {
      xj_opt <- get_xi(goods_elasticities, goods_cost, work_elasticities$alpha, tw_opt, Ec, w)
      colnames(xj_opt) <- free_goods
      opt <- cbind(opt, xj_opt)
    }

    obs <- as.matrix(apollo_inputs$database[, colnames(opt)] )
    err <- obs - opt

    if (!estimate_sig) {
      sig <- stats::cov(err, use = "pairwise.complete.obs")
      sig <- sqrt(diag(sig))
    }

    if (!(flag_times | flag_goods)) {
      mu = err/sig
      ll = -0.5*mu^2 -log(sig) -0.5*log(2*base::pi)
    } else {

      if (!estimate_rho) {
        if (estimate_sig) {
          rho <- diag(rep(1, length(sig)))
        } else {
          rho <- stats::cor(err, use = "pairwise.complete.obs")
        }
      }

      mu = sweep(err, MARGIN = 2, sig, "/")
      cond <- get_cond_err(mu, rho)
      cond_mu  <- cond$cond_mu
      cond_sd  <- cond$cond_sd

      term = -0.5*(cond_mu^2)
      ll = sweep(term ,  MARGIN = 2, log(sig * sqrt(cond_sd)), "-") - 0.5*log(2*base::pi)
    }

    if(is.matrix(ll)) ll <- rowSums(ll)
    L <- exp(ll)
    if (any(work_elasticities$Phi<0)) {L <- L*0.0001}
    return(L)
  }

  # ------------ #
  #### OUTPUT #### ---> En base a los datos que entrego me hace un reporte, puedo hacerlo pero aun no es necesario
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- mtuem_ab_likelihood(mtuem_settings, functionality="estimate")
    if (functionality == "output") {
      # Compute values of time and print
      tw_opt <- get_tw(work_elasticities, tau, Tc, Ec, w)
      vot <- get_values_of_time(tw_opt, work_elasticities, Tc, Ec, w)
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
    tw_opt <- get_tw(work_elasticities, tau, Tc, Ec, w)
    colnames(tw_opt) <- work_times
    opt <- tw_opt

    if (flag_times) {
      ti_opt <- get_ti(times_elasticities, work_elasticities$beta, tw_opt, tau, Tc)
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
      # TODO xj asume que precio del bien omitido es = 1
      xj_opt <- get_xi(goods_elasticities, goods_cost, work_elasticities$alpha, tw_opt, Ec, w)
      colnames(xj_opt) <- free_goods
      xj_other = matrix(w*tw_opt - Ec - rowSums(xj_opt), ncol=1)
      colnames(xj_other) <- c("Xfj")
      opt <- cbind(opt, xj_opt, xj_other)
    } else {
      xj = matrix(w*tw_opt - Ec, ncol=1)
      colnames(xj) <- c("X")
      opt <- cbind(opt, xj)
    }

    vot <- get_values_of_time(tw_opt, work_elasticities, Tc, Ec, w)
    colnames(vot) <- c("VoL", "VTAW")
    opt <- cbind(opt, vot)

    return(opt)
  }

  if(functionality=="gradient"){
    ec = Ec / (w * (tau-Tc))
    if (F) {
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
