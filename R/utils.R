#' Function for getting optimal time allocated to work, with Theta-Phi formulation
#' @keywords internal
#' @export
get_tw_thph <- function(work_elasticities, tau, Tc, Ec, w) {
  ec   <- Ec / (w * (tau-Tc))
  base <- work_elasticities$Theta + work_elasticities$Phi + work_elasticities$thw
  thetaphiec <- work_elasticities$Phi + work_elasticities$thw + (work_elasticities$Theta + work_elasticities$thw)* ec
  sqrt_term <-  sqrt(thetaphiec^2 - 4*work_elasticities$thw*ec*base)
  tw_opt <- 0.5 * (tau - Tc) * (thetaphiec + sqrt_term) / base
  tw_opt <- matrix(tw_opt, ncol = 1)
  return(tw_opt)
}

#' Function for getting optimal time allocated to activities, with Theta-Phi formulation
#' @keywords internal
#' @export
get_ti_thph <- function(times_elasticities, Theta, Tw, tau, Tc) {
  ti_opt <- sapply(unlist(times_elasticities), function(x) (x / Theta) * (tau - Tw - Tc))
  return(ti_opt)
}

#' Function for getting optimal allocation of expenses to goods, with Theta-Phi formulation
#' @keywords internal
#' @export
get_xi_thph <- function(goods_elasticities, goods_cost, Phi, Tw, Ec, w) {
  xj_opt <- sapply(unlist(goods_elasticities) / unlist(goods_cost), function(x) (x / Phi) * (w*Tw - Ec))
  return(xj_opt)
}

#' Function for getting values of time, with Theta-Phi formulation
#' @keywords internal
#' @export
get_values_of_time_thph <- function(tw, work_elasticities, Tc, Ec, w) {
    cteVoLi  <- (w*tw - Ec) / (168-tw-Tc)
    coefVoL <-  work_elasticities$Theta / work_elasticities$Phi

    #cteVTAWi <- (w*tw - Ec) / (tw)
    #coefVTAW <- work_elasticities$thw / work_elasticities$Phi

    VoL  <- cteVoLi * coefVoL
    VTAW <- VoL - w

    return(cbind(VoL, VTAW))
  }

#' Function for getting optimal time allocated to work, with alpha-beta formulation
#' @keywords internal
#' @export
get_tw_albe <- function(work_elasticities, tau, Tc, Ec, w) {
  ec   <- Ec / (w * (tau-Tc))
  betaalphaec <- work_elasticities$beta + work_elasticities$alpha*ec
  sqrt_term <-  sqrt(betaalphaec^2 - ec*(2*work_elasticities$alpha+2*work_elasticities$beta-1))
  tw_opt <- (tau - Tc) * (betaalphaec + sqrt_term)
  tw_opt <- matrix(tw_opt, ncol = 1)
  return(tw_opt)
}

#' Function for getting optimal time allocated to activities, with alpha-beta formulation
#' @keywords internal
#' @export
get_ti_albe <- function(times_elasticities, beta, Tw, tau, Tc) {
  ti_opt <- sapply(unlist(times_elasticities), function(x) (x / (1-2*beta)) * (tau - Tw - Tc))
  return(ti_opt)
}

#' Function for getting optimal allocation of expenses, with alpha-beta formulation
#' @keywords internal
#' @export
get_xi_albe <- function(goods_elasticities, goods_cost, alpha, Tw, Ec, w) {
  xj_opt <- sapply(unlist(goods_elasticities) / unlist(goods_cost), function(x) (x / (1-alpha)) * (w*Tw - Ec))
  return(xj_opt)
}

#' Function for getting values of time, with alpha-beta formulation
#' @keywords internal
#' @export
get_values_of_time_albe <- function(tw, work_elasticities, Tc, Ec, w) {
    cteVoLi  <- (w*tw - Ec) / (168-tw-Tc)
    coefVoL <-  (1-2*work_elasticities$beta) / (1-2*work_elasticities$alpha)

    #cteVTAWi <- (w*tw - Ec) / (tw)
    #coefVTAW <- (2*work_elasticities$beta + 2*work_elasticities$alpha - 1)/ (1-2*work_elasticities$alpha)

    VoL  <- cteVoLi * coefVoL
    VTAW <- VoL - w

    return(cbind(VoL, VTAW))
  }

#' Function for getting conditional normal errors
#' @keywords internal
#' @export
get_cond_err <- function(mu, rho) {
  neq <- ncol(mu)
  conditional_mu        <-  matrix(0, nrow = nrow(mu), ncol = neq)
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

