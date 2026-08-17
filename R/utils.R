#' Function for getting optimal time allocated to work, with Theta-Phi formulation
#' @keywords internal
#' @export
#'
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

#' Function for getting the aggregates of the additive quadratic formulation
#'
#' Computes the four sufficient statistics of the additive quadratic MTUEM over
#' the freely allocated activities and goods:
#' \eqn{S_\theta=\sum \theta_i/\beta_i}, \eqn{B=\sum 1/\beta_i},
#' \eqn{S_\phi=\sum P_j\phi_j/\eta_j} and \eqn{H=\sum P_j^2/\eta_j}.
#' The goods aggregates are price weighted because the monetary budget adds
#' expenditures \eqn{P_jX_j}, not quantities.
#'
#' @keywords internal
#' @export
#'
get_additive_quadratic_essentials <- function(times_elasticities, times_satiations, goods_elasticities, goods_satiations, goods_cost = 1) {
  theta <- unlist(times_elasticities)
  beta  <- unlist(times_satiations)
  phi   <- unlist(goods_elasticities)
  eta   <- unlist(goods_satiations)
  P     <- unlist(goods_cost)
  if (length(P) == 1) P <- rep(P, length(eta))

  values <- list(
    Sth = sum(theta / beta),
    Sph = sum(P * phi / eta),
    B   = sum(1 / beta),
    H   = sum(P^2 / eta)
  )
  return(values)
}

#' Function for getting optimal time allocated to work, with additive quadratic formulation
#' @keywords internal
#' @export
get_tw_additive_quadratic <- function(work_elasticity, work_satiation, Sth, Sph, B, H, tau, Tc, Ec, w) {
  num_1 <- unlist(work_elasticity)
  num_2 <- w*(Sph+Ec) / H
  num_3 <- (Sth-(tau-Tc)) / B
  den <- unlist(work_satiation) + w^2/H + 1/B
  tw_opt <- (num_1 + num_2 - num_3) / den
  tw_opt <- matrix(tw_opt, ncol = 1)
  return(tw_opt)
}

#' Function for getting optimal time allocated to activities, with additive quadratic formulation
#' @keywords internal
#' @export
get_ti_additive_quadratic <- function(times_elasticities, times_satiations, Sth, B, tau, Tw, Tc) {
  theta <- unlist(times_elasticities)
  beta  <- unlist(times_satiations)
  Tw    <- as.numeric(unlist(Tw))
  mu    <- (Sth - (tau - Tw - Tc)) / B  # time shadow price, one value per observation
  ti_opt <- outer(rep(1, length(mu)), theta / beta) - outer(mu, 1 / beta)
  return(ti_opt)
}

#' Function for getting optimal allocation of expenses to goods, with additive quadratic formulation
#' @keywords internal
#' @export
get_xi_additive_quadratic <- function(goods_elasticities, goods_satiations, goods_cost, Sph, H, Tw, Ec, w) {
  phi <- unlist(goods_elasticities)
  eta <- unlist(goods_satiations)
  P   <- unlist(goods_cost)
  if (length(P) == 1) P <- rep(P, length(eta))
  Tw  <- as.numeric(unlist(Tw))
  lambda <- (Sph - (w*Tw - Ec)) / H  # money shadow price, one value per observation
  xj_opt <- outer(rep(1, length(lambda)), phi / eta) - outer(lambda, P / eta)
  return(xj_opt)
}

#' Function for getting values of time, with additive quadratic formulation
#' @keywords internal
#' @export
get_values_of_time_additive_quadratic <- function(Sth, Sph, B, H, Tw, tau, Tc, Ec, w) {
    Tw <- as.numeric(unlist(Tw))
    mu = (Sth - (tau - Tw - Tc)) / B
    lambda = (Sph - (w*Tw - Ec)) / H
    VoL  <- mu / lambda
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
      inv_rho <- MASS::ginv(rho[i:1, i:1])
      conditional_mu[,j] <- rho[j,i:1] %*% inv_rho %*% t(mu[,i:1]) #(3|2,1)
      conditional_sd[j]  <- rho[j,j] - rho[j,i:1] %*% inv_rho %*% rho[i:1,j]
    }
  }

  conditional_mu = sweep(mu - conditional_mu, MARGIN = 2, sqrt(conditional_sd), "/")
  return(list(cond_mu = conditional_mu, cond_sd = conditional_sd))
}

