# Package-local environment used to stash the covariance matrices, posterior
# class weights and class shares recovered by the last successful call to
# mtuem_lc_profiled_likelihood(). Not exported.
.mtuem_profile_env <- new.env(parent = emptyenv())

#' Zero-mean multivariate normal density via Cholesky decomposition
#' @keywords internal
.mtuem_mvn_density <- function(E, Sigma) {
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  dg <- diag(R)
  if (any(!is.finite(dg)) || any(dg <= 1e-8)) return(NULL)

  M <- ncol(E)
  z <- backsolve(R, t(E), transpose = TRUE)
  q <- colSums(z^2)
  logdet <- 2 * sum(log(dg))
  ll <- -0.5 * (q + logdet + M * log(2 * base::pi))
  return(exp(ll))
}

#' Core fixed-point computation for the profiled latent class likelihood
#'
#' Concentrates the class-specific error covariance out of the likelihood via a
#' fixed point on the class-conditional density weights. Returns \code{NULL}
#' whenever the likelihood should be treated as degenerate for the current
#' parameter vector (non-finite residuals, a non-positive-definite covariance,
#' or class collapse), so the caller can fall back to \code{rep(0, N)}.
#' @keywords internal
.mtuem_profile_core <- function(class_settings, classProb, apollo_inputs, em_settings) {
  S <- length(class_settings)
  N <- nrow(apollo_inputs$database)

  E <- lapply(class_settings, mtuem_residuals, apollo_inputs = apollo_inputs)
  if (any(!is.finite(unlist(E)))) return(NULL)

  M <- ncol(E[[1]])
  minWeight <- em_settings$minWeight
  if (is.null(minWeight)) minWeight <- 2 * M

  Pi <- do.call(cbind, classProb)

  # Deterministic unweighted starting covariance, recomputed fresh on every
  # call. This must never be cached/warm-started from a previous evaluation:
  # doing so would make the objective path-dependent and break Apollo's
  # numerical gradients.
  Sigma <- lapply(E, function(e) stats::cov(e, use = "complete.obs"))

  dens <- matrix(NA_real_, nrow = N, ncol = S)

  for (iter in 1:em_settings$maxIterations) {
    for (s in 1:S) {
      d <- .mtuem_mvn_density(E[[s]], Sigma[[s]])
      if (is.null(d)) return(NULL)
      dens[, s] <- d
    }

    num <- Pi * dens
    den <- rowSums(num)
    if (any(!is.finite(den)) || any(den <= 0)) return(NULL)
    W <- num / den
    if (any(colSums(W) < minWeight)) return(NULL)

    Sigma_new <- vector("list", S)
    for (s in 1:S) {
      ws <- W[, s]
      Sigma_new[[s]] <- crossprod(E[[s]] * sqrt(ws)) / sum(ws)
    }

    delta <- max(abs(unlist(Sigma_new) - unlist(Sigma)))
    Sigma <- Sigma_new
    if (delta < em_settings$tol) break
  }

  # Recompute the class densities once with the converged Sigma
  for (s in 1:S) {
    d <- .mtuem_mvn_density(E[[s]], Sigma[[s]])
    if (is.null(d)) return(NULL)
    dens[, s] <- d
  }
  num <- Pi * dens
  P <- rowSums(num)
  if (any(!is.finite(P)) || any(P <= 0)) return(NULL)
  W <- num / P

  return(list(P = P, Sigma = Sigma, W = W, Pi = Pi))
}

#' Stash the converged covariance, posterior weights and class shares from a successful run
#' @keywords internal
.mtuem_stash_profile <- function(core) {
  assign("Sigma", core$Sigma, envir = .mtuem_profile_env)
  assign("W", core$W, envir = .mtuem_profile_env)
  assign("Pi", core$Pi, envir = .mtuem_profile_env)
}

#' Latent class likelihood with the error covariance concentrated out (Theta-Phi parameterization)
#'
#' Implements a latent class MTUEM likelihood in which the class-specific error covariance
#' matrices are profiled out of the likelihood via an inner fixed point, rather than estimated
#' as free parameters. Apollo therefore never sees a single covariance parameter for this model
#' component. \code{mtuem_lc_profiled_likelihood} already combines the classes into a single
#' mixture likelihood, so its result should be assigned directly to a component of \code{P} in
#' \code{apollo_probabilities} (do not pass it through \link[apollo]{apollo_lc}).
#'
#' @details
#' On every call, the class-specific error covariance matrices are re-estimated from scratch
#' starting from the deterministic unweighted covariance of the current residuals (never a warm
#' start from a previous call), so the likelihood surface stays path-independent for Apollo's
#' numerical gradients. The fixed point alternates between computing posterior class weights
#' from the current covariance matrices and re-estimating the covariance matrices as a
#' weighted second moment of the residuals, until the largest change in any covariance element
#' falls below \code{em_settings$tol} or \code{em_settings$maxIterations} is reached.
#'
#' @param lc_settings List of arguments. Must contain:
#'                      \itemize{
#'                        \item \strong{\code{class_settings}}: List of per-class \code{mtuem_settings} lists (same
#'                          shape as used by \link{mtuem_likelihood}), WITHOUT \code{sig}, \code{rho}, or \code{cholesky}.
#'                        \item \strong{\code{classProb}}: List of class probability vectors, as returned by
#'                          \link[apollo]{apollo_classAlloc}.
#'                        \item \strong{\code{em_settings}}: Optional list with \code{tol} (default \code{1e-10}),
#'                          \code{maxIterations} (default \code{50}), and \code{minWeight} (default \code{2 * M}, where
#'                          \code{M} is the number of equations).
#'                        \item \strong{\code{componentName}}: Optional character. Name given to model component.
#'                      }
#' @param functionality Character. Same meaning as in \link{mtuem_likelihood}.
#' @return For estimation functionalities, an array of likelihood for each individual. For
#'         \code{"prediction"}, a matrix with the per-class predicted allocations and the posterior
#'         class weights.
#' @export
mtuem_lc_profiled_likelihood <- function(lc_settings, functionality = "estimate") {
  apollo_inputs <- tryCatch(get('apollo_inputs', envir = parent.frame(), inherits = FALSE),
                             error = function(e) list(silent = FALSE))

  default <- list(
    em_settings = list(),
    componentName = NULL
  )
  tmp <- names(default)[!(names(default) %in% names(lc_settings))]
  for (i in tmp) lc_settings[[i]] <- default[[i]]
  rm(tmp)

  default_em <- list(tol = 1e-10, maxIterations = 50, minWeight = NULL)
  tmp <- names(default_em)[!(names(default_em) %in% names(lc_settings$em_settings))]
  for (i in tmp) lc_settings$em_settings[[i]] <- default_em[[i]]
  rm(tmp)

  class_settings <- lc_settings$class_settings
  classProb <- lc_settings$classProb
  em_settings <- lc_settings$em_settings
  S <- length(class_settings)
  N <- nrow(apollo_inputs$database)

  if (functionality == "preprocess") {
    preproc_settings <- list(componentName = "..", gradient = FALSE)
    if (!is.null(lc_settings$componentName)) {
      preproc_settings$componentName <- lc_settings$componentName
    }
    return(preproc_settings)
  }

  if (functionality %in% c("validate")) {
    return(invisible(rep(1, N)))
  }

  if (functionality == "zero_LL") {
    return(rep(NA, N))
  }

  if (functionality %in% c("estimate", "conditionals", "raw")) {
    core <- .mtuem_profile_core(class_settings, classProb, apollo_inputs, em_settings)
    if (is.null(core)) return(rep(0, N))
    .mtuem_stash_profile(core)
    return(core$P)
  }

  if (functionality %in% c("output", "report")) {
    ans <- mtuem_lc_profiled_likelihood(lc_settings, functionality = "estimate")
    if (functionality == "report") {
      for (s in 1:S) {
        cs <- class_settings[[s]]
        tw_opt <- get_tw_thph(cs$work_elasticities, cs$tau, cs$Tc, cs$Ec, cs$w)
        vot <- get_values_of_time_thph(tw_opt, cs$work_elasticities, cs$Tc, cs$Ec, cs$w)
        colnames(vot) <- c("VoL", "VTAW")
        vot <- colMeans(vot)
        cat("Class", s, "VoL:", vot[1], "\n")
        cat("Class", s, "VTAW:", vot[2], "\n")
      }
    }
    return(ans)
  }

  if (functionality == "prediction") {
    core <- .mtuem_profile_core(class_settings, classProb, apollo_inputs, em_settings)
    if (is.null(core)) return(rep(0, N))
    .mtuem_stash_profile(core)

    opt_list <- list()
    for (s in 1:S) {
      opt_s <- mtuem_likelihood(class_settings[[s]], functionality = "prediction")
      colnames(opt_s) <- paste0(colnames(opt_s), "_class", s)
      opt_list[[s]] <- opt_s
    }
    post <- core$W
    colnames(post) <- paste0("post_class", 1:S)

    out <- do.call(cbind, c(opt_list, list(post)))
    return(out)
  }

  # End of function
  stop("Invalid value of argument 'functionality'")
}

#' Retrieve the covariance recovered by the last profiled latent class likelihood evaluation
#'
#' Returns the class-specific error covariance matrices, standard deviations, correlation
#' matrices, posterior class weights, and class shares recovered by the inner fixed point of
#' the LAST call to \link{mtuem_lc_profiled_likelihood}. Because the covariance is profiled out
#' rather than estimated as a free parameter, this reflects whatever parameter vector Apollo last
#' evaluated \code{apollo_probabilities} at. To retrieve the covariance at the converged
#' estimates, re-evaluate \code{apollo_probabilities} at \code{model$estimate} before calling
#' this function.
#'
#' @return Named list:
#'         \itemize{
#'           \item \strong{\code{Sigma}}: List of class-specific error covariance matrices.
#'           \item \strong{\code{sig}}: List of class-specific standard deviation vectors, \code{sqrt(diag(Sigma_s))}.
#'           \item \strong{\code{rho}}: List of class-specific correlation matrices implied by \code{Sigma}.
#'           \item \strong{\code{posterior}}: N x S matrix of posterior class membership weights.
#'           \item \strong{\code{shares}}: Numeric vector with the average class allocation probability (prior) per class.
#'         }
#' @export
mtuem_get_profiled_sigma <- function() {
  if (!exists("Sigma", envir = .mtuem_profile_env)) {
    stop("No profiled covariance is available yet. Call apollo_probabilities (which in turn ",
         "calls mtuem_lc_profiled_likelihood) at least once before calling mtuem_get_profiled_sigma().")
  }

  Sigma <- get("Sigma", envir = .mtuem_profile_env)
  W     <- get("W", envir = .mtuem_profile_env)
  Pi    <- get("Pi", envir = .mtuem_profile_env)

  sig <- lapply(Sigma, function(s) sqrt(diag(s)))
  rho <- lapply(seq_along(Sigma), function(s) {
    D_inv_sqrt <- diag(1 / sig[[s]], nrow = length(sig[[s]]))
    r <- D_inv_sqrt %*% Sigma[[s]] %*% D_inv_sqrt
    diag(r) <- 1
    r
  })

  return(list(
    Sigma = Sigma,
    sig = sig,
    rho = rho,
    posterior = W,
    shares = colMeans(Pi)
  ))
}

#' Information criteria for a latent class model with a profiled-out error covariance
#'
#' Computes the effective number of free parameters and the corresponding AIC/BIC for a model
#' estimated with \link{mtuem_lc_profiled_likelihood}, counting the profiled-out covariance
#' parameters as estimated parameters (they are fitted by the inner fixed point at the optimum,
#' even though Apollo never sees them as entries of \code{apollo_beta}).
#'
#' @param model Object returned by \link[apollo]{apollo_estimate} (or \link{customMultiStart}'s
#'              \code{best_model}).
#' @param nClass Numeric. Number of latent classes.
#' @param nEq Numeric. Number of equations per class (the dimension of each class-specific
#'            covariance matrix).
#' @return Named list with \code{k} (effective number of free parameters), \code{AIC}, and \code{BIC}.
#' @export
mtuem_profile_ic <- function(model, nClass, nEq) {
  k <- length(model$estimate) - length(model$apollo_fixed) + nClass * nEq * (nEq + 1) / 2
  ll <- model$maximum
  n  <- model$nObs

  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + k * log(n)

  return(list(k = k, AIC = aic, BIC = bic))
}
