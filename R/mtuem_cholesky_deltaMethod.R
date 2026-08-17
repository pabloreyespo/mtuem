#' Function for applying the delta method to Cholesky parameters to retrieve standard deviations and correlations
#'
#' This function computes standard errors, t-ratios, and p-values for the reconstructed error standard deviations
#' (sigmas) and correlations (rhos) from the estimated Cholesky parameters using the delta method.
#'
#' @param model Object. Model returned by apollo estimation routines.
#' @param cholesky_names Character vector. Names of the estimated Cholesky parameters in the order they are used (L11, L21, L22, L31, L32, L33, ...).
#' @param eq_names Character vector. Optional. Names of the equations to label the output sigmas and rhos. If not provided, numeric indices are used.
#' @param deltaMethod_settings_extra List. Optional. Additional settings to pass to \code{apollo_deltaMethod} (e.g. \code{varcov}).
#' @return Data frame with the delta method results (estimates, standard errors, t-ratios, etc.) returned by \code{apollo_deltaMethod}.
#' @export
mtuem_cholesky_deltaMethod <- function(model, cholesky_names, eq_names = NULL, deltaMethod_settings_extra = NULL) {
  K <- length(cholesky_names)
  M <- round((-1 + sqrt(1 + 8 * K)) / 2)
  if (M * (M + 1) / 2 != K) {
    stop("The number of Cholesky parameter names must be of the form M * (M + 1) / 2.")
  }

  L_names <- matrix("", nrow = M, ncol = M)
  idx <- 1
  for (r in 1:M) {
    for (c in 1:r) {
      L_names[r, c] <- cholesky_names[idx]
      idx <- idx + 1
    }
  }

  expressions <- c()

  # First construct all sig expressions
  sig_exprs <- character(M)
  for (i in 1:M) {
    vars_expr <- paste0("(", paste0(L_names[i, 1:i], "^2", collapse = " + "), ")")
    sig_exprs[i] <- paste0("sqrt(", vars_expr, ")")

    label <- if (!is.null(eq_names)) paste0("sig_", eq_names[i]) else paste0("sig_", i)
    expressions[label] <- sig_exprs[i]
  }

  # Then construct all rho expressions
  for (i in 2:M) {
    for (j in 1:(i-1)) {
      cov_expr <- paste0("(", paste0(L_names[i, 1:j], " * ", L_names[j, 1:j], collapse = " + "), ")")
      rho_expr <- paste0(cov_expr, " / (", sig_exprs[i], " * ", sig_exprs[j], ")")

      label <- if (!is.null(eq_names)) paste0("rho_", eq_names[j], "_", eq_names[i]) else paste0("rho_", j, "_", i)
      expressions[label] <- rho_expr
    }
  }

  deltaMethod_settings <- list(expression = expressions)
  if (!is.null(deltaMethod_settings_extra)) {
    for (name in names(deltaMethod_settings_extra)) {
      deltaMethod_settings[[name]] <- deltaMethod_settings_extra[[name]]
    }
  }

  return(apollo::apollo_deltaMethod(model, deltaMethod_settings))
}
