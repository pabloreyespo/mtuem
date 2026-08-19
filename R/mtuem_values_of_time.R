#' Values of time for a Theta-Phi MTUEM, with class weights and standard errors
#'
#' Computes the value of leisure (VoL) and the value of assigning time to work (VTAW) from the
#' estimates of a model built with \link{mtuem_likelihood}, for a single-class model or for every
#' class of a latent class model.
#'
#' Both values of time are the product of an individual-specific constant and a ratio of
#' structural parameters,
#'
#' \deqn{VoL_i = \frac{w_i T_{wi} - E_{ci}}{\tau - T_{wi} - T_{ci}} \cdot \frac{\Theta}{\Phi},
#'       \qquad
#'       VTAW_i = \frac{w_i T_{wi} - E_{ci}}{T_{wi}} \cdot \frac{\theta_w}{\Phi}}
#'
#' where \eqn{T_{wi}} is the optimal work time implied by the estimates. The sample mean of each
#' constant is taken first, weighted by class membership, so the delta method only has to carry
#' the parameter ratio. In a latent class model each class therefore gets its own weighted
#' constant, and observations that the model assigns to another class barely count towards it.
#'
#' The weights default to the posterior class membership probabilities: from
#' \link[apollo]{apollo_lcConditionals} for a model estimated through \code{apollo_lc}, and from
#' \link{mtuem_get_profiled_sigma} for one estimated through
#' \link{mtuem_lc_profiled_likelihood}. Supply \code{weights} directly to override that choice,
#' for instance with the prior class shares.
#'
#' @param model Object returned by \link[apollo]{apollo_estimate} (or \code{customMultiStart}'s
#'              \code{best_model}).
#' @param apollo_inputs List grouping most common inputs. Created by
#'                      \link[apollo]{apollo_validateInputs}.
#' @param apollo_probabilities Function. Returns the probabilities of the model. Only needed when
#'                             \code{weights} is left to its default and the posterior class
#'                             probabilities have to be recovered.
#' @param nClass Numeric. Number of latent classes. Default: 1, i.e. no class suffix is appended
#'               to the parameter names.
#' @param weights Optional. Class membership weights. Either an N x nClass matrix (one column per
#'                class), a list of nClass numeric vectors of length N, or a single numeric vector
#'                of length N reused for every class. Defaults to the posterior class
#'                probabilities, and to equal weights for a single-class model.
#' @param tau Numeric. Time budget. Default: 168.
#' @param Tc,Ec,w Character. Column names in \code{apollo_inputs$database} holding committed time,
#'                committed expenditure and the wage rate. Pass the same columns used in the
#'                likelihood, e.g. \code{Ec = "EcI"}. When \code{class_specific_budgets} is TRUE,
#'                \code{Tc} and \code{Ec} are read per class as \code{"Tc_1"}, \code{"Ec_1"}, ...
#' @param class_specific_budgets Logical. Whether the committed time and expenditure differ by
#'                               class. Default: FALSE.
#' @param standard_errors Logical. Whether to apply \link[apollo]{apollo_deltaMethod} to the
#'                        parameter ratios. Default: TRUE.
#' @param conf_level Numeric. Confidence level for the reported interval. Default: 0.95.
#' @param silent Logical. Whether to suppress the printed report of means and confidence
#'               intervals. Default: FALSE.
#' @return Named list:
#'         \itemize{
#'           \item \strong{\code{summary}}: Data frame with one row per class and quantity
#'                 (\code{VoL}, \code{VTAW} and the observed \code{wage}), holding the weighted
#'                 mean, its standard error and the confidence interval.
#'           \item \strong{\code{individual}}: Data frame with the optimal work time, both values
#'                 of time and the weight, per observation and class.
#'           \item \strong{\code{delta}}: List with the raw \link[apollo]{apollo_deltaMethod}
#'                 output per class, or NULL when \code{standard_errors} is FALSE.
#'         }
#' @export
mtuem_values_of_time <- function(model, apollo_inputs, apollo_probabilities = NULL,
                                 nClass = 1, weights = NULL, tau = 168,
                                 Tc = "Tc", Ec = "Ec", w = "w",
                                 class_specific_budgets = FALSE,
                                 standard_errors = TRUE, conf_level = 0.95,
                                 silent = FALSE) {

  database <- apollo_inputs$database
  N <- nrow(database)
  b <- model$estimate

  W <- .mtuem_vot_weights(weights, model, apollo_probabilities, apollo_inputs, nClass, N)

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  wage <- as.numeric(database[[w]])

  summ <- NULL
  indiv <- NULL
  deltas <- vector("list", nClass)

  for (s in 1:nClass) {
    suffix <- if (nClass > 1) paste0("_", s) else ""
    pname <- function(x) paste0(x, suffix)

    work_elasticities <- list(
      Theta = .mtuem_vot_par(b, pname("Theta")),
      Phi   = .mtuem_vot_par(b, pname("Phi")),
      thw   = .mtuem_vot_par(b, pname("thw"))
    )

    Tc_s <- as.numeric(database[[if (class_specific_budgets) paste0(Tc, "_", s) else Tc]])
    Ec_s <- as.numeric(database[[if (class_specific_budgets) paste0(Ec, "_", s) else Ec]])

    Tw <- as.numeric(get_tw_thph(work_elasticities, tau, Tc_s, Ec_s, wage))

    cteVoLi  <- (wage * Tw - Ec_s) / (tau - Tw - Tc_s)
    cteVTAWi <- (wage * Tw - Ec_s) / Tw

    coefVoL  <- work_elasticities$Theta / work_elasticities$Phi
    coefVTAW <- work_elasticities$thw / work_elasticities$Phi

    wt <- W[, s]
    cteVoL  <- stats::weighted.mean(cteVoLi,  wt, na.rm = TRUE)
    cteVTAW <- stats::weighted.mean(cteVTAWi, wt, na.rm = TRUE)

    indiv <- rbind(indiv, data.frame(
      class = s,
      Tw = Tw,
      VoL = coefVoL * cteVoLi,
      VTAW = coefVTAW * cteVTAWi,
      weight = wt
    ))

    wage_mean <- stats::weighted.mean(wage, wt, na.rm = TRUE)
    wage_se <- .mtuem_vot_wtd_se(wage, wt)

    if (standard_errors) {
      delta <- apollo::apollo_deltaMethod(model, list(expression = c(
        VoL  = paste0(cteVoL,  "*", pname("Theta"), "/", pname("Phi")),
        VTAW = paste0(cteVTAW, "*", pname("thw"),   "/", pname("Phi"))
      )))
      deltas[[s]] <- delta
      se <- .mtuem_vot_se_column(delta)
      vol_se <- se[1]
      vtaw_se <- se[2]
    } else {
      vol_se <- NA_real_
      vtaw_se <- NA_real_
    }

    summ <- rbind(summ, data.frame(
      class = s,
      quantity = c("VoL", "VTAW", "wage"),
      mean = c(coefVoL * cteVoL, coefVTAW * cteVTAW, wage_mean),
      s.e. = c(vol_se, vtaw_se, wage_se),
      check.names = FALSE
    ))
  }

  summ$lower <- summ$mean - z * summ[["s.e."]]
  summ$upper <- summ$mean + z * summ[["s.e."]]
  rownames(summ) <- NULL

  if (!silent) .mtuem_vot_report(summ, nClass, conf_level)

  return(list(
    summary = summ,
    individual = indiv,
    delta = if (standard_errors) deltas else NULL
  ))
}

#' Read one structural parameter from a model, by name
#' @keywords internal
.mtuem_vot_par <- function(b, name) {
  if (!(name %in% names(b))) {
    stop("Parameter '", name, "' is not among the model estimates. Check nClass and the ",
         "parameter naming: mtuem_values_of_time expects the Theta-Phi parameterization, with ",
         "a '_s' suffix per class when nClass > 1.")
  }
  as.numeric(b[[name]])
}

#' Weighted standard error of a mean
#' @keywords internal
.mtuem_vot_wtd_se <- function(x, wt) {
  ok <- !is.na(x) & !is.na(wt)
  x <- x[ok]; wt <- wt[ok]
  m <- stats::weighted.mean(x, wt)
  v <- sum(wt * (x - m)^2) / (sum(wt) - 1)
  sqrt(v / sum(wt))
}

#' Locate the standard error column of an apollo_deltaMethod result
#' @keywords internal
.mtuem_vot_se_column <- function(delta) {
  col <- grep("^s\\.e\\.", colnames(delta))
  if (length(col) == 0) col <- 3
  as.numeric(delta[, col[1]])
}

#' Resolve the class membership weights used to average the values of time
#' @keywords internal
.mtuem_vot_weights <- function(weights, model, apollo_probabilities, apollo_inputs, nClass, N) {
  if (is.null(weights)) {
    if (nClass == 1) return(matrix(1, nrow = N, ncol = 1))
    weights <- .mtuem_vot_posterior(model, apollo_probabilities, apollo_inputs, nClass)
  }

  if (is.list(weights) && !is.data.frame(weights)) weights <- do.call(cbind, weights)
  if (is.data.frame(weights)) {
    weights <- weights[, !(colnames(weights) %in% c("ID", "Observation")), drop = FALSE]
    weights <- as.matrix(weights)
  }
  if (is.vector(weights)) weights <- matrix(weights, nrow = length(weights), ncol = nClass)

  if (nrow(weights) != N || ncol(weights) != nClass) {
    stop("weights must have ", N, " rows and ", nClass, " columns, but has ",
         nrow(weights), " x ", ncol(weights), ".")
  }
  weights
}

#' Recover posterior class membership probabilities
#' @keywords internal
.mtuem_vot_posterior <- function(model, apollo_probabilities, apollo_inputs, nClass) {
  post <- NULL
  if (!is.null(apollo_probabilities)) {
    post <- tryCatch(
      apollo::apollo_lcConditionals(model, apollo_probabilities, apollo_inputs),
      error = function(e) NULL)
  }

  if (is.null(post)) {
    post <- tryCatch({
      profiled <- mtuem_get_profiled_sigma()
      if (ncol(profiled$posterior) == nClass) profiled$posterior else NULL
    }, error = function(e) NULL)
  }

  if (is.null(post)) {
    stop("Could not recover the posterior class probabilities. Pass apollo_probabilities so ",
         "they can be computed with apollo_lcConditionals, or supply them through the ",
         "'weights' argument.")
  }
  post
}

#' Print the values of time as one line per quantity, with confidence intervals
#' @keywords internal
.mtuem_vot_report <- function(summ, nClass, conf_level) {
  labels <- c(VoL = "VoL              ", VTAW = "VTAW             ",
              wage = "observed wage    ")
  cat(sprintf("
Values of time, %.0f%% confidence intervals
", 100 * conf_level))
  for (s in 1:nClass) {
    if (nClass > 1) cat(sprintf("
CLASS %d:
", s))
    rows <- summ[summ$class == s, ]
    for (i in seq_len(nrow(rows))) {
      label <- labels[[rows$quantity[i]]]
      if (is.null(label)) label <- rows$quantity[i]
      cat(sprintf("  %s %9.3f  [ %9.3f ; %9.3f ]
",
                  label, rows$mean[i], rows$lower[i], rows$upper[i]))
    }
  }
  cat("
")
}
