#' Procedure to estimate from multiple starting points.
#'
#' The implementation is a wrapper of apollo.
#'
#' @param apollo_beta Named numeric vector. Names and values for parameters. Used by \code{apollo_estimate}
#' @param apollo_fixed Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation. Used by \code{apollo_estimate}
#' @param apollo_probabilities Function. Returns probabilities of the model to be estimated. Used by \code{apollo_estimate}. Must receive three arguments:
#'                          \itemize{
#'                            \item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
#'                            \item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
#'                            \item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
#'                          }
#' @param apollo_inputs List grouping most common inputs. Created by function \link{apollo_validateInputs}.
#' @param customMultistart_settings List of settings for multi start. Mimics \code{apollo_multiStart}. The following settuings may be set:
#'                                \itemize{
#'                                  \item \strong{\code{apolloBetaMax}}: Named numeric vector with the maximum value to be sampled for each apollo_beta. Default: apollo_beta + 1
#'                                  \item \strong{\code{apolloBetaMin}}: Named numeric vector with the minimum value to be sampled for each apollo_beta. Default: apollo_beta - 1
#'                                  \item \strong{\code{nCandidates}}: Numneric Scalar. Number of candidates to be sampled. Default: 100
#'                          }
#' @param estimation_settings List. Contains settings for \code{apollo_estimate}. If not specified, default will be used.
#' @param first_em Boolean. If TRUE a few iterations of Expectation-Maximization algorithm will be run. Only use for latent class models. Default: FALSE
#' @param em_iter_max Numeric vector. Number of iterations of EM algorithm to be run. Only valid if \code{first_em} is set to TRUE. Default: 5.
#' @param non_saddle Boolean. If TRUE the procedure will inspect if the estimation result is at least a local optimum to be accepted as a best solution. If FALSE saddle points will also be accepted. Default: TRUE
#' @param verbose Boolean. If TRUE the estimation procedure will be printed in the console. Default FALSE.
#' @param normalization Named list. If not NULL (default) normalizes the parameters. Each list must contains the name of the respective parameters inside \code{work_elasticities}, \code{times_elasticities}, \code{good_elasticities}, \code{sig} and \code{rho} as well as "normalization" = "theta_phi" or "alpha_beta".
#' @param nClass Named list. Default 1. Number of classes.
#'
#' @return a named list containing the best model, the best log-likelihood and all estimated models.
#'        \itemize{
#'          \item \strong{\code{"best_model"}}: The apollo object of the best candidate model.
#'          \item \strong{\code{"best_ll"}}: The value of the log-likelihood of the best candidate model.
#'          \item \strong{\code{"models"}}: A named list with all candidate models
#'        }
#'
#' @export
customMultiStart <- function(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                              customMultistart_settings = NA,
                              estimation_settings = NA,
                              first_em = F,
                              em_iter_max = 5,
                              non_saddle = T,
                              verbose = F,
                              normalization = NA,
                             nClass = 1) {

  default_customMultistart_settings  = list(
      apolloBetaMax = apollo_beta + 1,
      apolloBetaMin = apollo_beta - -1,
      nCandidates = 100
  )

  default_normalization = list(
    normalization = "none",
    work_elasticities = c(),
    times_elasticities = c(),
    goods_elasticities = c(),
    sig = c(),
    rho = c()
  )

  if (length(customMultistart_settings) == 1 && is.na(customMultistart_settings)) customMultistart_settings <- default_customMultistart_settings
  tmp <- names(default_customMultistart_settings)[!(names(default_customMultistart_settings) %in% names(customMultistart_settings))]
  for (i in tmp) customMultistart_settings[[i]] <- default_customMultistart_settings[[i]]
  rm(tmp)

  if (any(!is.na(normalization))) {
    tmp <- names(default_normalization)[!(names(default_normalization) %in% names(normalization))]
    for (i in tmp) normalization[[i]] <- default_normalization[[i]]
    rm(tmp)
  }

  default_estimation_settings = list(
    writeIter = FALSE,
    maxIterations=500,
    estimationRoutine = "bgw",
    hessianRoutine = "maxLik",
    bgw_settings = list(maxFunctionEvals = 1000)
  )

  if (length(estimation_settings) == 1 && is.na(estimation_settings)) estimation_settings <- default_estimation_settings
  tmp <- names(default_estimation_settings)[!(names(default_estimation_settings) %in% names(estimation_settings))]
  for (i in tmp) estimation_settings[[i]] <- default_estimation_settings[[i]]
  rm(tmp)
  
  estimation_settings[["scaleAfterConvergence"]] <- FALSE

  apolloBetaMax = customMultistart_settings[["apolloBetaMax"]]
  apolloBetaMin = customMultistart_settings[["apolloBetaMin"]]
  nCandidates   = customMultistart_settings[["nCandidates"]]
  estimationRoutine = estimation_settings[["estimationRoutine"]]

  apollo_beta_original <- apollo_beta
  beta_matrix <- t(matrix(rep(apollo_beta_original, nCandidates), ncol = nCandidates))
  colnames(beta_matrix) <- names(apollo_beta_original)

  apollo_variables <- names(apollo_beta_original)[!names(apollo_beta_original) %in% apollo_fixed]
  rand <- stats::runif(nCandidates * length(apollo_variables), apolloBetaMin, apolloBetaMax)
  rand <- matrix(rand, nrow = nCandidates, ncol = length(apollo_variables))
  beta_matrix[,apollo_variables] <- beta_matrix[,apollo_variables] + rand

  normalize <- function(beta_matrix, normalization, work_elasticities, times_elasticities, goods_elasticities, sig, rho) {
    if (normalization== "alpha_beta") {
      tmp <- apollo_variables[(apollo_variables %in% work_elasticities)]
      mask <- beta_matrix[, tmp] >= 0.5
      beta_matrix[, tmp][mask] <- 0.49
      rm(tmp)

      tmp <- apollo_variables[apollo_variables %in% c(times_elasticities, goods_elasticities, sig)]
      mask <- beta_matrix[, tmp] <= 0
      beta_matrix[, tmp][mask] <- 0.001
      rm(tmp)

      # TODO esto no funciona con clases aún

      al <- work_elasticities[startsWith("alpha", work_elasticities)]
      be <- work_elasticities[startsWith("beta", work_elasticities)]


      cap_times_elasticities <- 1-2*beta_matrix[,be]
      cap_goods_elasticities <- 1-2*beta_matrix[,al]

     } else {
      tmp <- apollo_variables[
        (apollo_variables %in% c(work_elasticities, times_elasticities, goods_elasticities, sig)) & apollo_variables != "thw"
      ]
      mask <- beta_matrix[, tmp] <= 0
      beta_matrix[, tmp][mask] <- 0.001
      rm(tmp)

      th <- work_elasticities[startsWith(work_elasticities, "Theta" )]
      ph <- work_elasticities[startsWith(work_elasticities, "Phi" )]
      if (th %in% apollo_variables) {
        cap_times_elasticities <- beta_matrix[,th]
      } else {
        cap_times_elasticities <- 1
      }

      if (ph %in% apollo_variables) {
        cap_goods_elasticities <- beta_matrix[,ph]
      } else {
        cap_goods_elasticities <- 1
      }
    }

    if (length(times_elasticities)>0) {
        rsum <- rowSums(beta_matrix[,times_elasticities])
        mask <- rsum > cap_times_elasticities
        rand <- stats::runif(sum(mask),0,1)
        beta_matrix[mask, times_elasticities] <- sweep(
          beta_matrix[mask, times_elasticities],
          1,
          cap_times_elasticities / (rsum[mask] + rand),
          "*" )
    }

    if (length(goods_elasticities)>0) {
      rsum <- rowSums(beta_matrix[,goods_elasticities])
      mask <- rsum > cap_goods_elasticities
      rand <- runif(sum(mask),0,1)
      beta_matrix[mask, goods_elasticities] <- sweep(
        beta_matrix[mask, goods_elasticities],
        1,
        cap_goods_elasticities[mask] / (rsum[mask] + rand),
        "*" )
    }

    tmp <- apollo_variables[startsWith(apollo_variables,"ph" )]
    beta_matrix[, tmp] <- 0
    rm(tmp)

    return(beta_matrix)
  }

  if (normalization$normalization %in% c("alpha_beta", "theta_phi")) {
    norm = normalization$normalization
    work_elasticities = normalization$work_elasticities
    times_elasticities = normalization$times_elasticities
    goods_elasticities = normalization$goods_elasticities
    sig = normalization$sig
    rho = normalization$rho

    if (nClass <= 1) {
      # TODO seguir desde aqui
      beta_matrix <- normalize(
        beta_matrix,
        norm,
        work_elasticities,
        times_elasticities,
        goods_elasticities,
        sig,
        rho
      )
    } else {
      for (s in 1:nClass) {
        if (!is.null(work_elasticities)) {work_elasticities <- paste0(work_elasticities,"_",s)}
        if (!is.null(times_elasticities)) {times_elasticities <- paste0(times_elasticities,"_",s)}
        if (!is.null(goods_elasticities)) {goods_elasticities <- paste0(goods_elasticities,"_",s)}
        if (!is.null(sig)) {sig <- paste0(sig,"_",s)}
        if (!is.null(rho)) {rho <- paste0(rho,"_",s)}
        beta_matrix <- normalize(
          beta_matrix,
          norm,
          work_elasticities,
          times_elasticities,
          goods_elasticities,
          sig,
          rho
        )
      }
    }
  }

  # initial testing
  works <- rep(TRUE, nCandidates)
  indexes <- 1:nCandidates

  cat('Testing candidates...')
  for (i in indexes) {
    suppressWarnings({ P <- apollo_probabilities(beta_matrix[i, ], apollo_inputs, functionality="raw") })
    ll <-sum(log(unlist(P)))
    if(is.na(ll)) {works[i] <- F}
  }
  cat('Done...', sum(works), 'candidates available\n')

  indexes <- indexes[works]
  models <- list()
  best_model <- NA
  best_ll <- -Inf

  
  for (i in indexes) {
    if (verbose) {
      ##### VERBOSE START ####
      cat("Candidate",i,"...\n")
      apollo_beta <- beta_matrix[i,]
      flag <- F
      suppressWarnings({ try({
        if (first_em) {
          cat("Starting EM...\n")

          model <- apollo::apollo_lcEM(
            apollo_beta, 
            apollo_fixed, 
            apollo_probabilities, 
            apollo_inputs,
            lcEM_settings = list(EMmaxIterations = em_iter_max, postEM = 0),
            estimate_settings = list(scaleAfterConvergence = F))

          apollo_beta <- model$estimate
        }
        cat("Starting solver...\n")
        models[[i]] <- apollo::apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                                               estimate_settings= estimation_settings)
        cat("LL:", -round(models[[i]]$maximum, 2),"...\n")
        flag <- T
        cat("Done...\n")
      }, silent = TRUE)})
      ##### VERBOSE END ####
    } else {
      ### NOT VERBOSE START###
      cat("Candidate",i,"...")
      apollo_beta <- beta_matrix[i,]
      flag <- F
      suppressWarnings({ try({
        if (first_em) {
          cat("Starting EM...")
          invisible(utils::capture.output({
            model <- apollo::apollo_lcEM(
              apollo_beta, 
              apollo_fixed, 
              apollo_probabilities, 
              apollo_inputs,
              lcEM_settings = list(EMmaxIterations = em_iter_max, postEM = 0),
              estimate_settings = list(scaling = F))
          }))
          apollo_beta <- model$estimate
          cat("LL:", -round(model$maximum, 2),"...")
        }
        cat("Starting solver...")
        invisible(utils::capture.output({
          models[[i]] <- apollo::apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                                                 estimate_settings= estimation_settings)
        }))
        
        cat("LL:", -round(models[[i]]$maximum, 2),"...")
        flag <- T
        cat("Done...")
      }, silent = TRUE)})
      ### NOT VERBOSE END###
    }
  
    if (flag) {
      code = models[[i]]$code
      is_bgw = estimationRoutine=="bgw"
      is_bfgs = estimationRoutine=="bfgs"
      eigValue = models[[i]]$eigValue
      if (length(eigValue) == 0) { eigValue <- Inf }
      ll = models[[i]] $maximum

      if (non_saddle) {
        flag_expression <- (code==4&is_bgw|code==0&is_bfgs) & (eigValue<0) & (ll>best_ll)
      } else {
        flag_expression <- (code==4&is_bgw|code==0&is_bfgs) & (ll>best_ll)
      }

      if (flag_expression) {
        best_ll <- models[[i]]$maximum
        best_model <- models[[i]]
        cat("Best found")
      }
      cat('\n')
    } else {
      cat("Failed\n")
      models[[i]] <- NA
    }
  }

  out = list(
    best_model = best_model,
    best_ll = best_ll,
    models = models
  )

  return(out)
}



