% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/customMultiStart.R
\name{customMultiStart}
\alias{customMultiStart}
\title{Procedure to estimate from multiple starting points.}
\usage{
customMultiStart(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  customMultistart_settings = NA,
  estimation_settings = NA,
  first_em = F,
  em_iter_max = 5,
  non_saddle = T,
  verbose = F
)
}
\arguments{
\item{apollo_beta}{Named numeric vector. Names and values for parameters. Used by \code{apollo_estimate}}

\item{apollo_fixed}{Character vector. Names (as defined in \code{apollo_beta}) of parameters whose value should not change during estimation. Used by \code{apollo_estimate}}

\item{apollo_probabilities}{Function. Returns probabilities of the model to be estimated. Used by \code{apollo_estimate}. Must receive three arguments:
\itemize{
\item \strong{\code{apollo_beta}}: Named numeric vector. Names and values of model parameters.
\item \strong{\code{apollo_inputs}}: List containing options of the model. See \link{apollo_validateInputs}.
\item \strong{\code{functionality}}: Character. Can be either \strong{\code{"components"}}, \strong{\code{"conditionals"}}, \strong{\code{"estimate"}} (default), \strong{\code{"gradient"}}, \strong{\code{"output"}}, \strong{\code{"prediction"}}, \strong{\code{"preprocess"}}, \strong{\code{"raw"}}, \strong{\code{"report"}}, \strong{\code{"shares_LL"}}, \strong{\code{"validate"}} or \strong{\code{"zero_LL"}}.
}}

\item{apollo_inputs}{List grouping most common inputs. Created by function \link{apollo_validateInputs}.}

\item{customMultistart_settings}{List of settings for multi start. Mimics \code{apollo_multiStart}. The following settuings may be set:
\itemize{
\item \strong{\code{apolloBetaMax}}: Named numeric vector with the maximum value to be sampled for each apollo_beta. Default: apollo_beta + 1
\item \strong{\code{apolloBetaMin}}: Named numeric vector with the minimum value to be sampled for each apollo_beta. Default: apollo_beta - 1
\item \strong{\code{nCandidates}}: Numneric Scalar. Number of candidates to be sampled. Default: 100
}}

\item{estimation_settings}{List. Contains settings for \code{apollo_estimate}. If not specified, default will be used.}

\item{first_em}{Boolean. If TRUE a few iterations of Expectation-Maximization algorithm will be run. Only use for latent class models. Default: FALSE}

\item{em_iter_max}{Numeric vector. Number of iterations of EM algorithm to be run. Only valid if \code{first_em} is set to TRUE. Default: 5.}

\item{non_saddle}{Boolean. If TRUE the procedure will inspect if the estimation result is at least a local optimum to be accepted as a best solution. If FALSE saddle points will also be accepted. Default: TRUE}

\item{verbose}{Boolean. If TRUE the estimation procedure will be printed in the console. Default FALSE.}
}
\value{
a named list containing the best model, the best log-likelihood and all estimated models.
\itemize{
\item \strong{\code{"best_model"}}: The apollo object of the best candidate model.
\item \strong{\code{"best_ll"}}: The value of the log-likelihood of the best candidate model.
\item \strong{\code{"models"}}: A named list with all candidate models
}
}
\description{
The implementation is a wrapper of apollo.
}
