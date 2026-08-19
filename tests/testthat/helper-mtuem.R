# Shared fixtures for the mtuem test battery.
#
# The likelihood functions look up `apollo_inputs` in their calling frame, so
# every test that calls one must define `apollo_inputs` locally. The builders
# below return it as an element of the fixture; tests assign it explicitly.

# Deterministic Theta-Phi fixture: observed allocations are the optimal
# allocations plus normal noise with a known covariance.
mtuem_fixture_thph <- function(n = 200,
                               work_elasticities = list(Theta = 1, Phi = 1, thw = 0),
                               times_elasticities = list(th1 = 0.3, th2 = 0.2),
                               goods_elasticities = list(),
                               tau = 168,
                               seed = 1,
                               sd = 2) {
  set.seed(seed)
  Tc <- stats::runif(n, 40, 80)
  Ec <- stats::runif(n, 50, 150)
  w  <- stats::runif(n, 10, 25)

  tw_opt <- get_tw_thph(work_elasticities, tau, Tc, Ec, w)
  opt <- as.matrix(tw_opt)
  nms <- "Tw"

  if (length(times_elasticities) > 0) {
    ti_opt <- get_ti_thph(times_elasticities, work_elasticities$Theta, tw_opt, tau, Tc)
    opt <- cbind(opt, ti_opt)
    nms <- c(nms, paste0("Tf", seq_along(times_elasticities)))
  }
  if (length(goods_elasticities) > 0) {
    xj_opt <- get_xi_thph(goods_elasticities, 1, work_elasticities$Phi, tw_opt, Ec, w)
    opt <- cbind(opt, xj_opt)
    nms <- c(nms, paste0("Ef", seq_along(goods_elasticities)))
  }
  colnames(opt) <- nms

  M <- ncol(opt)
  obs <- opt + matrix(stats::rnorm(n * M, sd = sd), ncol = M)
  colnames(obs) <- nms

  database <- data.frame(PeID = 1:n, Tc = Tc, EcI = Ec, w = w)
  database <- cbind(database, as.data.frame(obs))

  mtuem_settings <- list(
    work_times = "Tw",
    free_times = if (length(times_elasticities) > 0) nms[2:(1 + length(times_elasticities))] else list(),
    free_goods = if (length(goods_elasticities) > 0) utils::tail(nms, length(goods_elasticities)) else list(),
    work_elasticities = work_elasticities,
    times_elasticities = times_elasticities,
    goods_elasticities = goods_elasticities,
    Tc = Tc, Ec = Ec, w = w, tau = tau
  )

  list(database = database,
       apollo_inputs = list(database = database),
       mtuem_settings = mtuem_settings,
       opt = opt, obs = obs,
       Tc = Tc, Ec = Ec, w = w, tau = tau, n = n)
}

# Reference zero-mean multivariate normal density, computed with solve() and
# determinant() so it is independent of the Cholesky/conditional-decomposition
# machinery under test.
mtuem_dmvnorm_ref <- function(E, Sigma) {
  E <- as.matrix(E)
  M <- ncol(E)
  Sinv <- solve(Sigma)
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  q <- rowSums((E %*% Sinv) * E)
  exp(-0.5 * (q + logdet + M * log(2 * base::pi)))
}

# Lower-triangular Cholesky factor of Sigma, flattened row-major
# (L11, L21, L22, L31, L32, L33, ...) exactly as the likelihoods expect.
mtuem_chol_vec <- function(Sigma) {
  M <- ncol(Sigma)
  L <- t(chol(Sigma))
  rc <- do.call(rbind, lapply(1:M, function(r) cbind(r, 1:r)))
  as.numeric(L[rc])
}

# Minimal object that passes apollo_deltaMethod's model validity check.
mtuem_fake_apollo_model <- function(estimate, varcov) {
  dimnames(varcov) <- list(names(estimate), names(estimate))
  list(estimate = estimate,
       varcov = varcov,
       robvarcov = varcov,
       apollo_fixed = character(0),
       apollo_control = list(modelName = "fake", noValidation = TRUE,
                             HB = FALSE, mixing = FALSE, panelData = FALSE))
}
