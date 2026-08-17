library(testthat)
tryCatch(devtools::load_all(".", quiet = TRUE), error = function(e) library(mtuem))

test_that("get_tw_additive_quadratic computes exact analytical work time", {
  # Arrange
  tau <- 168
  Tc <- 80
  Ec <- 200
  w <- 15
  
  theta_w <- 0
  beta_w  <- 10
  
  theta_1 <- 800; beta_1 <- 20
  theta_i <- 3000; beta_i <- 300
  
  phi_1 <- 100; eta_1 <- 1
  phi_i <- 80;  eta_i <- 0.5
  
  times_elasticities <- list(theta_1 = theta_1, theta_i = theta_i)
  times_satiations   <- list(beta_1 = beta_1, beta_i = beta_i)
  goods_elasticities <- list(phi_1 = phi_1, phi_i = phi_i)
  goods_satiations   <- list(eta_1 = eta_1, eta_i = eta_i)
  goods_cost         <- list(E1 = 1, Ei = 1)
  
  # Act
  tw_calc <- get_tw_additive_quadratic(
    work_elasticity = theta_w,
    work_satiation  = beta_w,
    times_elasticities = times_elasticities,
    times_satiations   = times_satiations,
    goods_elasticities = goods_elasticities,
    goods_satiations   = goods_satiations,
    goods_cost         = goods_cost,
    Tc = Tc, Ec = Ec, w = w, tau = tau
  )
  
  # Manual analytical formula
  Sth <- theta_1 / beta_1 + theta_i / beta_i
  B   <- 1 / beta_1 + 1 / beta_i
  Sph <- (1 * phi_1 / eta_1) + (1 * phi_i / eta_i)
  H   <- (1^2 / eta_1) + (1^2 / eta_i)
  
  tw_expected <- (theta_w + w * (Sph + Ec) / H - (Sth - (tau - Tc)) / B) / (beta_w + w^2 / H + 1 / B)
  
  # Assert
  expect_equal(tw_calc, tw_expected, tolerance = 1e-9)
  expect_true(tw_calc > 0)
  expect_true(tw_calc < (tau - Tc))
})

test_that("Optimal times and expenditures satisfy exact budget constraints", {
  # Arrange
  tau <- 168
  Tc <- 70
  Ec <- 150
  w <- 12
  
  theta_w <- 0; beta_w <- 8
  times_elasticities <- list(theta_1 = 600, theta_i = 2500)
  times_satiations   <- list(beta_1 = 20,  beta_i = 250)
  goods_elasticities <- list(phi_1 = 120, phi_i = 70)
  goods_satiations   <- list(eta_1 = 1,   eta_i = 0.4)
  goods_cost         <- list(E1 = 1, Ei = 1)
  
  # Act
  tw <- get_tw_additive_quadratic(
    work_elasticity = theta_w, work_satiation = beta_w,
    times_elasticities = times_elasticities, times_satiations = times_satiations,
    goods_elasticities = goods_elasticities, goods_satiations = goods_satiations,
    goods_cost = goods_cost, Tc = Tc, Ec = Ec, w = w, tau = tau
  )
  
  ti <- get_ti_additive_quadratic(
    Tw = tw, times_elasticities = times_elasticities,
    times_satiations = times_satiations, Tc = Tc, tau = tau
  )
  
  xi <- get_xi_additive_quadratic(
    Tw = tw, goods_elasticities = goods_elasticities,
    goods_satiations = goods_satiations, goods_cost = goods_cost,
    Ec = Ec, w = w
  )
  
  # Assert Adding-up Identities
  time_total <- tw + ti$T1 + ti$Ti
  money_total <- xi$E1 + xi$Ei
  
  expect_equal(time_total, tau - Tc, tolerance = 1e-9)
  expect_equal(money_total, w * tw - Ec, tolerance = 1e-9)
})

test_that("Values of time satisfy fundamental identity VoL = w + VTAW", {
  # Arrange
  tau <- 168
  Tc <- 75
  Ec <- 180
  w <- 14
  
  theta_w <- 0; beta_w <- 11
  times_elasticities <- list(theta_1 = 700, theta_i = 3000)
  times_satiations   <- list(beta_1 = 22,  beta_i = 280)
  goods_elasticities <- list(phi_1 = 110, phi_i = 75)
  goods_satiations   <- list(eta_1 = 1,   eta_i = 0.45)
  goods_cost         <- list(E1 = 1, Ei = 1)
  
  tw <- get_tw_additive_quadratic(
    work_elasticity = theta_w, work_satiation = beta_w,
    times_elasticities = times_elasticities, times_satiations = times_satiations,
    goods_elasticities = goods_elasticities, goods_satiations = goods_satiations,
    goods_cost = goods_cost, Tc = Tc, Ec = Ec, w = w, tau = tau
  )
  
  # Act
  vot <- get_values_of_time_additive_quadratic(
    Tw = tw, work_elasticity = theta_w, work_satiation = beta_w,
    times_elasticities = times_elasticities, times_satiations = times_satiations,
    goods_elasticities = goods_elasticities, goods_satiations = goods_satiations,
    goods_cost = goods_cost, Tc = Tc, Ec = Ec, w = w, tau = tau
  )
  
  # Assert
  expect_true(!is.na(vot$VoL))
  expect_true(!is.na(vot$VTAW))
  expect_equal(vot$VoL, w + vot$VTAW, tolerance = 1e-9)
})
