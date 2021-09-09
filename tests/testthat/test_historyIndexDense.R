library(testthat)
test_that('error: please provide Lag with appropriate value', {
  set.seed(1)
  ### functional covariate X(t) ###
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 150
  N <- 101
  Xi <- matrix(rnorm(2*n), nrow = n, ncol = 2)
  denseLt <- list()
  denseLy <- list()
  t0 <- seq(0, 15, length.out = N)
  for (i in 1:n) {
    denseLt[[i]] <- t0
    denseLy[[i]] <- lambdaX[1]*Xi[i, 1]*phi1(t0) + lambdaX[2]*Xi[i, 2]*phi2(t0)
  }
  denseX0 <- list(Ly = denseLy, Lt = denseLt)

  ### generate coefficient function gamma(u), beta(u) ###
  Lag <- 5
  u0 <- t0[t0<=Lag]
  t0_out <- t0[t0>=Lag]
  gamma_u <- function(u) sqrt(2/5) * cos(pi * u /5)
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2

  ### functional response Y(t), t in t0_out ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    Xt <- denseX0$Ly[[i]]
    Xtu <- t(sapply((1:N)[t0>=Lag], function(j){
      rev(Xt[(j-length(u0)+1):j])  #history index for X[t-u:t]
    }))
    IntGammaXtu <- apply(Xtu, 1, function(v){
      fdapace::trapzRcpp(u0, gamma_u(u0) * v)
    })
    #append 0 in the first length(u0)-1 element(useless info. in our modeling)
    denseLy[[i]] <- beta_0(t0_out) + IntGammaXtu * beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1)
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)

  ### functional predictor X(t) (adjust for t0_out) ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n){
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- denseX0$Ly[[i]][t0>=Lag]
  }
  denseX <- list(Ly = denseLy,
                 Lt = denseLt)
  expect_error(historyIndexDense(Y = denseY, X = list(X = denseX)), 
               'please provide Lag with appropriate value')
})

test_that('error: Lag should be smaller than the length of the time interval', {
  set.seed(1)
  ### functional covariate X(t) ###
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 150
  N <- 101
  Xi <- matrix(rnorm(2*n), nrow = n, ncol = 2)
  denseLt <- list()
  denseLy <- list()
  t0 <- seq(0, 15, length.out = N)
  for (i in 1:n) {
    denseLt[[i]] <- t0
    denseLy[[i]] <- lambdaX[1]*Xi[i, 1]*phi1(t0) + lambdaX[2]*Xi[i, 2]*phi2(t0)
  }
  denseX0 <- list(Ly = denseLy, Lt = denseLt)
  
  ### generate coefficient function gamma(u), beta(u) ###
  Lag <- 5
  u0 <- t0[t0<=Lag]
  t0_out <- t0[t0>=Lag]
  gamma_u <- function(u) sqrt(2/5) * cos(pi * u /5)
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2
  
  ### functional response Y(t), t in t0_out ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    Xt <- denseX0$Ly[[i]]
    Xtu <- t(sapply((1:N)[t0>=Lag], function(j){
      rev(Xt[(j-length(u0)+1):j])  #history index for X[t-u:t]
    }))
    IntGammaXtu <- apply(Xtu, 1, function(v){
      fdapace::trapzRcpp(u0, gamma_u(u0) * v)
    })
    #append 0 in the first length(u0)-1 element(useless info. in our modeling)
    denseLy[[i]] <- beta_0(t0_out) + IntGammaXtu * beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1)
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)
  
  ### functional predictor X(t) (adjust for t0_out) ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n){
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- denseX0$Ly[[i]][t0>=Lag]
  }
  denseX <- list(Ly = denseLy,
                 Lt = denseLt)
  expect_error(historyIndexDense(Y = denseY, X = list(X = denseX), Lag = 15), 
               'Lag should be smaller than the length of the time interval')
})

test_that('functional history index model simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  ### functional covariate X(t) ###
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 150
  N <- 101
  Xi <- matrix(rnorm(2*n), nrow = n, ncol = 2)
  denseLt <- list()
  denseLy <- list()
  t0 <- seq(0, 15, length.out = N)
  for (i in 1:n) {
    denseLt[[i]] <- t0
    denseLy[[i]] <- lambdaX[1]*Xi[i, 1]*phi1(t0) + lambdaX[2]*Xi[i, 2]*phi2(t0)
  }
  denseX0 <- list(Ly = denseLy, Lt = denseLt)
  
  ### generate coefficient function gamma(u), beta(u) ###
  Lag <- 5
  u0 <- t0[t0<=Lag]
  t0_out <- t0[t0>=Lag]
  gamma_u <- function(u) sqrt(2/5) * cos(pi * u /5)
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2
  
  ### functional response Y(t), t in t0_out ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    Xt <- denseX0$Ly[[i]]
    Xtu <- t(sapply((1:N)[t0>=Lag], function(j){
      rev(Xt[(j-length(u0)+1):j])  #history index for X[t-u:t]
    }))
    IntGammaXtu <- apply(Xtu, 1, function(v){
      fdapace::trapzRcpp(u0, gamma_u(u0) * v)
    })
    #append 0 in the first length(u0)-1 element(useless info. in our modeling)
    denseLy[[i]] <- beta_0(t0_out) + IntGammaXtu * beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1)
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)
  
  ### functional predictor X(t) (adjust for t0_out) ###
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n){
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- denseX0$Ly[[i]][t0>=Lag]
  }
  denseX <- list(Ly = denseLy,
                 Lt = denseLt)
  fit <- historyIndexDense(Y = denseY, X = list(X = denseX), Lag = Lag)
  eps <- sum((fit$beta0 - beta_0(fit$workGridY))^2)/sum((beta_0(fit$workGridY))^2)
  expect_equal(eps<=1e-3, TRUE)
})
