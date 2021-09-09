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
  denseX0 <- list(Ly = denseLy,Lt = denseLt)
  
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
  
  ### sparify the dense data ###
  SparseY <- fdapace::Sparsify(samp = t(sapply(denseY$Ly, function(v) v)), 
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE)
  )
  SparseX <- fdapace::Sparsify(samp = t(sapply(denseX$Ly, function(v) v)),
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE))
  
  ### model fitting for sparse case ###
  expect_error(historyIndexSparse(Y = SparseY, X = list(X = SparseX)), 
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
  denseX0 <- list(Ly = denseLy,Lt = denseLt)
  
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
  
  ### sparify the dense data ###
  SparseY <- fdapace::Sparsify(samp = t(sapply(denseY$Ly, function(v) v)), 
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE)
  )
  SparseX <- fdapace::Sparsify(samp = t(sapply(denseX$Ly, function(v) v)),
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE))
  
  ### model fitting for sparse case ###
  expect_error(historyIndexSparse(Y = SparseY, X = list(X = SparseX), Lag = 30), 
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
  denseX0 <- list(Ly = denseLy,Lt = denseLt)
  
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
  
  ### sparify the dense data ###
  SparseY <- fdapace::Sparsify(samp = t(sapply(denseY$Ly, function(v) v)), 
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE)
  )
  SparseX <- fdapace::Sparsify(samp = t(sapply(denseX$Ly, function(v) v)),
                               pts = t0_out,
                               sparsity = sample(10:20, n, replace = TRUE))
  
  ### model fitting for sparse case ###
  Fit_result <- historyIndexSparse(Y = SparseY, X = list(X = SparseX), Lag = Lag)
  eps1=sum((beta_1(Fit_result$workGridY) - Fit_result$beta1[[1]])^2) / sum((beta_1(Fit_result$workGridY))^2)
  eps2=sum((beta_0(Fit_result$workGridY) - Fit_result$beta0)^2) / sum((beta_0(Fit_result$workGridY))^2)
  eps3=sum((gamma_u(Fit_result$workGridLag[[1]]) - Fit_result$gamma[[1]])^2) / sum((gamma_u(Fit_result$workGridLag[[1]]))^2)
  
  expect_equal(eps1<=1e-1 & eps2<=1e-1 & eps3<=1e-1, TRUE)
})
