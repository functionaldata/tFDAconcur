library(testthat)
test_that('error: please provide Lag with appropriate value', {
  set.seed(1)
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 50
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
  
  Lag <- c(3)
  t0_out <- t0[t0>=max(Lag)]
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2
  ### functional response Y(t), t in t0_out ### 
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- beta_0(t0_out) +  denseX0$Ly[[i]][1:length(t0_out)]* beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1) 
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)
  expect_error(ConcReg_Lag(Y=denseY, X=list(X1 = denseX0)), 
               'please provide Lag with appropriate value')
})

test_that('error: Lag should be smaller than the length of the time interval', {
  set.seed(1)
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 50
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
  
  Lag <- c(3)
  t0_out <- t0[t0>=max(Lag)]
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2
  ### functional response Y(t), t in t0_out ### 
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- beta_0(t0_out) +  denseX0$Ly[[i]][1:length(t0_out)]* beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1) 
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)
  expect_error(ConcReg_Lag(Y=denseY, X=list(X1 = denseX0), Lag=30),
               'Lag should be smaller than the length of the time interval')
})

test_that('functional concurrent regression model with lag simulated setting works (accurate estimate to the true target)', {
  set.seed(1)
  phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
  phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
  lambdaX <- c(10, 5)
  n <- 50
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
  
  Lag <- c(3)
  t0_out <- t0[t0>=max(Lag)]
  beta_1 <- function(t) 5*sin(pi*t/10)
  beta_0 <- function(t) t^2/2
  ### functional response Y(t), t in t0_out ### 
  denseLt <- list()
  denseLy <- list()
  for (i in 1:n) {
    denseLt[[i]] <- t0_out
    denseLy[[i]] <- beta_0(t0_out) +  denseX0$Ly[[i]][1:length(t0_out)]* beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1) 
  }
  denseY <- list(Ly = denseLy, Lt = denseLt)
  model = ConcReg_Lag(Y=denseY, X=list(X1 = denseX0), Lag=Lag)
  
  eps1 = sum((model$beta-beta_1(model$workGridY))^2)
  eps0 = sum((model$beta0-beta_0(model$workGridY))^2)
  
  expect_equal(eps1<=2e-3 & eps0<=2e-2 , TRUE)
})
