#library(testthat)
test_that("Works with objects returned by ptFCReg", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 0.5
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  EYgivenX <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2]
  }))
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  res <- ptFCReg(tGrid = tGrid, dat = dat)
  smres <- smPtFCRegCoef(res, bw = 2.5 / (nGridIn-1), kernel_type = 'epan')
  fit_res <- fitted_ptFCReg(res)
  expect_lt( max(sqrt(apply( (fit_res - EYgivenX)^2, 1, pracma::trapz, x = tGrid))), 0.15)
})

test_that("Works with objects returned by smPtFCRegCoef", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 0.5
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  EYgivenX <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2]
  }))
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  res <- ptFCReg(tGrid = tGrid, dat = dat)
  smres <- smPtFCRegCoef(res, bw = 2.5 / (nGridIn-1), kernel_type = 'epan')
  fit_smres <- fitted_ptFCReg(smres)
  expect_lt( max(sqrt(apply( (fit_smres - EYgivenX)^2, 1, pracma::trapz, x = tGrid))), 0.15)
})
