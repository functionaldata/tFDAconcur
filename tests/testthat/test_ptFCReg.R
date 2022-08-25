# library(testthat)
test_that("The field in dat for the response is not a matrix.", {
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  Z1 <- tGrid
  beta0 <- 0
  dat <- list(Y = beta0 + Z1 + rnorm(nGridIn), Z1=Z1)
  expect_error(ptFCReg(tGrid = tGrid, dat = dat),
               "The field in dat for the response is not a matrix.")
})
test_that("Fields of dat are not vectors or matrices.", {
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  Z1 <- seq(0,1,length.out = n)
  beta0 <- 0
  beta1 <- cos(tGrid)
  sigma <- 1
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta1 * Z1[i] + epsilon[i]
  }))
  dat <- list(Y = list(t = tGrid, Y = Y), Z1=Z1)
  expect_error(ptFCReg(tGrid = tGrid, dat = dat),
               "The field in dat for the response is not a matrix.")
})
test_that("Columns in the matrices corresponding to the functional variables do not match with tGrid.", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  X_1 <- X_1[,-1]
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  expect_error(
    ptFCReg(tGrid = tGrid, dat = dat),
    "Columns in the matrices corresponding to the functional variables do not match with tGrid."
  )
})
test_that("Numbers of rows in the matrices corresponding to the functional variables are not all the same.", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  X_1 <- X_1[-1,]
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  expect_error(
    ptFCReg(tGrid = tGrid, dat = dat),
    "Numbers of rows in the matrices corresponding to the functional variables are not all the same."
  )
})
test_that("Lengths of the vectors corresponding to scalar covariates are not all the same as the number of subjects in the functional response.", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  Z <- Z[-1,]
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  expect_error(
    ptFCReg(tGrid = tGrid, dat = dat),
    "Lengths of the vectors corresponding to scalar covariates are not all the same as the number of subjects in the functional response."
  )
})
test_that("Works for functional covariates alone", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- cos(tGrid)
  Z <- rnorm(n, 0, sqrt(2))
  X_1 <- Z %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta * X_1[i, ]+ epsilon[i]
  }))
  dat <- list(X1=X_1,  Y=Y)
  res <- ptFCReg(tGrid = tGrid, dat = dat)
  expect_lt( max(sqrt(apply( (res$beta - beta)^2, 1, pracma::trapz, x = tGrid))), sigma/20)
})
test_that("Works for scalar covariates alone", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- 1.5 + sin(tGrid)
  Z <- rnorm(n, 0, sqrt(2))
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta * Z[i] + epsilon[i]
  }))
  dat <- list(Z1=Z, Y=Y)
  res <- ptFCReg(tGrid = tGrid, dat = dat)
  expect_lt( max(sqrt(apply( (res$beta - beta)^2, 1, pracma::trapz, x = tGrid))), sigma/20)
})
test_that("Works for functional and scalar covariates", {
  set.seed(1)
  n <- 50
  nGridIn <- 101
  tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
  muX1 <- tGrid * 2 # mean function for X_1
  sigma <- 1
  beta0 <- 0
  beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- t(sapply(seq_len(n), function(i) {
    beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
  }))
  dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
  res <- ptFCReg(tGrid = tGrid, dat = dat)
  expect_lt( max(sqrt(apply( (res$beta - beta)^2, 1, pracma::trapz, x = tGrid))), sigma/20)
})

