

#library(fdapace)
#library(testthat)

test_that("Works with objects returned by ConcurReg", {
  set.seed(1)
  n <- 75
  nGridIn <- 150
  sparsity <- 5:10 # Sparse data sparsity
  Sup <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
  bw <- 0.1
  outGrid <- round(seq(min(Sup), 1, by=0.05), 2)
  # Simulate functional data
  mu <- Sup * 2 # mean function for X_1
  sigma <- 1
  beta_0 <- 0
  beta_1 <- 1
  beta_2 <- 1
  Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
  X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
  epsilon <- rnorm(n, sd=sigma)
  Y <- matrix(NA, n, nGridIn)
  for (i in seq_len(n)) {
    Y[i, ] <- beta_0 + beta_1 * X_1[i, ] + beta_2 * Z[i, 2] + epsilon[i]
  }
  # Sparsify functional data
  X_1sp <- fdapace::Sparsify(X_1, Sup, sparsity)
  Ysp <- fdapace::Sparsify(Y, Sup, sparsity)
  vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
  withError2D <- ConcurReg(vars, outGrid)
  l1 <- sqrt(mean((withError2D$beta0)^2))
  l2 <- sqrt(mean((withError2D$beta[1,]-1)^2))
  l3 <- sqrt(mean((withError2D$beta[2,]-1)^2))
  expect_lt( max(c(l1, l2, l3)), 1.5)
})








