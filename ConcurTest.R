
#usethis::create_package("mypackage")
#RcppArmadillo.package.skeleton("mypackage")
#devtools::load_all("tFDAconcur-FCReg/")
#devtools::document("tPACETemp/")
#install.packages("tPACETemp",repos=NULL,type="source",dependencies = TRUE)
#devtools::document("mypackage")
#library(RcppArmadillo)
#library(mypackage)

library(fdapace)

set.seed(1)
n <- 75
nGridIn <- 150
sparsity <- 5:10 # Sparse data sparsity
T <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
bw <- 0.1
outGrid <- round(seq(min(T), 1, by=0.05), 2)
# Simulate functional data
mu <- T * 2 # mean function for X_1
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
set.seed(1)
X_1sp <- Sparsify(X_1, T, sparsity)
set.seed(1)
Ysp <- Sparsify(Y, T, sparsity)
vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
withError2D <- ConcurReg(vars, outGrid)









