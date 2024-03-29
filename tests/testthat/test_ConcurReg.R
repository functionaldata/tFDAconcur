#library(testthat)
# devtools::load_all()

test_that('Simple dense case works', {
  # X(t) = Y(t), and they are contaminated with errors.
  # beta_0 = 0; beta_1 = t
  set.seed(1)
  n <- 100
  p <- 20
  sigma <- 1
  pts <- seq(0, 1, length.out=p)
  T <- 1 + pts * p
  X <- fdapace::Wiener(n, pts) + rnorm(n)
  Y <- t(apply(X, 1, `*`, e2=pts))
  Xn <- X + rnorm(n * p, sd=sigma)
  Yn <- Y + rnorm(n * p, sd=sigma)
  vars <- list(X=fdapace::MakeFPCAInputs(tVec=T, yVec=Xn), 
               Y=fdapace::MakeFPCAInputs(tVec=T, yVec=Yn))
  
  bw <- 0.25 * diff(range(T))
  kern <- 'epan'
  res = ConcurReg(vars, T, measurementError=FALSE, diag1D='none')
  resN = ConcurReg(vars, T, measurementError=TRUE, diag1D='none')
  # plot(pts, res$beta); abline(a=0, b=1)
  # plot(pts, resN$beta); abline(a=0, b=1)
  res$beta = as.numeric(res$beta)
  resN$beta = as.numeric(resN$beta)
  pos = !is.na(res$beta)
  posN = !is.na(resN$beta)
  res$beta = res$beta[complete.cases(res$beta)]
  resN$beta = resN$beta[complete.cases(resN$beta)]
  expect_equal(res$beta, pts[pos], scale=1, tolerance=0.1)
  expect_equal(resN$beta, pts[posN], scale=1, tolerance=0.1)
  # Assume noise is better than no noise
  expect_lt(mean((resN$beta - pts[posN])^2), mean((res$beta - pts[pos])^2)) 
})

# Y(t) = \beta_0(t) + \beta_1(t) X_1(t) + \beta_2(t) X_2(t) + \beta_3 Z_3 + \beta_4 Z_4 + \epsilon
# X_1(t) = \mu_{X_1}(t) + Z_1 \phi_1(t); X_2(t) = \mu_{X_2}(t) + Z_2 \phi_2(t);
# \mu_{X_1}(t) = 2t, \mu_{X_2}(t) = -2t 
# (Z_1, \dots, Z_4) \sim N(0, \Sigma). \Sigma = 1/4 + diag(3/4)
# \phi_1(t) = 1, \phi_2(t) = sqrt(2) * cos(\pi t), t \in [0, 1]. So $cov(X_1(s), X_2(t)) = sqrt(2)/4 \cos(\pi t).$
# \epsilon \sim N(0, \sigma^2).
# \beta_1 = \beta_2 = \beta_3 = \beta_4 = 1.

set.seed(1)
n <- 100
nGridIn <- 200
sparsity <- 5:10 # must have length > 1
bw <- 0.2
kern <- 'gauss'
T <- seq(0, 1, length.out=nGridIn)

Sigma <- 1 / 4 + diag(3 / 4, 4)
mu <- T * 2
sigma <- 1

beta_0 <- 0
beta_1 <- 1
beta_2 <- 1
beta_3 <- 1
beta_4 <- 1

Z <- MASS::mvrnorm(n, rep(0, 4), Sigma)
X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
X_2 <- Z[, 2, drop=FALSE] %*% matrix(sqrt(2) * cos(pi * T), 1, nGridIn) - matrix(mu, n, nGridIn, byrow=TRUE)
# tmp <- cov(X_1, X_2); colMeans(tmp)
epsilon <- rnorm(n, sd=sigma)
Y <- matrix(NA, n, nGridIn)
for (i in seq_len(n)) {
  Y[i, ] <- beta_0 + beta_1 * X_1[i, ] + beta_2 * X_2[i, ] + beta_3 * Z[i, 3] + beta_4 * Z[i, 4] + epsilon[i]
}
denseMuY <- colMeans(Y)
denseCovY <- cov(Y)

set.seed(1)
X_1sp <- fdapace::Sparsify(X_1, T, sparsity)
set.seed(1)
X_2sp <- fdapace::Sparsify(X_2, T, sparsity)
set.seed(1)
Ysp <- fdapace::Sparsify(Y, T, sparsity)
outGrid <- round(seq(min(T), 1, by=0.05), 2)
vars <- list(X_1=X_1sp, X_2=X_2sp, Z_3=Z[, 3], Z_4=Z[, 4], Y=Ysp)

test_that('Scalar-scalar cov', 
          expect_equal(as.numeric(uniCov(Z[, 3], Z[, 4])), cov(Z[, 3], Z[, 4]))
)

test_that('Scalar-function cov = Function-scalar cov', {
  expect_equal(uniCov(X_1sp, Z[, 3], bw, outGrid), 
               t(uniCov(Z[, 3], X_1sp, bw, outGrid)))
})

cov12 <- uniCov(X_1sp, X_2sp, rep(bw,2), outGrid, kern)
# cov12_rd <- uniCov(X_1sp, X_2sp, bw, outGrid, kern, rmDiag=TRUE)
cov12_1D <- uniCov(X_1sp, X_2sp, rep(bw,2), outGrid, kern, use1D=TRUE)
cov13 <- uniCov(X_1sp, Z[, 3], rep(bw,2), outGrid, kern)
cov21 <- uniCov(X_2sp, X_1sp, rep(bw,2), outGrid, kern)
# cov11 <- uniCov(X_1sp, X_1sp, bw, outGrid, kern)
# cov11_rd <- uniCov(X_1sp, X_1sp, bw, outGrid, kern, rmDiag=TRUE)
# cov22 <- uniCov(X_2sp, X_2sp, bw, outGrid, kern)
# cov22_rd <- uniCov(X_2sp, X_2sp, bw, outGrid, kern, rmDiag=TRUE)
# rgl::persp3d(outGrid, outGrid, cov12)
test_that('Function-function cov works', {
  expect_equal(cov12, t(cov21))
  
  # x-direction is close to constant.
  expect_true(mean(apply(cov12, 2, sd), trim=0.1) < 0.08)
  
  # y-direction is close to 1/4 * (1.5 + cos(pi t))
  expect_true(sqrt(mean((colMeans(cov12) - sqrt(2) / 4 * cos(pi * outGrid))^2, trim=0.1)) < 0.133)
  
  # 1D and 2D smoother is similar
  expect_equal(diag(cov12), diag(cov12_1D), tolerance=0.2)
})
covAll <- MvCov(vars, rep(bw,5), outGrid, kern,center= TRUE)
covAllNoError <- MvCov(vars, rep(bw,5), outGrid, kern, center= TRUE,measurementError=FALSE)

test_that('Multi-function/scalar cov works', {
  # The cov(x, y) and cov(y, x) is symmetric.
  expect_equal(as.numeric(cov12), as.numeric(covAll[, , 1, 2]))
  expect_equal(covAll[, , 3, 1], t(covAll[, , 1, 3]))
  expect_equal(covAll[, , 4, 3], t(covAll[, , 3, 4]))
  expect_equal(covAll[, , 1, 2], t(covAll[, , 2, 1]))
  expect_equal(covAll[, , 1, 2], covAllNoError[, , 1, 2])
  tmp <- max(abs(covAll[, , 1, 1] - covAllNoError[, , 1, 1]))
  expect_true(tmp > 0.01 && tmp < 0.35)
  
  # rgl::persp3d(outGrid, outGrid, covAll[, , 1, 1], col='blue', xlab='X(s)', ylab='Y(t)')
  # rgl::persp3d(outGrid, outGrid, covAll[, , 1, 1], col='blue', xlab='X(s)', ylab='Y(t)')
})

demeanRes <- demean(vars, rep(bw,5), kern)
varsDemean <- demeanRes[['xList']]
muDemean <- demeanRes[['muList']]
covAllDemean <- MvCov(varsDemean, rep(bw,5), outGrid, kern, center=FALSE)

test_that('demean works', {
  expect_equal(demeanRes[['muList']][['Z_3']], mean(vars[['Z_3']]))
  expect_equal(demeanRes[['muList']][['Z_4']], mean(vars[['Z_4']]))
  expect_equal(uniCov(varsDemean[['X_1']], varsDemean[['X_2']], rep(bw,2), outGrid, center=FALSE), cov12)
  expect_equal(uniCov(varsDemean[['X_1']], varsDemean[['Z_3']], rep(bw,2), outGrid, center=FALSE), cov13)
  expect_equal(covAll, covAllDemean)
})

# withError2D <- ConcurReg(vars,outGrid, bw,bw, measurementError = TRUE)
# withError1D <- ConcurReg(vars, outGrid, bw,bw, measurementError = TRUE, diag1D='cross')
noError2D <- ConcurReg(vars, outGrid, rep(bw,5), rep(bw,5), measurementError=FALSE, diag1D = "none")
noError1D <- ConcurReg(vars, outGrid, bw, bw, measurementError=FALSE, diag1D = 'all')

# matplot(outGrid, t(withError2D$beta), 'l')
# matplot(outGrid, t(noError2D$beta), 'l')
# matplot(outGrid, t(noError1D$beta), 'l')

expect_error(ConcurReg(vars, outGrid, bw, bw,  measurementError=TRUE, diag1D='all'), "Cannot assume measurement error when diag1D == 'all'")

# # Minimal eigenvalues sometimes smaller than 0.
# minLambda <- sapply(seq_along(outGrid), function(i) {
# min(eigen(noError1D[['cov']][i, i, 1:4, 1:4])[['values']])
# })
# minLambda2D <- sapply(seq_along(outGrid), function(i) {
# min(eigen(noError2D[['cov']][i, i, 1:4, 1:4])[['values']])
# })

# since we drop diag1D = "cross" case, we do not need to test it.
test_that('1D and 2D covariance estimates are similar', {
  # expect_true(sqrt(mean(
  #   (withError2D[['beta']][,-which(colSums(is.na(withError2D[['beta']]))>0)] -
  #      withError1D[['beta']][,-which(colSums(is.na(withError1D[['beta']]))>0)]
  #   )^2,
  #   trim=0.2)) < 0.2)
  expect_equal(noError2D[['beta']], noError1D[['beta']], tolerance=0.1)
})

# withError2DRect <- ConcurReg(vars, outGrid,bw, bw,  kern='rect')
# withError1DRect <- ConcurReg(vars,  outGrid,bw,bw, diag1D='cross', kern='rect')
noError2DRect <- ConcurReg(vars, outGrid, rep(bw,5), rep(bw,5), measurementError=FALSE, diag1D = "none", kern='rect')
noError1DRect <- ConcurReg(vars,  outGrid, bw, bw, measurementError=FALSE, diag1D='all', kern='rect')
# withError2DEpan <- ConcurReg(vars,  outGrid,bw,bw,  kern='epan')
# withError1DEpan <- ConcurReg(vars, outGrid, bw,bw,diag1D='cross', kern='epan')
noError2DEpan <- ConcurReg(vars, outGrid, rep(bw,5), rep(bw,5),  measurementError=FALSE, diag1D = "none", kern='epan')
noError1DEpan <- ConcurReg(vars, outGrid, bw, bw, measurementError=FALSE, diag1D='all', kern='epan')

test_that('Different kernel type works', {
  # expect_true(sqrt(mean(
  #   (withError2DRect[['beta']][,-which(colSums(is.na(withError2DRect[['beta']]))>0)] - 
  #      withError1DRect[['beta']][,-which(colSums(is.na(withError1DRect[['beta']]))>0)])^2, 
  #   trim=0.2)) < 0.2)
  expect_equal(noError2DRect[['beta']], noError2D[['beta']], tolerance=0.2)
  expect_equal(noError1DRect[['beta']], noError1D[['beta']], tolerance=0.2)
  # expect_equal(withError2DRect[['beta']], withError2DEpan[['beta']], tolerance=0.2)
  # expect_equal(withError1DRect[['beta']], withError1DEpan[['beta']], tolerance=0.2)
  expect_equal(noError2DRect[['beta']], noError2DEpan[['beta']], tolerance=0.2)
  expect_equal(noError1DRect[['beta']], noError1DEpan[['beta']], tolerance=0.2)
})

fpcaX1 <- fdapace::FPCA(X_1sp[['Ly']], X_1sp[['Lt']], list(userBwMu=bw, userBwCov=bw))
fpcaX2 <- fdapace::FPCA(X_2sp[['Ly']], X_2sp[['Lt']], list(userBwMu=bw, userBwCov=bw))
fpcaY <- fdapace::FPCA(Ysp[['Ly']], Ysp[['Lt']], list(userBwMu=bw, userBwCov=bw))
FPCAlist <- list(Y=fpcaY, X_1=fpcaX1, X_2=fpcaX2)
# imputeRes <- imputeConReg(FPCAlist, Z[, 3:4], outGrid)
# test_that('imputation and 1D beta estimates are similar', {
#   expect_equal(imputeRes[['beta0']], noError1D[['beta0']], scale=1, tolerance=0.5)
#   expect_equal(unname(imputeRes[['beta']]), unname(noError1D[['beta']]), scale=1, tolerance=0.5)
# })

## subsetVars
test_that('subseting covariates is fine', {
  subVars <- subsetVars(vars, 1)
  subVarsTF <- subsetVars(vars, c(TRUE, rep(FALSE, n - 1)))
  expect_equal(subVars, subVarsTF)
  expect_equal(subVars[['X_1']][['Lt']][[1]], X_1sp[['Lt']][[1]])
  expect_equal(length(subVars[['X_1']][['Lt']]), 1)
  expect_equal(length(subVars[['Z_3']]), 1)
})

test_that('Test based on previous implementation: simple concurrent regression works fine', {
  set.seed(123);  N = 1001;  M = 101;
  # Define the continuum
  s = seq(0,10,length.out = M)
  # Define the mean and 2 eigencomponents
  meanFunct <- function(s) 0.2*s + 2.0*exp(-(s-5)^2)
  eigFunct1 <- function(s) +cos(2*s*pi/10) / sqrt(5)
  eigFunct2 <- function(s) 1* -sin(2*s*pi/10) / sqrt(5)
  # Create FPC scores
  Ksi = matrix(rnorm(N*2), ncol=2);
  Ksi = apply(Ksi, 2, scale)
  Ksi = Ksi %*% diag(c(5,2))
  # Create X_covariate
  xTrue = Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2))
  # Create beta_Func
  betaFunc1 = c(2,2) %*%  t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2))
  z1 <- rnorm(N,sd=1)
  # Create scalar dep. variable a
  y = matrix(rep(0,N*M), ncol = M); 
  yTrue = matrix(rep(0,N*M), ncol = M); 
  for (i in 1:N) { 
    y[i,] = rnorm(sd=1.99,0, n=M) + z1[i] * 2.5 +  0.2 * s +(xTrue[i,]) * t(betaFunc1); 
    yTrue[i,] = rnorm(sd=0.0011,0, n=M) + z1[i] * 2.5 +  0.2 * s +(xTrue[i,]) * t(betaFunc1); 
  }
  sparsitySchedule = 1:16;
  set.seed(1)
  Yf <- fdapace::Sparsify(y, s, sparsitySchedule)
  set.seed(1)
  Xf <- fdapace::Sparsify(xTrue, s, sparsitySchedule)
  
  outGrid <- s
  vars <- list(X = Xf, Z = z1, Y = Yf)
  
  Q <- ConcurReg(vars, outGrid, 0.5,0.5, kern = 'epan', measurementError=FALSE)
  
  expect_equal( 2.5, mean(Q$beta[2,][!is.na(Q$beta[2,])]) , tol= 0.04 )
  expect_gt( cor( Q$beta0[!is.na(Q$beta0)], 0.2*s[!is.na(Q$beta0)]), 0.95) # this should be change to beta0 at some point.
  expect_gt( cor(Q$beta[1,][!is.na(Q$beta[1,])], as.vector(betaFunc1)[!is.na(Q$beta[1,])]), 0.99)
})

test_that("Works with bandwidth selection returned by ConcurReg", {
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
  
  withError2D <- ConcurReg(vars, outGrid, measurementError = TRUE, diag1D = "none")
  l1 <- sqrt(mean((withError2D$beta0)^2, na.rm = TRUE))
  l2 <- sqrt(mean(((withError2D$beta[1,]-1)^2), na.rm = TRUE))
  l3 <- sqrt(mean(((withError2D$beta[2,]-1)^2), na.rm = TRUE))
  expect_lt( max(c(l1, l2, l3)), 1)
  
  # Sparsify functional data
  set.seed(1)
  X_1sp <- fdapace::Sparsify(X_1, Sup, sparsity)
  set.seed(1)
  Ysp <- fdapace::Sparsify(Y, Sup, sparsity)
  vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
  noError1D <- ConcurReg(vars, outGrid, measurementError = FALSE, diag1D = "all")
  l1 <- sqrt(mean((noError1D$beta0)^2, na.rm = TRUE))
  l2 <- sqrt(mean(((noError1D$beta[1,]-1)^2), na.rm = TRUE))
  l3 <- sqrt(mean(((noError1D$beta[2,]-1)^2), na.rm = TRUE))
  expect_lt( max(c(l1, l2, l3)), 1)
})









