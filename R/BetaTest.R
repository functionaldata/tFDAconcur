rm(list = ls())
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

print(model$workGridY) # workGrid for Y between 6 to 15
plot(beta_1(model$workGridY)) # workGrid for X1 between 3 to 12 
plot(beta_0(model$workGridY))

plot(model$beta[1,])
plot(model$beta0)

