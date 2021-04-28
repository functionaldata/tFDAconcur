#' Smooth the concurrent effects functions in a \code{ptFCReg} object using local linear regression.
#' The local linear regression is implemented using the function \code{\link[fdapace]{Lwls1D}}.
#' @param object An object of class \code{ptFCReg} returned by the function \code{\link{ptFCReg}}.
#' @param bw Scalar holding the bandwidth.
#' @param kernel_type Character holding the kernel type (see \code{\link[fdapace]{Lwls1D}} for supported kernels).
#' @return An object of class \code{ptFCReg}, where the fields \code{beta0} and \code{beta} 
#' hold the smoothed intercept functions and concurrent effects functions, respectively. 
#' See \code{\link{ptFCReg}} for a complete list of the fields. 
#' @examples
#' set.seed(1)
#' n <- 50
#' nGridIn <- 101
#' tGrid <- seq(0, 1, length.out=nGridIn) # Functional data support
#' muX1 <- tGrid * 2 # mean function for X_1
#' sigma <- 1
#' beta0 <- 0
#' beta <- rbind(cos(tGrid), 1.5 + sin(tGrid))
#' Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
#' X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(muX1, n, nGridIn, byrow=TRUE)
#' epsilon <- rnorm(n, sd=sigma)
#' Y <- t(sapply(seq_len(n), function(i) {
#'   beta0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
#' }))
#' dat <- list(X1=X_1, Z1=Z[, 2], Y=Y)
#' res <- ptFCReg(tGrid = tGrid, dat = dat)
#' smres <- smPtFCRegCoef(res, bw = 2.5 / (nGridIn-1), kernel_type = 'epan')
#' @export
smPtFCRegCoef <- function(object, bw, kernel_type) {
  smCoefs <- apply(
    rbind(object$beta0,object$beta), 1, fdapace::Lwls1D,
    bw = bw, kernel_type = kernel_type, win = rep(1,length(object$tGrid)),
    xin = object$tGrid, xout = object$tGrid
  )
  smres <- object
  smres$beta0 <- smCoefs[,1]
  smres$beta <- t(smCoefs[,-1, drop = FALSE])
  return(smres)
}
