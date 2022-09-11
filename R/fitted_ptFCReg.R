#' Fitted functional responses from a \code{ptFCReg} object.
#' @param object An object of class \code{ptFCReg}, returned by \code{\link{ptFCReg}} or \code{\link{smPtFCRegCoef}}.
#' @return An \eqn{n}-by-\eqn{m} matrix, where each row corresponds to one subject,
#' and each column corresponds to a time point in \code{object$tGrid}.
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
#' fit_res <- stats::fitted(res)
#' fit_smres <- stats::fitted(smres)
#' @export
fitted_ptFCReg <- function(object) {
  fit <- sapply(seq_along(object$tGrid), function(j) {
    df <- object$Ldf[[j]]
    as.matrix(df[names(df)!='Y']) %*% object$beta[,j]
  })
  fit <- fit + matrix(object$beta0, nrow = nrow(fit), ncol = ncol(fit), byrow = TRUE)
  return(fit)
}
