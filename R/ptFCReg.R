#' Functional concurrent regression using pointwise multiple linear regression.
#' @param tGrid A vector of length \eqn{m} with the input time points.
#' @param dat A list of input functional/scalar covariates. 
#' Each field corresponds to a functional (a matrix) or scalar (a vector) variable. 
#' The last entry is assumed to be the functional response if no entry is names \code{'Y'}. 
#' If a field corresponds to a functional variable, it should be an \eqn{n}-by-\eqn{m} matrix,
#' where each row holds the observations for one subject on the common grid \code{tGrid}. 
#' If a field corresponds to a scalar covariate, it should be a vector of length \eqn{n}.
#' @return A list containing the following fields: 
#' \describe{
#' \item{beta0}{A vector containing the time-varying intercept evaluated on \code{tGrid}.}
#' \item{beta}{A matrix for the concurrent regression effects, 
#' where rows correspond to different predictors and columns to different time points in \code{tGrid}.}
#' \item{tGrid}{The input \code{tGrid}.}
#' \item{R2}{A vector of the time-varying \eqn{R^2(t)}, evaluated at \eqn{t} in \code{tGrid}.}
#' \item{Ldf}{A list holding the input data, each element of which is a data frame holding 
#' the data observed at one element of \code{tGrid}.}
#' }
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
#' @export
ptFCReg <- function(tGrid, dat) {
  if ( !'Y' %in% names(dat) ) {
    names(dat)[length(dat)] <- 'Y'
  }
  for(i in 1:(length(dat)-1)){
    if(is.na(names(dat)[i])){
      names(dat)[i] = paste0('X',sprintf("%d",i))
    }
  }
  if ( !is.matrix(dat$Y) ) {
    stop("The field in dat for the response is not a matrix.")
  }
  nSubj <- nrow(dat$Y)
  lenGrid <- length(tGrid)
  if (!all(sapply(dat, function(var) { is.vector(var) | is.matrix(var)}))) {
    stop("Fields of dat are not vectors or matrices.")
  }
  
  Z <- sapply(dat, is.vector)
  # indices for scalar covariates in dat:
  indZ <- which(Z)
  # indices for functional covariates in dat:
  indX <- which(sapply(dat[names(dat) != 'Y'], is.matrix)) 
  
  if ( length(indX) > 0 ) {
    if ( any( abs( sapply(dat[indX], ncol) - lenGrid) > 0 ) ) {
      stop("Columns in the matrices corresponding to the functional variables do not match with tGrid.")
    }
    if ( any( abs( sapply(dat[indX], nrow) - nSubj) > 0 ) ) {
      stop("Numbers of rows in the matrices corresponding to the functional variables are not all the same.")
    }
  }
  if ( length(indZ) > 0 ) {
    if ( any( abs( sapply(dat[indZ], length) - nSubj) > 0 ) ) {
      stop("Lengths of the vectors corresponding to scalar covariates are not all the same as the number of subjects in the functional response.")
    }
  }
  
  # A list of data frames, each corresponding to the data observed at one time point
  Ldf <- lapply(seq_len(lenGrid), function(j) {
    df <- as.data.frame(lapply(dat, function(var) {
      if (is.vector(var)) {
        return(var)
      } else return(var[,j])
    }))
    names(df) <- names(dat)
    df
  })
  
  Lmod <- lapply(Ldf, function(df) {
    stats::lm( Y ~ ., data = df )
  })
  coef <- sapply(Lmod, coef)
  R2 <- sapply(Lmod, function(mod) summary(mod)$r.sq)
  res <- list(beta0 = coef[1,], beta = coef[-1, , drop = FALSE], tGrid = tGrid, R2 = R2, Ldf = Ldf)
  rownames(res$beta) <- names(dat[names(dat) != 'Y'])
  class(res) <- "ptFCReg"
  return(res)
}


# ptFCReg <- function(tGrid, dat) {
#   if ( !'Y' %in% names(dat) ) {
#     names(dat)[length(dat)] <- 'Y'
#   }
#   if ( !is.matrix(dat$Y) ) {
#     stop("The field in dat for the response is not a matrix.")
#   }
#   nSubj <- nrow(dat$Y)
#   lenGrid <- length(tGrid)
#   if (!all(sapply(dat, function(var) { is.vector(var) | is.matrix(var)}))) {
#     stop("Fields of dat are not vectors or matrices.")
#   }
#   
#   Z <- sapply(dat, is.vector)
#   # indices for scalar covariates in dat:
#   indZ <- which(Z)
#   # indices for functional covariates in dat:
#   indX <- which(sapply(dat[names(dat) != 'Y'], is.matrix)) 
#   
#   if ( length(indX) > 0 ) {
#     if ( any( abs( sapply(dat[indX], ncol) - lenGrid) > 0 ) ) {
#       stop("Columns in the matrices corresponding to the functional variables do not match with tGrid.")
#     }
#     if ( any( abs( sapply(dat[indX], nrow) - nSubj) > 0 ) ) {
#       stop("Numbers of rows in the matrices corresponding to the functional variables are not all the same.")
#     }
#   }
#   if ( length(indZ) > 0 ) {
#     if ( any( abs( sapply(dat[indZ], length) - nSubj) > 0 ) ) {
#       stop("Lengths of the vectors corresponding to scalar covariates are not all the same as the number of subjects in the functional response.")
#     }
#   }
#   
#   # A list of data frames, each corresponding to the data observed at one time point
#   Ldf <- lapply(seq_len(lenGrid), function(j) {
#     df <- as.data.frame(lapply(dat, function(var) {
#       if (is.vector(var)) {
#         return(var)
#       } else return(var[,j])
#     }))
#     names(df) <- names(dat)
#     df
#   })
#   
#   Lmod <- lapply(Ldf, function(df) {
#     stats::lm( Y ~ ., data = df )
#   })
#   coef <- sapply(Lmod, coef)
#   R2 <- sapply(Lmod, function(mod) summary(mod)$r.sq)
#   res <- list(beta0 = coef[1,], beta = coef[-1, , drop = FALSE], tGrid = tGrid, R2 = R2, Ldf = Ldf)
#   rownames(res$beta) <- names(dat[names(dat) != 'Y'])
#   class(res) <- "ptFCReg"
#   return(res)
# }
