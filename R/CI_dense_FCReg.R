#' @title Bootstrap pointwise confidence intervals for the coefficient functions in functional concurrent regression for densely observed data.
#' @param tGrid A vector of length \eqn{m} with the input time points.
#' @param dat A list of input functional/scalar covariates. 
#' Each field corresponds to a functional (a matrix) or scalar (a vector) variable. 
#' The last entry is assumed to be the functional response if no entry is names \code{'Y'}. 
#' If a field corresponds to a functional variable, it should be an \eqn{n}-by-\eqn{m} matrix,
#' where each row holds the observations for one subject on the common grid \code{tGrid}. 
#' If a field corresponds to a scalar covariate, it should be a vector of length \eqn{n}.
#' @param bw Scalar holding the bandwidth.
#' @param kernel_type Character holding the kernel type (see \code{\link[fdapace]{Lwls1D}} for supported kernels).
#' @param R An integer holding the number of bootstrap replicates. Default: 999.
#' @param level A number taking values in [0,1] determing the confidence level. Default: 0.95.
#'
#' @return A list containing the following fields: 
#' \describe{
#' \item{CI_beta0}{CI for the intercept function --- A data frame holding three variables: 
#' \code{CIgrid} --- the time grid where the CIs are evaluated;
#' \code{CI_beta0.lower} and \code{CI_beta0.upper} --- the lower and upper bounds of the CIs
#' for the intercept function on \code{CIgrid}.}
#' 
#' \item{CI_beta}{ A list containing CIs for the slope functions --- the length of
#' the list is same as the number of covariates. Each list contains the following fields:
#' A data frame holding three variables: \code{CIgrid} --- the time grid where the CIs are evaluated,
#' \code{CI_beta_j.lower} and \code{CI_beta_j.upper} --- the lower and upper bounds of the CIs 
#' for the intercept function on \code{CIgrid} for \eqn{j = 1,2,\dots}.}
#' 
#' \item{CI_R2}{CI the time-varying \eqn{R^2(t)} --- A data frame holding three variables: 
#' \code{CIgrid} --- the time grid where the CIs are evaluated,
#' \code{CI_R2.lower} and \code{CI_R2.upper} --- the lower and upper bounds of the CIs 
#' for the time-varying \eqn{R^2(t)} on \code{CIgrid}.}
#' 
#' \item{level}{The confidence level of the CIs.}
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
#' smres <- smPtFCRegCoef(res, bw =  2.5 / (nGridIn-1), kernel_type = 'epan')
#' res_CI = GetCI_Dense(dat, tGrid, level = 0.95, R = 10, bw = 2.5 / (nGridIn-1), kernel_type = 'epan')
#' beta1 = res_CI$CI_beta[[1]] ##extracting CI for beta1
#' beta1a = beta1$CI_beta1.lower
#' beta1b = beta1$CI_beta1.upper
#' true_beta = beta[1,]  ##extracting true coef beta1 in the simulation setting
#' est_beta = smres$beta[1,] ## ##extracting estimated coef beta1 from 
#' ###fitting the concurrent regression model
#' plot(beta1$CIgrid, beta1a, type= 'l', ylim = c(0,2)) ##plot of lower CI for beta1
#' lines(beta1$CIgrid, beta1b) ##plot of lower CI for beta1
#' lines(beta1$CIgrid, true_beta, col ='red')  
#' ##plot of true coef beta1 in the simulation setting
#' lines(beta1$CIgrid, est_beta, col ='blue') 
#' ##plot of estimated coef beta1 from fitting the concurrent regression model
#' @export

#utils::globalVariables(c("prob", "section", "y"))
GetCI_Dense <- function(dat, tGrid, level = 0.95, R = 10, bw, kernel_type){
  if (length(level) > 1) {
    level = level[1]
    warning("The input level has more than 1 element; only the first one is used.")
  }
  if (level < 0 | level > 1) {
    stop("Invalid input value of level.")
  }
  if (R %% 1 != 0 | R < 0) {
    stop("R should be an positive integer.")
  }
  
  n <- nrow(dat$Y)
  betaMat <- lapply(1:R, function(b) {
    ind <- sample(x = seq_len(n), size = n, replace = TRUE)
    for(j in 1:length(dat)){
      if ( is.matrix(dat[[j]]) ) {
        dat[[j]] = dat[[j]][ind,]
      }else{
        dat[[j]] = dat[[j]][ind]
      }
    }
    res <- ptFCReg(tGrid = tGrid, dat = dat)
    smres <- smPtFCRegCoef(res, bw = bw, kernel_type = 'epan')
    return(list(beta0 = smres$beta0, beta = smres$beta, R2 = smres$R2))
  })
  CI_beta0 = apply(t(sapply(1:R, function(b){
    betaMat[[b]]$beta0
  })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))
  CI_beta0 = data.frame(CI_beta0.lower = CI_beta0[1,], CI_beta0.upper = CI_beta0[2,], CIgrid = tGrid)
  CI_beta = lapply(1:(length(dat)-1), function(j){
    ci_beta_df =  data.frame( t(apply(t(sapply(1:R, function(b){
      betaMat[[b]]$beta[j,]
      })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))))
    names(ci_beta_df) =  c( sprintf("CI_beta%d.lower", j)  ,sprintf("CI_beta%d.upper", j)) 
    ci_beta_df$CIgrid = tGrid
    return(ci_beta_df)
  })
  CI_R2 = apply(t(sapply(1:R, function(b){
    betaMat[[b]]$R2 
  })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))
  CI_R2 = data.frame(CI_R2.lower = CI_R2[1,], CI_R2.upper = CI_R2[2,], CIgrid = tGrid)
  return(list(CI_beta0 = CI_beta0, CI_beta = CI_beta, CI_R2 = CI_R2,
              level = level))
}
  
 