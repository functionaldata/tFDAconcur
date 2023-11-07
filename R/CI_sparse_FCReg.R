#' @title Bootstrap pointwise confidence intervals for the coefficient functions in functional concurrent regression for sparsely observed data.
#' 
#' @param vars A list of input functional/scalar covariates.
#' Each field corresponds to a functional (a list) or scalar (a vector) covariate. 
#' The last entry is assumed to be the response if no entry is named 'Y'.
#' If a field corresponds to a functional covariate, it should have two fields: 'Lt', a list of time points, and 'Ly', a list of functional values.
#' @param outGrid A vector with the output time points, which need to be within 5\% and 95\% of the range of functional covariates. If NULL, outGrid will be generated from 5\% to 95\% of the range of functional covariates with 51 grids for fitting. Default: NULL
#' @param level A number taking values in [0,1] determining the confidence level. Default: 0.95.
#' @param R An integer holding the number of bootstrap replicates. Default: 999.
#' @param userBwMu A scalar/vector bandwidth used for smoothing the mean function. Each entry in the vector represents the bandwidth used for the corresponding covariate in vars. For the scalar covariates, you can input 0 as a placeholder. If you only input a scalar, the function will use the same bandwidth to smooth all mean functions. --- a scalar/vector of positive numeric -
#' default: NULL --- if no scalar/vector value is provided, the bandwidth value for the smoothed mean function is chosen using 'GCV'; 
#' @param userBwCov A scalar/vector bandwidth used for smoothing the auto or cross-covariances. If you use 1D smoothing for the diagonal line of the covariance (diag1D="all"), only one scalar input is needed. If you use 2D smoothing for the covariance (diag1D="none"), a vector of bandwidth is required. Each entry in the vector represents the bandwidth used for the corresponding covariate in vars. For the scalar covariates, you can input 0 as a placeholder. --- a scalar/vector of positive numeric - 
#' default: NULL --- if no scalar/vector is provided, the bandwidth value for the smoothed cross-covariance function is chosen using 'GCV';
#' @param kern Smoothing kernel choice, common for mu and covariance; 
#' "rect", "gauss", "epan", "gausvar", "quar" (default: "gauss")
#' @param measurementError Assume measurement error in the data; logical - default: FALSE. If TRUE the diagonal raw covariance will be removed when smoothing.
#' @param diag1D  A string specifying whether to use 1D smoothing for the diagonal line of the covariance. 
#' 'none': don't use 1D smoothing; 'all': use 1D for both auto- and cross-covariances. (default : 'all')
#' @param useGAM Use GAM smoothing instead of local linear smoothing (semi-parametric option);  logical - default: FALSE.
#' 
#' @details If measurement error is assumed, the diagonal elements of the raw covariance will be removed. This could result in highly unstable estimate 
#' if the design is very sparse, or strong seasonality presents. 
#' WARNING! For very sparse functional data, setting measurementError = TRUE is not recommended.
#' 
#' @return A list containing the following fields: 
#' \describe{
#' \item{CI_beta0}{CI for the intercept function --- A data frame holding three variables: 
#' \code{CIgrid} --- the time grid where the CIs are evaluated,
#' \code{CI_beta0.lower} and \code{CI_beta0.upper} --- the lower and upper bounds of the CIs 
#' for the intercept function on \code{CIgrid}.}
#' 
#' \item{CI_beta}{ A list containing CIs for the slope functions --- the length of
#' the list is same as the number of covariates. Each list contains the following fields:
#' A data frame holding three variables: \code{CIgrid} --- the time grid where the CIs are evaluated,
#' \code{CI_beta_j.lower} and \code{CI_beta_j.upper} --- the lower and upper bounds of the CIs 
#' for the coefficient functions on \code{CIgrid} for \eqn{j = 1,2,\dots}.}
#' 
#' \item{CI_R2}{CI the time-varying \eqn{R^2(t)} --- A data frame holding three variables: 
#' \code{CIgrid} --- the time grid where the CIs are evaluated,
#' \code{CI_R2.lower} and \code{CI_R2.upper} --- the lower and upper bounds of the CIs 
#' for the time-varying \eqn{R^2(t)} on \code{CIgrid}.}
#' 
#' \item{level}{The confidence level of the CIs.}
#' }
#' @examples 
#'
#' set.seed(1)
#' n <- 30
#' nGridIn <- 100
#' sparsity <- 5:10 # Sparse data sparsity
#' T <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
#' bw <- 0.1
#' outGrid <- round(seq(min(T), 1, by=0.05), 2)
#' # Simulate functional data
#' mu <- T * 2 # mean function for X_1
#' sigma <- 1
#' 
#' beta_0 <- 0
#' beta <- rbind(cos(T), 1.5 + sin(T))
#' beta_2 <- 1
#' 
#' Z <- MASS::mvrnorm(n, rep(0, 2), diag(2))
#' X_1 <- Z[, 1, drop=FALSE] %*% matrix(1, 1, nGridIn) + matrix(mu, n, nGridIn, byrow=TRUE)
#' epsilon <- rnorm(n, sd=sigma)
#' Y <- matrix(NA, n, nGridIn)
#' for (i in seq_len(n)) {
#'   Y[i, ] <- beta_0 + beta[1,] * X_1[i, ] + beta[2,] * Z[i, 2] + epsilon[i]
#' }
#' 
#' # Sparsify functional data
#' set.seed(1)
#' X_1sp <- fdapace::Sparsify(X_1, T, sparsity)
#' Ysp <- fdapace::Sparsify(Y, T, sparsity)
#' vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
#' res <-  GetCI_Sparse(vars, outGrid[-c(1,21)], level = 0.95, R = 2,
#'                      userBwMu = c(.1,.1,.1), userBwCov = c(.1,.1,.1),
#'                      kern='gauss', measurementError=TRUE, diag1D='none',
#'                      useGAM = FALSE)
#' @export

GetCI_Sparse = function(vars, outGrid = NULL, level = 0.95, R = 999, userBwMu = NULL, userBwCov = NULL,  
                        kern = "gauss", measurementError = FALSE, diag1D = "all", useGAM = FALSE){
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
  
  n <- lengthVars(vars)
  p <- length(vars) - 1
  
  if (is.null(names(vars))){
    names(vars) <- c(paste0('X', seq_len(length(vars) - 1)), 'Y')
  }
  
  if ('Y' %in% names(vars)) {
    vars <- c(vars[names(vars) != 'Y'], vars['Y'])
    
  } else if (names(vars)[length(vars)] == '') {
    names(vars)[length(vars)] <- 'Y'
  }
  
  Yname <- names(vars)[length(vars)]
  
  # Handle NaN, int to double
  vars[sapply(vars, is.list)] <- lapply(
    vars[sapply(vars, is.list)], 
    function(v) HandleNumericsAndNAN(v[['Ly']], v[['Lt']])
  )
  # outGrid <- as.numeric(outGrid)
  
  # handle ourGrid range
  temp <- lapply(
    vars[sapply(vars, is.list)], 
    function(v) {
      return(range(unlist(v[['Lt']])))
    }
  )

  l <- max(unlist(lapply(temp, function(v){return(v[1])})))
  u <- min(unlist(lapply(temp, function(v){return(v[2])})))
  
  l5 = l+0.05*(u-l)
  u95 = u-0.05*(u-l)
  if(is.null(outGrid)){
    outGrid = seq(l5, u95, length.out=51)
  }else{
    if(min(outGrid)<l5 | max(outGrid)>u95){
      stop('The input "outGrid" is out of the 5% and 95% of the range of Lt, the list of time points in the "vars".')
    }
  }
  
  betaMat <- lapply(1:R, function(b) {
    OutOfRange = TRUE
    while (OutOfRange) {
      ind = sample(x = seq_len(n), size = n, replace = TRUE)
      temp_ind = lapply(vars[sapply(vars, is.list)], function(v) {return(range(unlist(v[['Lt']][ind])))})
      l_ind = max(unlist(lapply(temp_ind, function(v){return(v[1])})))
      u_ind = min(unlist(lapply(temp_ind, function(v){return(v[2])})))
      if(l_ind <= l5 & u_ind >= u95){
        OutOfRange = FALSE
      }
    }
    vars_ind = vars
    for(j in 1:(p+1)){
      if ( is.list(vars[[j]]) ) {
        vars_ind[[j]]$Lt = vars[[j]]$Lt[ind]
        vars_ind[[j]]$Ly = vars[[j]]$Ly[ind]
      }else{
        vars_ind[[j]] = vars[[j]][ind]
      }
    }
    #res = ConcurReg(vars_ind, outGrid, userBwMu = .5, userBwCov=.5,  kern='gauss', measurementError=TRUE, diag1D='none', useGAM = FALSE, returnCov=TRUE)
    res = ConcurReg(vars_ind, outGrid, userBwMu = userBwMu, userBwCov = userBwCov,  kern = kern,
           measurementError = measurementError, diag1D = diag1D, useGAM = useGAM)
    length(res$beta0)
    return(list(beta0 = res$beta0, beta = res$beta, R2 = res$R2, outGrid = res$outGrid))
  })

  CI_beta0 = apply(t(sapply(1:R, function(b){
    betaMat[[b]]$beta0[ stats::complete.cases(betaMat[[b]]$beta0)]
  }, simplify = TRUE)), 2,
    stats::quantile, c((1-level)/2, 1-(1-level)/2))
  
  
  CI_beta0 = data.frame(CI_beta0.lower = CI_beta0[1,], CI_beta0.upper = CI_beta0[2,], CIgrid = betaMat[[1]]$outGrid)
  CI_beta = lapply(1:p, function(j){
    ci_beta_df =  data.frame( t(apply(t(sapply(1:R, function(b){
      betaMat[[b]]$beta[j,][ stats::complete.cases(betaMat[[b]]$beta[j,])]
    })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))))
    names(ci_beta_df) =  c( sprintf("CI_beta%d.lower", j)  ,sprintf("CI_beta%d.upper", j)) 
    ci_beta_df$CIgrid = betaMat[[1]]$outGrid
    return(ci_beta_df)
  })
  names(CI_beta) = sapply(1:p, function(j) { sprintf("CI_beta%d", j)})
  
  CI_R2 = apply(t(sapply(1:R, function(b){
    betaMat[[b]]$R2[stats::complete.cases(betaMat[[b]]$R2)]
  })), 2, stats::quantile, c((1-level)/2, 1-(1-level)/2))
  CI_R2 = data.frame(CI_R2.lower = CI_R2[1,], CI_R2.upper = CI_R2[2,], CIgrid = betaMat[[1]]$outGrid)
  return(list(CI_beta0 = CI_beta0, CI_beta = CI_beta, CI_R2 = CI_R2,
              outGrid = outGrid, level = level))
  
}


