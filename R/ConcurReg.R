#' Functional Concurrent Regression using 2D smoothing 
#' 
#' Functional concurrent regression with dense or sparse functional data for scalar or functional dependent variables. Note: function-to-scalar regression can also be handled using the VCAM function in fdapace. 
#' 
#' @param vars A list of input functional/scalar covariates.
#'  Each field corresponds to a functional (a list) or scalar (a vector) covariate. 
#'  The last entry is assumed to be the response if no entry is named 'Y'.
#'  If a field corresponds to a functional covariate, it should have two fields: 'Lt', a list of time points, and 'Ly', a list of functional values.
#' @param outGrid A vector of output time points.
#' @param userBwMu A scalar/vector bandwidth used for smoothing the mean function. Each entry in the vector represents the bandwidth used for the corresponding covariate in vars. For the scalar covariates, you can input 0 as a placeholder. If you only input a scalar, the function will use the same bandwidth to smooth all mean functions. --- a scalar/vector of positive numeric -
#' default: NULL --- if no scalar/vector value is provided, the bandwidth value for the smoothed mean function is chosen using 'GCV'; 
#' @param userBwCov A scalar/vector bandwidth used for smoothing the auto or cross-covariances. If you use 1D smoothing for the diagonal line of the covariance (diag1D="all"), only one scalar input is needed. If you use 2D smoothing for the covariance (diag1D="none"), a vector of bandwidth is required. Each entry in the vector represents the bandwidth used for the corresponding covariate in vars. For the scalar covariates, you can input 0 as a placeholder. --- a scalar/vector of positive numeric - 
#' default: NULL --- if no scalar/vector is provided, the bandwidth value for the smoothed cross-covariance function is chosen using 'GCV';
#' @param kern Smoothing kernel choice, common for mu and covariance; "rect", "gauss", "epan", "gausvar", "quar" - default: "gauss".
#' @param measurementError Assume measurement error in the data; logical - default: FALSE. If TRUE the diagonal raw covariance will be removed when smoothing.
#' @param diag1D  A string specifying whether to use 1D smoothing for the diagonal line of the covariance. 
#' 'none': don't use 1D smoothing; 'all': use 1D for both auto- and cross-covariances. (default : 'all')
#' @param useGAM Use GAM smoothing instead of local linear smoothing (semi-parametric option);  logical - default: FALSE.
#' @param returnCov Return the covariance surfaces, which is a four dimensional array. The first two dimensions correspond to outGrid
#'  and the last two correspond to the covariates and the response, i.e. (i, j, k, l) entry being Cov(X_k(t_i), X_l(t_j));  logical - default: FALSE.
#' 
#' @details If measurement error is assumed, the diagonal elements of the raw covariance will be removed. This could result in highly unstable estimate if the design is very sparse, or strong seasonality presents. 
#' WARNING! For very sparse functional data, setting measurementError = TRUE is not recommended.
#' @return A list containing the following fields:
#' \item{beta}{A matrix for the concurrent regression effects, where rows correspond to different predictors and columns to different time points.}
#' \item{beta0}{A vector containing the time-varying intercept.}
#' \item{outGrid}{A vector of the output time points.}
#' \item{cov}{A 4-dimensional array for the (cross-)covariance surfaces, with the (i, j, k, l) entry being Cov(X_k(t_i), X_l(t_j))}
#' \item{R2}{A vector of the time-varying R2.}
#' \item{n}{The sample size.}
#' @examples 
#' # Y(t) = beta_0(t) + beta_1(t) X_1(t) + beta_2(t) Z_2 + epsilon
#' #
#' # Settings
#' set.seed(1)
#' n <- 30
#' nGridIn <- 50
#' sparsity <- 5:10 # Sparse data sparsity
#' T <- round(seq(0, 1, length.out=nGridIn), 4) # Functional data support
#' bw <- 0.1
#' outGrid <- round(seq(min(T), 1, by=0.05), 2)
#'  outGrid <- seq(min(T), max(T), by=0.05)
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
#'   Y[i, ] <- beta_0 + beta[1,]*X_1[i, ] + beta[2,]*Z[i, 2] + epsilon[i]
#' }
#' 
#' # Sparsify functional data
#' set.seed(1)
#' X_1sp <- fdapace::Sparsify(X_1, T, sparsity)
#' set.seed(1)
#' Ysp <- fdapace::Sparsify(Y, T, sparsity)
#' vars <- list(X_1=X_1sp, Z_2=Z[, 2], Y=Ysp)
#' res2 <- ConcurReg(vars, outGrid, userBwMu=c(.5,0,.5), userBwCov=c(.5,0,.5), kern='gauss',
#'  measurementError=TRUE, diag1D='none', useGAM = FALSE, returnCov=TRUE)
#' @references
#' \itemize{
#' \item \cite{Yao, F., Müller, H.G., Wang, J.L. "Functional Linear Regression Analysis for Longitudinal Data." Annals of Statistics 33, (2005): 2873-2903.(Dense data)} 
#' \item \cite{Sentürk, D., Müller, H.G. "Functional varying coefficient models for longitudinal data." J. American Statistical Association, 10, (2010): 1256--1264.}
#' \item \cite{Sentürk, D., Nguyen, D.V. "Varying Coefficient Models for Sparse Noise-contaminated Longitudinal Data", Statistica Sinica 21(4), (2011): 1831-1856. (Sparse data)} 
#' }
#' @export


ConcurReg <- function(vars, outGrid, userBwMu=NULL, userBwCov=NULL,
                      kern='gauss', measurementError=FALSE, diag1D='all', useGAM = FALSE, returnCov=TRUE) {
  
  n <- lengthVars(vars)
  p <- length(vars) - 1
  if (p == 0){
    stop('Too few covariates.')
  }
  if(!is.null(userBwCov)){
    if(diag1D == "all" & length(userBwCov) != 1){
      stop('For 1D smoothing of the covariance diagonal line, userBw only takes one scalar value.')
    }
    
  }
  
  if(!is.null(userBwCov) & length(vars) != length(userBwCov)){
    if(diag1D == "all" & length(userBwCov) == 1){
      userBwCov = rep(userBwCov, length(vars))
    }else{
      if(length(userBwCov) == 1){
        userBwCov = rep(userBwCov, length(vars))
        warning('The length of the input userBwCov is not consistent with the length of the input vars. Use the same userBwCov for all input vars.') 
      }else{
        stop('The length of the input userBwCov is not consistent with the length of the input vars.') 
      }
    }
    
  }
  if(!is.null(userBwMu)){
    if(length(userBwMu)==1){
      userBwMu = rep(userBwMu, length(vars))
    }
  }
  
  if (is.null(names(vars))){
    names(vars) <- c(paste0('X', seq_len(length(vars) - 1)), 'Y')
  }
    
  if ('Y' %in% names(vars)) {
    if(!is.null(userBwCov)){
      userBwCov <- c(userBwCov[names(vars) != 'Y'], userBwCov[names(vars) == 'Y'])
    }
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
  grid.index <- which((outGrid >= l) & (outGrid <= u))
  grid.full <- outGrid
  outGrid <- outGrid[grid.index]
  
  # raise a warning, outgrid is out of the input Lt
  if(length(grid.full) > length(outGrid)){
    warning('The input "outGrid" is out of the range of Lt, the list of time points in the "vars", in which case output "beta0" and "beta" will contain NA on the boundary.')
  }
  # De-mean.
  demeanedRes <- demean(vars, userBwMu, kern)
  vars <- demeanedRes[['xList']]
  muList <- demeanedRes[['muList']]
  
  allCov <- MvCov(vars, userBwCov, outGrid, kern, 
                  measurementError, center=FALSE, 
                  diag1D)
  beta <- sapply(seq_len(dim(allCov)[1]), function(i) {
    tmpCov <- allCov[i, i, , ]
    beta_ti <- qr.solve(tmpCov[1:p, 1:p], tmpCov[1:p, p + 1])
    beta_ti
  })
  if (is.null(nrow(beta)))
    beta <- matrix(beta, 1)
  rownames(beta) <- names(vars)[-length(vars)]
  
  # coefficient of determination: 
  #   R2 = stats::cov(X, Y)' var(X)^{-1} stats::cov(X, Y) / var(Y)
  R2 <- sapply(seq_len(dim(allCov)[1]), function(i) {
    tmpCov <- allCov[i, i, , ]
    tmpCov[p + 1, 1:p, drop=FALSE] %*% beta[, i, drop=FALSE] / tmpCov[p + 1, p + 1]
  })
  
  R2<-sapply(R2,function(x){max(min(x,1),0)}) #capping the R2 value between 0 and 1.
  
  muBeta <- sapply(seq_len(p), function(j) {
    if (!is.function(muList[[j]])) { # scalar mean
      beta[j, ] * rep(muList[[j]], length(outGrid))
    } else { # functional mean
      beta[j, ] * muList[[j]](outGrid)
    }
  })
  beta0 <- muList[[Yname]](outGrid) - colSums(t(muBeta))
  
  ## enlarge output
  ## 
  beta0.full <- rep(NA,length(grid.full))
  beta0.full[grid.index] <- beta0
  if(is.vector(beta)){
    beta.full <- rep(NA, length(grid.full))
    beta.full[grid.index] <- beta
  }else{
    beta.full <- matrix(NA, dim(beta)[1], length(grid.full))
    beta.full[, grid.index] <- beta
  }
  allCov.full <- array(NA, c(length(grid.full), length(grid.full), dim(allCov)[3], dim(allCov)[4]))
  allCov.full[grid.index, grid.index,,] <- allCov
  
  #res <- list(beta=beta, beta0 = beta0, outGrid=outGrid, cov=allCov, R2=R2, n=n)
  res <- list(beta=beta.full, beta0 = beta0.full, outGrid=grid.full, cov=allCov.full, R2=R2, n=n)
  if (!returnCov)
    res[['cov']] <- NULL
  res
}

demean <- function(vars, userBwMu, kern) {
  tmp <- lapply(1:length(vars), function(i) {
    x = vars[[i]]
    if(is.null(userBwMu)){
      userBwMu_i = NULL
    }else{
      userBwMu_i = userBwMu[i]
    }
    
    if (is.numeric(x)) { # scalar
      xmu <- mean(x)
      x <- x - xmu
    } else if (is.list(x)) { # functional
      Tin <- sort(unique(unlist(x[['Lt']])))
      
      if(is.null(userBwMu_i)){ # bandwidth choice for mean function is using GCV
        optns <- fdapace::SetOptions(x[['Ly']], x[['Lt']], list(userBwMu=userBwMu_i, methodBwMu ='GCV', kernel=kern))
        xmu <- GetSmoothedMeanCurve(x[['Ly']], x[['Lt']] , Tin, Tin[1], optns)[['mu']]
      } else{
        xmu <- GetSmoothedMeanCurve(x[['Ly']], x[['Lt']] , Tin, Tin[1],
                                              list(userBwMu=userBwMu_i, kernel=kern))[['mu']]
      }
      
      #muFun <- stats::approxfun(Tin, xmu)
      muFun <- stats::approxfun(Tin, xmu, rule=2)
      x[['Ly']] <- lapply(1:length(x[['Ly']]), function(i)
        x[['Ly']][[i]]- muFun(x[['Lt']][[i]]))
      xmu <- muFun
    }
    
    list(x=x, mu=xmu)
  })
  
  xList <- lapply(tmp, `[[`, 'x')
  muList <- lapply(tmp, `[[`, 'mu')
  
  names(xList) = names(vars)
  names(muList) = names(vars)
  
  list(xList = xList, muList = muList)
}

## Multivariate function/scalar covariance.
# INPUTS: same as FCReg
# Output: a 4-D array containing the covariances. The first two dimensions corresponds to 
# time s and t, and the last two dimensions correspond to the variables taken covariance upon.

MvCov <- function(vars, userBwCov, outGrid,
                  kern, measurementError=TRUE, 
                  center, diag1D='none') {
  if (!is.list(vars) || length(vars) < 1)
    stop('`vars` needs to be a list of length >= 1')
  
  if (diag1D == 'all' && measurementError) {
    stop("Cannot assume measurement error when diag1D == 'all'")
  }
  isFuncVars <- sapply(vars, is.list)
  p <- length(isFuncVars)
  pFunc <- sum(isFuncVars)
  pScaler <- sum(!isFuncVars)
  
  if (any(isFuncVars)) {
    tAll <- do.call(c, lapply(vars[isFuncVars], function(x) unlist(x[['Lt']])))
    Tin <- sort(unique(tAll))
    
    if (missing(outGrid))
      outGrid <- Tin
    lenoutGrid <- length(outGrid)
    
  } else {
    stop('No functional observation found')
  }
  
  # First two dimensions are for s, t, and the last two dimensions are for matrix of # random variables.
  res <- array(NA, c(lenoutGrid, lenoutGrid, p, p))
  for (j in seq_len(p)) {
    for (i in seq_len(p)) {
      #print(c(i,j))
      if (j <= i) {
        use1D <- diag1D == 'all'
        if(is.null(userBwCov)){
          covRes <- uniCov(X = vars[[i]], Y = vars[[j]], 
                           userBwCov,
                           outGrid,
                           kern, 
                           rmDiag = (i == j) && measurementError,
                           center, use1D)
        }else{
          covRes <- uniCov(X = vars[[i]], Y = vars[[j]], 
                           userBwCov[c(i,j)],
                           outGrid,
                           kern, 
                           rmDiag = (i == j) && measurementError,
                           center, use1D)
        }
        
        if (attr(covRes, 'covType') %in% c('FF', 'SS')){
          res[, , i, j] <- covRes
          } else {
          if (nrow(covRes) == 1){   # stats::cov(scalar, function){
            res[, , i, j] <- matrix(covRes, lenoutGrid, lenoutGrid, byrow=TRUE)
          } else {                    # stats::cov(function, scalar)
            res[, , i, j] <- matrix(covRes, lenoutGrid, lenoutGrid, byrow=FALSE)
          }
        }
      } else { # fill up the symmetric stats::cov(y, x)
        res[, , i, j] <- t(res[, , j, i])
      }
    }
  }
  
  return(res)
}

## Univariate function/scalar covariance.
# rmDiag: whether to remove the diagonal of the raw covariance. Ignored if 1D smoother is used.
# center: whether to center the covariates before calculate covariance.
# use1D: whether to use 1D smoothing for estimating the diagonal covariance.

#X <- vars[[3]]
#Y <- vars[[2]]
uniCov <- function(X, Y, userBwCov, outGrid, kern='gauss', 
                   rmDiag=FALSE, center=TRUE, use1D=FALSE) {
  flagScalerFunc <- FALSE
  # Force X to be a function in the scalar-function case.
  if (!is.list(X) && is.list(Y)) {
    flagScalerFunc <- TRUE
    tmp <- X
    X <- Y
    Y <- tmp
    
    if(!is.null(userBwCov)){
      userBwCov = rev(userBwCov)
    }
  }
  
  # Scalar-scalar
  if (!is.list(X) && !is.list(Y)) {
    res <- stats::cov(X, Y)
    attr(res, 'covType') <- 'SS'
    
    # Scalar-function    
  } else if (is.list(X) && !is.list(Y)) {
    Tin <- sort(unique(unlist(X[['Lt']])))
    if (center) {
      
      if(is.null(userBwCov)){ # bandwidth choice for mean function is using GCV
        optns <- fdapace::SetOptions(X[['Ly']], X[['Lt']], list(userBwMu=userBwCov, methodBwMu ='GCV', kernel=kern))
        Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], Tin, Tin[1], optns)[['mu']]
      } else{
        Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], Tin, Tin[1],
                                    list(userBwMu=userBwCov[1], kernel=kern))[['mu']]
      }
      Ymu <- mean(Y)
    } else {
      Xmu <- rep(0, length(Tin))
      Ymu <- 0
    }
    res <- fdapace::GetCrCovYZ(userBwCov[1], Y, Ymu, X[['Ly']], X[['Lt']], Xmu, Tin, kern)[['smoothedCC']]
    res <- as.matrix(fdapace::ConvertSupport(Tin, outGrid, mu=res))
    #res <- as.matrix(res)
    if (flagScalerFunc) 
      res <- t(res)
    
    attr(res, 'covType') <- 'FS'
    
    # function-function  
  } else {
    TinX <- sort(unique(unlist(X[['Lt']])))
    TinY <- sort(unique(unlist(Y[['Lt']])))
    noutGrid <- length(outGrid)
    if (center) {
      if (min(TinX) > min(outGrid) || min(TinY) > min(outGrid) || 
          max(TinY) < max(outGrid) || max(TinX) < max(outGrid))
        stop('Observation time points coverage too low')
      
      if(is.null(userBwCov)){ # bandwidth choice for mean function is using GCV
        optns <- fdapace::SetOptions(X[['Ly']], X[['Lt']], list(userBwMu=userBwCov, methodBwMu ='GCV', kernel=kern))
        Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], TinX, TinX[1], optns)[['mu']]
      } else{
        Xmu <- GetSmoothedMeanCurve(X[['Ly']], X[['Lt']], TinX, TinX[1],
                                              list(userBwMu=userBwCov[1], kernel=kern))[['mu']]
      }
      
      if(is.null(userBwCov)){
        optns <- fdapace::SetOptions(Y[['Ly']], Y[['Lt']], list(userBwMu=userBwCov, methodBwMu ='GCV', kernel=kern))
        Ymu <- GetSmoothedMeanCurve(Y[['Ly']], Y[['Lt']], TinX, TinX[1], optns)[['mu']]
      }else{
        Ymu <- GetSmoothedMeanCurve(Y[['Ly']], Y[['Lt']], TinY, TinY[1], list(userBwMu=userBwCov[2], kernel=kern))[['mu']]
      }
    } else {
      Xmu <- rep(0, length(TinX))
      Ymu <- rep(0, length(TinY))
    }
    names(Xmu) <- TinX
    names(Ymu) <- TinY
    
    if (use1D) { 
      Xvec <- unlist(X[['Ly']])
      Yvec <- unlist(Y[['Ly']])
      tvecX <- unlist(X[['Lt']])
      tvecY <- unlist(Y[['Lt']])
      if (!identical(tvecX, tvecY)){
        stop('Cannot use 1D covariance smoothing if the observation time points for X and Y are different')
      }
      
      ord <- order(tvecX)
      tvecX <- tvecX[ord]
      Xvec <- Xvec[ord]
      Yvec <- Yvec[ord]
      Xcent <- Xvec - Xmu[as.character(tvecX)]
      Ycent <- Yvec - Ymu[as.character(tvecX)]
      
      if(is.null(userBwCov)){
        optns <- fdapace::SetOptions(X[['Ly']], X[['Lt']], list(userBwMu=userBwCov, methodBwMu ='GCV', kernel=kern))
        bw_mu =  unlist(GCVLwls1D1(yy = Xcent * Ycent, tt = tvecX, kernel = kern, npoly = 1, nder = 0, dataType = optns$dataType) )[1] 
        covXY <- fdapace::Lwls1D(bw=bw_mu, kern, npoly=1L, nder=0L, xin=tvecX, yin=Xcent * Ycent, win=rep(1, length(tvecX)), xout=outGrid)
      }else{
        covXY <- fdapace::Lwls1D(bw=userBwCov[1], kern, npoly=1L, nder=0L, xin=tvecX, yin=Xcent * Ycent, win=rep(1, length(tvecX)), xout=outGrid)
      }
      
      res <- matrix(NA, noutGrid, noutGrid)
      diag(res) <- covXY
    } else { # use 2D smoothing
      
      if(is.null(userBwCov)){
        tmp <- fdapace::GetCrCovYX(bw1 = userBwCov, bw2 = userBwCov, 
                                   Ly1 = X[['Ly']], Lt1 = X[['Lt']], Ymu1 = Xmu, 
                                   Ly2 = Y[['Ly']], Lt2 = Y[['Lt']], Ymu2 = Ymu, 
                                   rmDiag=rmDiag, kern = kern, bwRoutine = 'grid-search')
      }else{
        tmp <- fdapace::GetCrCovYX(bw1 = userBwCov[1], bw2 = userBwCov[2], 
                                   Ly1 = X[['Ly']], Lt1 = X[['Lt']], Ymu1 = Xmu, 
                                   Ly2 = Y[['Ly']], Lt2 = Y[['Lt']], Ymu2 = Ymu, 
                                   rmDiag=rmDiag, kern = kern, bwRoutine = 'grid-search')
      }
      
      # if(snippet){
      #   tmp <- fdapace::GetCrCovYX(userBwCov, userBwCov, X[['Ly']], X[['Lt']], Xmu,
      #                   Y[['Ly']], Y[['Lt']], Ymu, rmDiag=rmDiag, kern=kern, bwRoutine = 'grid-search')
      # }else{
      #   tmp <- fdapace::GetCrCovYX(userBwCov, userBwCov, X[['Ly']], X[['Lt']], Xmu,
      #                     Y[['Ly']], Y[['Lt']], Ymu, rmDiag=rmDiag, kern=kern)
      # }
      
      gd <- tmp[['smoothGrid']]
      res <- matrix(
        fdapace:::interp2lin(as.numeric(gd[, 1]),
                             as.numeric(gd[, 2]),
                             matrix(as.numeric(tmp[['smoothedCC']]),
                                    nrow(tmp[['smoothedCC']]),
                                    ncol(tmp[['smoothedCC']])),
                             rep(as.numeric(outGrid), times=noutGrid),
                             rep(as.numeric(outGrid), each=noutGrid)),
        noutGrid, noutGrid)
      
      # rawCC <- GetRawCrCovFuncFunc(Ly1 = X[['Ly']], Lt1 = X[['Lt']], Ymu1 = Xmu, Ly2 = Y[['Ly']], Lt2 = Y[['Lt']], Ymu2 = Ymu)
      # if (rmDiag) {
      #   diagInd <- rawCC$tpairn[, 1] == rawCC$tpairn[, 2]
      #   rawCC$tpairn <- rawCC$tpairn[!diagInd, , drop=FALSE]
      #   rawCC$rawCCov <- rawCC$rawCCov[!diagInd]
      # }
      # res <- Lwls2D(userBwCov, kern, rawCC[['tpairn']], rawCC[['rawCCov']], 
      #               xout1=outGrid, xout2=outGrid, crosscov=TRUE)
      
    }
    attr(res, 'covType') <- 'FF'
  }
  
  return(res)
}

## Concurrent functional regression by imputation. This does not provide consistent estimates.
## FPCAlist: a list of functional covariates and response. Each field corresponds to a covariate. 
#            The last entry is assumed to be the response if no entry is names 'Y'.
imputeConReg <- function(FPCAlist, Z, outGrid) {
  
  if (is.null(names(FPCAlist)))
    names(FPCAlist) <- c(paste0('X', seq_len(length(FPCAlist) - 1)), 'Y')
  
  if ('Y' %in% names(FPCAlist)) {
    Yname <- 'Y'
    FPCAlist <- c(FPCAlist[names(FPCAlist) != 'Y'], FPCAlist['Y'])
  } else 
    Yname <- names(FPCAlist)[length(FPCAlist)]
  
  imputeCurves <- sapply(FPCAlist, function(x) 
    apply(stats::fitted(x), 1, function(fit) 
      stats::approx(x[['workGrid']], fit, outGrid)[['y']]),
    simplify='array')
  alphaBeta <- apply(imputeCurves, 1, function(XYt) {
    Yt <- XYt[, ncol(XYt)]
    designMat <- cbind(1, XYt[, -ncol(XYt), drop=FALSE], Z)
    beta_t <- qr.solve(designMat, Yt)
    return(beta_t)
  })
  beta0 <- alphaBeta[1, ]
  beta <- alphaBeta[-1, , drop=FALSE]
  
  return(list(beta0 = beta0, beta = beta, outGrid = outGrid))
}

## subset a list of covariates and responses.
subsetVars <- function(vars, subset) {
  sapply(vars, function(x) {
    if (is.list(x)) {
      sapply(x, `[`, subset, drop=FALSE, simplify=FALSE)
    } else if (is.numeric(x)) {
      x[subset, drop=FALSE]
    } else {
      stop('Cannot subset variable')
    }
  }, simplify=FALSE)
}

## get the number of subjects for a list of covariates and responses.
lengthVars <- function(vars, subset) {
  lenEach <- sapply(vars, function(x) {
    if (is.list(x)) {
      sapply(x, length)
    } else if (is.numeric(x)) {
      length(x)
    } else {
      stop('Cannot subset variable')
    }
  }, simplify=FALSE)
  len <- unique(unlist(lenEach))
  if (length(len) != 1) {
    stop('Length of variables are not the same!')
  }
  
  return(len)
}






