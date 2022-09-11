#' @title Functional History Index Model
#' @description Functional history index model for functional responses and functional predictors.
#' @param Y A list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., list(Ly = LY,Lt = LT)).
#' @param X A list of lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}.
#' @param Lag A length \code{length(X)} vector denoting the lags for all predictors.
#' @param optnsY A list of options control parameters for the response specified by \code{list(name=value)}. See `Details' in  \code{FPCA}.
#' @param optnsX A list of options control parameters for the predictors specified by \code{list(name=value)}. See `Details' in  \code{FPCA}.
#' @details The functional history index model is defined as 
#' \eqn{E[Y(t)|X_1(t), \cdots, X_p(t)] = \beta_0(t) + \sum_{j=1}^p\beta_j(t)\int_0^{\Delta_j}\gamma_j(s)X_j(t-s)ds} 
#' for \eqn{t\in[\max_j\{\Delta_j\}, T]} with a suitable \eqn{T>0}. 
#' Write \eqn{\alpha_j(t, s)=\beta_j(t)\gamma_j(s)}. It becomes 
#' \eqn{E[Y(t)|X_1(t), \cdots, X_p(t)] = \beta_0(t) + \sum_{j=1}^p\int_0^{\Delta_j}\alpha_j(t, s)X_j(t-s)ds}. 
#' For more details we refer to 
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' @return A list of the following:
#' \item{beta0}{a vector of \code{length(workGridY)} representing the fitted \eqn{\beta_0(t)}.}
#' \item{beta1}{a list of vectors with the \eqn{j}-th element representing the fitted \eqn{\beta_j(t)}.}
#' \item{gamma}{a list of vectors with the \eqn{j}-th element representing the fitted \eqn{\gamma_j(s)}.}
#' \item{workGridY}{a vetor representing the working grid for the response.}
#' \item{workGridLag}{a list of vectors with the \eqn{j}-th element representing the working grid in \eqn{[0, \Delta_j]}.}
#' @examples 
#' set.seed(1)
#' ### functional covariate X(t) ###
#' phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
#' phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
#' lambdaX <- c(10, 5)
#' n <- 150
#' N <- 101
#' Xi <- matrix(rnorm(2*n), nrow = n, ncol = 2)
#' denseLt <- list()
#' denseLy <- list()
#' t0 <- seq(0, 15, length.out = N)
#' for (i in 1:n) {
#'   denseLt[[i]] <- t0
#'   denseLy[[i]] <- lambdaX[1]*Xi[i, 1]*phi1(t0) + lambdaX[2]*Xi[i, 2]*phi2(t0)
#' }
#' denseX0 <- list(Ly = denseLy,Lt = denseLt)
#'
#' ### generate coefficient function gamma(u), beta(u) ###
#' Lag <- 5
#' u0 <- t0[t0<=Lag]
#' t0_out <- t0[t0>=Lag]
#' gamma_u <- function(u) sqrt(2/5) * cos(pi * u /5)
#' beta_1 <- function(t) 5*sin(pi*t/10)
#' beta_0 <- function(t) t^2/2
#'
#' ### functional response Y(t), t in t0_out ### 
#' denseLt <- list()
#' denseLy <- list()
#' for (i in 1:n) {
#'   denseLt[[i]] <- t0_out
#'   Xt <- denseX0$Ly[[i]]
#'   Xtu <- t(sapply((1:N)[t0>=Lag], function(j){
#'     rev(Xt[(j-length(u0)+1):j])  #history index for X[t-u:t]
#'   }))
#'   IntGammaXtu <- apply(Xtu, 1, function(v){
#'     fdapace::trapzRcpp(u0, gamma_u(u0) * v)
#'   })
#'   #append 0 in the first length(u0)-1 element(useless info. in our modeling)
#'   denseLy[[i]] <- beta_0(t0_out) + IntGammaXtu * beta_1(t0_out) + rnorm(length(t0_out), 0, 0.1) 
#' }
#' denseY <- list(Ly = denseLy, Lt = denseLt)
#'
#' ### functional predictor X(t) (adjust for t0_out) ###
#' denseLt <- list()
#' denseLy <- list()
#' for (i in 1:n){
#'   denseLt[[i]] <- t0_out
#'   denseLy[[i]] <- denseX0$Ly[[i]][t0>=Lag]
#' }
#' denseX <- list(Ly = denseLy,
#'                Lt = denseLt)
#'                
#' ### sparify the dense data ###
#' SparseY <- fdapace::Sparsify(samp = t(sapply(denseY$Ly, function(v) v)), 
#'                     pts = t0_out,
#'                     sparsity = sample(10:20, n, replace = TRUE)
#' )
#' SparseX <- fdapace::Sparsify(samp = t(sapply(denseX$Ly, function(v) v)),
#'                     pts = t0_out,
#'                     sparsity = sample(10:20, n, replace = TRUE))
#'
#' ### model fitting for sparse case ###
#' Fit_result <- historyIndexSparse(Y = SparseY, X = list(X = SparseX), Lag = Lag)
#' Fit_result$beta0
#' Fit_result$beta1
#' Fit_result$gamma
#' @references
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' \cite{Yao, F., Müller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. Annals of Statistics 33, 2873--2903.}
#' \cite{Hall, P., Horowitz, J.L. (2007). Methodology and convergence rates for functional linear regression. The Annals of Statistics, 35(1), 70--91.}
#' @export

historyIndexSparse <- function(Y, X, Lag = NULL, optnsY = NULL, optnsX = NULL){
  # require(fdapace)
  # require(MASS)
  d0 = length(X)
  
  #time range
  timeRange = range(unlist(Y$Lt)) 
  for(j in 1:d0){
    timeRange[1] <- max(timeRange[1], min(unlist(X[[j]]$Lt)))
    timeRange[2] <- min(timeRange[2], max(unlist(X[[j]]$Lt)))
  }# common time range for response and predictors
  
  if(is.null(Lag)) stop('please provide Lag with appropriate value')
  maxLag <- max(Lag)
  if(maxLag>=diff(timeRange)){
    stop('Lag should be smaller than the length of the time interval')
  }
  if(!is.null(optnsX) & !is.null(optnsY)){
    if(optnsX$nRegGrid != optnsY$nRegGrid){
      stop("must have the same output Grid")
    }
  }
  
  #prune X
  buff = .Machine$double.eps * max(abs(timeRange)) * 10 #correct float inaccuracy
  for(j in 1:d0){
    XjFilter <- lapply(X[[j]]$Lt, function(Lti) Lti>=timeRange[1]-buff & Lti<=timeRange[2]+buff)
    X[[j]]$Lt <- lapply(1:length(X[[j]]$Lt), function(i) X[[j]]$Lt[[i]][XjFilter[[i]]])
    X[[j]]$Ly <- lapply(1:length(X[[j]]$Ly), function(i) X[[j]]$Ly[[i]][XjFilter[[i]]])
  }
  
  # filter response not in timeRange
  YFilter <- lapply(Y$Lt, function(Lti) Lti>=timeRange[1]-buff & Lti<=timeRange[2]+buff)
  Y$Lt <- lapply(1:length(Y$Lt), function(i) Y$Lt[[i]][YFilter[[i]]])
  Y$Ly <- lapply(1:length(Y$Ly), function(i) Y$Ly[[i]][YFilter[[i]]])
  
  ### FPCA for Y ###
  YFPCA <- fdapace::FPCA(Y$Ly, Y$Lt, optns = optnsY)
  Gyy = YFPCA$smoothedCov
  workGridY = YFPCA$workGrid
  
  #initialize output: gamma and smGrid
  beta0 = YFPCA$mu[workGridY>=maxLag+timeRange[1]-buff]
  beta1 = list()
  tGrid = workGridY[workGridY>=(maxLag+timeRange[1]-buff)]
  gamma = list()
  smGrid = list()
  
  for (j in 1:d0){ #calculate gamma(u) for each X_j
    # FPCA for X_j
    XFPCA = fdapace::FPCA(X[[j]]$Ly, X[[j]]$Lt, optns = optnsX)
    Gxx = XFPCA$smoothedCov
    muX = XFPCA$mu
    workGridX = XFPCA$workGrid
    workGridXAdj = workGridX - min(workGridX)
    #grid for s
    sGrid = workGridXAdj[workGridXAdj <= Lag[j]]
    gammas = numeric(length(sGrid))
    Lj = length(sGrid)
    
    # Gxy  ####(Be careful here!, may need more discussion!)#####
    # Covxy = GetCrCovYX(Ly1 = X[[j]]$Ly, 
    #                      Lt1 = X[[j]]$Lt, 
    #                      Ymu1 = XFPCA$mu, 
    #                      Ly2 = Y$Ly, 
    #                      Lt2 = Y$Lt,
    #                      Ymu2 = YFPCA$mu)
    # Gxy = Covxy$smoothedCC
    # Gxy = t(apply(Gxy, 1, function(v){
    # ConvertSupport(fromGrid = Covxy$smoothGrid[,2], 
    #                toGrid = workGridY,
    #                mu = v)
    #       }))
    Covxy = fdapace::GetCrCovYX(Ly1 = XFPCA$xiEst %*% t(XFPCA$phi),
                       Ly2 = YFPCA$xiEst %*% t(YFPCA$phi))
    Gxy = Covxy$rawCC$rawCCov
    # workGridY
    # alpha tm for each t in the workGridY
    for (i in (1:length(workGridY))[workGridY>=maxLag+timeRange[1]-buff]){
      t = workGridY[i]
      #if (!(t %in% workGridX)) next
      #workIndXs = workGridX >= t-Lag[j]-buff & workGridX <= t+buff
      #workIndexs = (1:length(workGridX))[workIndXs]
      workIndexs  = (i-Lj+1):i #length Lj, start from i-Lj+1 ends at i
      workGridXs = workGridX[workIndexs]
      alphats = numeric(length(workGridXs))
  
      #G(t-s, t-s) 
      Gxxs = Gxx[rev(workIndexs), rev(workIndexs)]
      #obtain eigenfunction and eigenvalues
      EigenGxxs = GetEigenAnalysisResults(Gxxs, sGrid, optns = list(maxK = 10, FVEthreshold = 0.99, verbose = FALSE))
      phis = EigenGxxs$phi
      lambdas = EigenGxxs$lambda
      #obtain alpham and alphats
      Gxys = Gxy[rev(workIndexs), i]
      for (l in 1:length(lambdas)){
        alpham = 1 / lambdas[l] * fdapace::trapzRcpp(sGrid, Gxys * phis[,l])
        alphats = alphats + alpham * phis[,l]
        }
      gammas = gammas  + alphats
      }
    #obtain gamma s
    gammas = gammas / fdapace::trapzRcpp(sGrid, gammas^2) ^ (1/2) * ifelse(gammas[1] < 0, -1, 1)
    gamma[[j]] = gammas
    smGrid[[j]] = sGrid
    
    #obtain coefficient beta(t)
    #Gxx(t-s,t)
    beta_j = numeric(length(tGrid))
    k = 1
    for (i in (1:length(workGridY))[workGridY>=maxLag+timeRange[1]-buff]){
      Gxyt = Gxy[i,i]
      Gxxt = Gxx[rev((i-Lj+1):i), i] 
      beta_j[k] = Gxyt / fdapace::trapzRcpp(sGrid, gammas * Gxxt)
      beta0[k] = beta0[k] - beta_j[k] * fdapace::trapzRcpp(sGrid, gammas * muX[rev((i-Lj+1):i)])
      k = k+1
    }
    beta1[[j]] = beta_j
  }
  return(list(beta0 = beta0, 
              beta1 = beta1,
              gamma = gamma, 
              workGridY = tGrid, 
              workGridLag = smGrid))
  
}