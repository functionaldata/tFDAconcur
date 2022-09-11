#' @title Functional History Index Model
#' @description Functional history index model for dense functional responses and dense functional predictors.
#' @param Y a list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., list(Ly = LY,Lt = LT)).
#' @param X a list of lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}.
#' @param Lag a length \code{length(X)} vector denoting the lags for all predictors.
#' @param optnsY a list of options control parameters for the response specified by \code{list(name=value)}. See `Details' in \code{fdapace::FPCA}.
#' @param optnsX a list of options control parameters for the predictors specified by \code{list(name=value)}. See `Details' in \code{fdapace::FPCA}.
#' @details The functional history index model is defined as 
#' \eqn{E[Y(t)|X_1(t), \cdots, X_p(t)] = \beta_0(t) + \sum_{i=1}^p\beta_i(t)\int_0^{\Delta_i}\gamma_i(s)X_i(t-s)ds} 
#' for \eqn{t\in[\max_i\{\Delta_i\}, T]} with a suitable \eqn{T>0}. 
#' Write \eqn{\alpha_i(t, s)=\beta_i(t)\gamma_i(s)}. It becomes 
#' \eqn{E[Y(t)|X_1(t), \cdots, X_p(t)] = \beta_0(t) + \sum_{i=1}^p\int_0^{\Delta_i}\alpha_i(t, s)X_i(t-s)ds}. 
#' For more details we refer to 
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' @return A list of the following:
#' \item{beta0}{a vector of \code{length(workGridY)} representing the fitted \eqn{\beta_0(t)}}
#' \item{alpha}{a list of matrices with the \eqn{i}-th element representing the fitted \eqn{\alpha_i(t, s)}}
#' \item{yHat}{an n by \code{length(workGridY)} matrix of fitted \eqn{Y_i(t)}'s.}
#' \item{workGridY}{a vetor representing the working grid for the response.}
#' \item{workGridLag}{a list of vectors with the \eqn{i}-th element representing the working grid in \eqn{[0, \Delta_i]}.}
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
#' denseX0 <- list(Ly = denseLy, Lt = denseLt)
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
#' fit <- historyIndexDense(Y = denseY, X = list(X = denseX), Lag = Lag)
#' fit$beta0
#' fit$alpha[[1]]
#' @references 
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' \cite{Yao, F., Müller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. Annals of Statistics 33, 2873--2903.}
#' \cite{Hall, P., Horowitz, J.L. (2007). Methodology and convergence rates for functional linear regression. The Annals of Statistics, 35(1), 70--91.}
#' @export

historyIndexDense <- function(Y, X, Lag = NULL, optnsY = NULL, optnsX = NULL){
  d0 <- length(X)
  timeRange <- range(unlist(Y$Lt))
  for(j in 1:d0){
    timeRange[1] <- max(timeRange[1], min(unlist(X[[j]]$Lt)))
    timeRange[2] <- min(timeRange[2], max(unlist(X[[j]]$Lt)))
  }# common time range for response and predictors
  
  if(is.null(Lag)) stop('please provide Lag with appropriate value')
  maxLag <- max(Lag)
  if(maxLag >= diff(timeRange)){
    stop('Lag should be smaller than the length of the time interval')
  }
  # filter observations not in timeRange
  YFilter <- lapply(Y$Lt, function(Lti) Lti>=timeRange[1] & Lti<=timeRange[2])
  Y$Lt <- lapply(1:length(Y$Lt), function(i) Y$Lt[[i]][YFilter[[i]]])
  Y$Ly <- lapply(1:length(Y$Ly), function(i) Y$Ly[[i]][YFilter[[i]]])
  for(j in 1:d0){
    XjFilter <- lapply(X[[j]]$Lt, function(Lti) Lti>=timeRange[1] & Lti<=timeRange[2])
    X[[j]]$Lt <- lapply(1:length(X[[j]]$Lt), function(i) X[[j]]$Lt[[i]][XjFilter[[i]]])
    X[[j]]$Ly <- lapply(1:length(X[[j]]$Ly), function(i) X[[j]]$Ly[[i]][XjFilter[[i]]])
  }
  
  if(length(unique(sapply(Y$Lt, length)))==1){# if each subject has the same number of observations: N
    LtMatrix <- sapply(Y$Lt, sort)# N by n matrix
    if(any(abs(apply(LtMatrix, 1, diff))>0)){# if any difference found in Lt for any pair of subjects --> irregular design
      YFPCA <- fdapace::FPCA(Y$Ly, Y$Lt, optns = optnsY)
      LyMatrix <- YFPCA$mu + YFPCA$phi %*% t(YFPCA$xiEst)# nWorkGrid by n matrix
      Y$Ly <- lapply(seq_len(ncol(LyMatrix)), function(i) LyMatrix[, i])# convert to a list
      Y$Lt <- rep(list(YFPCA$workGrid), length(Y$Lt))
    }# if the response follows irregular design, invoke FPCA() function to calculate the smoothed observation on workGrid --> regular design
  }
  else{
    YFPCA <- fdapace::FPCA(Y$Ly, Y$Lt, optns = optnsY)
    LyMatrix <- YFPCA$mu + YFPCA$phi %*% t(YFPCA$xiEst)# nWorkGrid by n matrix
    Y$Ly <- lapply(seq_len(ncol(LyMatrix)), function(i) LyMatrix[, i])# convert to a list
    Y$Lt <- rep(list(YFPCA$workGrid), length(Y$Lt))
  }# similarly invoke FPCA() function
  
  timeRange[1] <- timeRange[1] + maxLag
  workIndicator <- Y$Lt[[1]]>=timeRange[1] & Y$Lt[[1]]<=timeRange[2]
  workIndex <- (1:length(Y$Lt[[1]]))[workIndicator]
  workGridY <- Y$Lt[[1]][workIndicator]
  gamma <- vector('list', length = d0)
  alpha <- rep.int(0, length(workIndex))
  yHat <- matrix(nrow = length(Y$Ly), ncol = length(workIndex))
  for(i in 1:length(workIndex)){
    index <- workIndex[i]
    t <- Y$Lt[[1]][index]
    Yi <- sapply(Y$Ly, function(Lyi) Lyi[index])
    Xi <- vector('list', length = d0)
    for(j in 1:d0){
      XjFilter <- lapply(X[[j]]$Lt, function(Lti) Lti>=t-Lag[j] & Lti<=t)# close interval here
      Xi[[j]]$Lt <- lapply(1:length(X[[j]]$Lt), function(k) X[[j]]$Lt[[k]][XjFilter[[k]]])
      Xi[[j]]$Ly <- lapply(1:length(X[[j]]$Ly), function(k) X[[j]]$Ly[[k]][XjFilter[[k]]])
    }
    flmi <- fdapace::FLM1(Y = Yi, X = Xi, optnsListY = optnsY, optnsListX = optnsX)
    beta <- flmi$betaList
    for(j in 1:d0){
      gamma[[j]] <- rbind(gamma[[j]], rev(beta[[j]]))# different part of the predictors should follow the same design. Otherwise for different t, beta[[j]] may be of different length.
    }
    alpha[i] <- flmi$alpha
    yHat[, i] <- flmi$yHat
  }
  workGridLag <- lapply(flmi$workGridX, function(i) i-t)
  return(list(beta0 = alpha, alpha = gamma, yHat = yHat, workGridY = workGridY, workGridLag = workGridLag))
}