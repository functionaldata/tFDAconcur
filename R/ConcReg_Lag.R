#' @title Functional Concurrent Regression with Lag Model
#' @description Functional concurrent regression model with lag for dense functional responses and dense functional predictors.
#' @param Y a list which contains functional responses in the form of a list LY and the time points LT at which they are observed (i.e., list(Ly = LY,Lt = LT)).
#' @param X a list of lists which contains the observed functional predictors list Lxj and the time points list Ltj at which they are observed. It needs to be of the form \code{list(list(Ly = Lx1,Lt = Lxt1),list(Ly = Lx2,Lt = Lxt2),...)}.
#' @param Lag a length \code{length(X)} vector denoting the lags for all predictors.
#' @param optnsY a list of options control parameters for the response specified by \code{list(name=value)}. See `Details' in \code{fdapace::FPCA}.
#' @param optnsX a list of options control parameters for the predictors specified by \code{list(name=value)}. See `Details' in \code{fdapace::FPCA}.
#' @details The functional concurrent regression model with lag is defined as 
#' \deqn{E[Y(t)|X_1(t), \cdots, X_p(t)] = \beta_0(t) + \sum_{j=1}^p\beta_j(t)X_j(t-Lag[j])} 
#' For more details we refer to 
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' @return A list of the following:
#' \item{beta0}{a vector of \code{length(workGridY)} representing the fitted \eqn{\beta_0(t)}}
#' \item{beta}{A matrix for the concurrent regression effects, where rows correspond to different predictors and columns to different time points.}
#' \item{workGridY}{a vetor representing the working grid for the response.}
#' \item{Lag}{a vector representing the lags for all predictors}
#' @examples 
#' phi1 <- function(t) sin(pi*t / 5) / sqrt(5)
#' phi2 <- function(t) cos(pi*t / 5) / sqrt(5)
#' lambdaX <- c(10, 5)
#' n <- 50
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
#' Lag <- c(3)
#' t0_out <- t0[t0>=max(Lag)]
#' beta_1 <- function(t) 5*sin(pi*t/10)
#' beta_0 <- function(t) t^2/2
#' ### functional response Y(t), t in t0_out ### 
#' denseLt <- list()
#' denseLy <- list()
#' for (i in 1:n) {
#'   denseLt[[i]] <- t0_out
#'   denseLy[[i]] <- beta_0(t0_out) +  
#'     denseX0$Ly[[i]][1:length(t0_out)]* beta_1(t0_out) + 
#'     rnorm(length(t0_out), 0, 0.1) 
#' }
#' denseY <- list(Ly = denseLy, Lt = denseLt)
#' model = ConcReg_Lag(Y=denseY, X=list(X1 = denseX0), Lag=Lag)
#' 
#' print(model$workGridY) # workGrid for Y between 6 to 15
#' plot(beta_1(model$workGridY)) # workGrid for X1 between 3 to 12 
#' plot(beta_0(model$workGridY))
#' 
#' plot(model$beta[1,])
#' plot(model$beta0)
#' @references 
#' \cite{Şentürk, D. and Müller, H.G., (2010). Functional varying coefficient models for longitudinal data. Journal of the American Statistical Association, 105(491), pp.1256-1264.}
#' \cite{Yao, F., Müller, H.G., Wang, J.L. (2005). Functional linear regression analysis for longitudinal data. Annals of Statistics 33, 2873--2903.}
#' \cite{Hall, P., Horowitz, J.L. (2007). Methodology and convergence rates for functional linear regression. The Annals of Statistics, 35(1), 70--91.}
#' @export

ConcReg_Lag <- function(Y, X, Lag = NULL, optnsY = NULL, optnsX = NULL){
  d0 <- length(X)
  timeRange <- range(unlist(Y$Lt))
  for(j in 1:d0){
    timeRange[1] <- max(timeRange[1], min(unlist(X[[j]]$Lt)))
    timeRange[2] <- min(timeRange[2], max(unlist(X[[j]]$Lt)))
  }# common time range for response and predictors
  buff = .Machine$double.eps * max(abs(timeRange)) * 10 #correct float inaccuracy
  if(is.null(Lag)) stop('please provide Lag with appropriate value')
  maxLag <- max(Lag)
  if(maxLag >= diff(timeRange)){
    stop('Lag should be smaller than the length of the time interval')
  }
  # filter observations not in timeRange
  YFilter <- lapply(Y$Lt, function(Lti) Lti>=timeRange[1]-buff & Lti<=timeRange[2]+buff)
  Y$Lt <- lapply(1:length(Y$Lt), function(i) Y$Lt[[i]][YFilter[[i]]])
  Y$Ly <- lapply(1:length(Y$Ly), function(i) Y$Ly[[i]][YFilter[[i]]])
  for(j in 1:d0){
    XjFilter <- lapply(X[[j]]$Lt, function(Lti) Lti>=timeRange[1]-buff & Lti<=timeRange[2]+buff)
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
  workIndicator <- Y$Lt[[1]]>=timeRange[1]-buff & Y$Lt[[1]]<=timeRange[2]+buff
  workIndex <- (1:length(Y$Lt[[1]]))[workIndicator]
  workGridY <- Y$Lt[[1]][workIndicator]

  #create dat list 
  dat = list()
  for(j in 1:d0){
    XjFilter <- lapply(X[[j]]$Lt, function(Lti) Lti>=timeRange[1]-Lag[j]-buff & Lti<=timeRange[2]-Lag[j]+buff)# close interval here
    dat[[j]] <- t(sapply(1:length(X[[j]]$Ly), function(k) X[[j]]$Ly[[k]][XjFilter[[k]]]))
  }
  names(dat) = paste0("X", 1:d0)
  dat$Y = t(sapply(1:length(Y$Ly), function(k) Y$Ly[[k]][workIndicator]))
  model = ptFCReg(workGridY, dat)

  return(list(beta0 = model$beta0, 
              beta = model$beta,
              workGridY = model$tGrid,
              Lag = Lag,
              R2 = model$R2
              ))
}


