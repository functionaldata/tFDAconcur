# #‘ \cite{This function is copied from the fdapace package}

CVLwls1D <- function(y, t, kernel, npoly, nder, dataType, kFolds = 5, useBW1SE = FALSE ){

  # If 'y' and 't' are vectors "cheat" and break them in a list of 10 elements
  if ( is.vector(y) && is.vector(t) && !is.list(t) && !is.list(y) ){
    if (length(t) < 21) {
      stop("You are trying to use a local linear weight smoother in a vector with less than 21 values.\n")
    }
    myPartition =   c(1:10, sample(10, length(t)-10, replace=TRUE));
    y = split(y, myPartition)
    t = split(t, myPartition)
    dataType = 'Sparse';
  } 
  
  # Make everything into vectors
  ncohort = length(t);
  tt  = unlist(t);
  yy  = unlist(y);
  ind = unlist(lapply( 1:ncohort, function(j) rep(j, times=length(t[[j]]))));
  yyn = yy[order(tt)];
  ind = ind[order(tt)];
  ttn = sort(tt);
  
  # Get minimum reasonable bandwidth
  a0=ttn[1];
  b0=ttn[length(ttn)];
  rang = b0-a0;
  dstar = Minb(tt, npoly+2);
  if (dataType != 'Dense'){
   h0 = 2.5*dstar;
  } else {
   h0 = dstar;
  }
  if (h0 > rang/4){
    h0 = h0*.75;
    warning(sprintf("Warning: the min bandwith choice is too big, reduce to %f !", (h0)  ))
  }    
  
  # Get the candidate bandwidths
  nbw = 11;
  bw = rep(0,nbw-1);
  n = length(unique(tt));
  for (i in 1:(nbw-1)){
    bw[i]=2.5*rang/n*(n/2.5)^((i-1)/(nbw-1)); # Straight from MATLAB
  }
  bw = bw-min(bw)+h0;
  
  #ave = rep(0, length(t[[1]]));
  #
  #if (dataType == 'Dense'){
  #  for (i in 1:ncohort){
  #    ave = ave + t[[i]]/ncohort;
  #  }
  #}

  cv = matrix(0, ncol = length(bw), nrow = kFolds);
  #count = c();
  theFolds =  CreateFolds(unique(ind), k= kFolds)

  for (j in 1:(nbw-1)){
    # cv[j]=0;
    # count[j]=0;
    #for (i in 1:ncohort){
    for (i in 1:kFolds){
      
      xout= ttn[ ind %in% theFolds[[i]]];
      obs = yyn[ ind %in% theFolds[[i]]];
      xin = ttn[!ind %in% theFolds[[i]]];
      yin = yyn[!ind %in% theFolds[[i]]];
      
      win=rep(1,length(yin));
      #win[ind==i] = NA;        
      #if(dataType=='Dense') {
      #  yyn=(ave*ncohort-t[[i]])/(ncohort-1);
      #  ttn=t[[1]];
      #  win=pracma::ones(1,length(t[[1]]));    
      #  yyn = yyn[order(ttn)]
      #  ttn = sort(ttn)           
      #}  
 
      mu = tryCatch(
        fdapace::Lwls1D(bw= bw[j], kernel_type = kernel, npoly=npoly, nder= nder, xin = xin, yin= yin, xout=xout, win = win), 
        error=function(err) {
        warning('Invalid bandwidth during CV. Try enlarging the window size.')
        return(Inf)
      })
        
      cv[i,j] = sum((obs-mu)^2)
      # print(cv)
      if(is.na(cv[i,j])){
        cv[i,j] = Inf;
      }
      #count[j] = count[j]+1;
    }
  }
  #cv = cv[(count/ncohort>0.90)];
  #bw = bw[(count/ncohort>0.90)];
  if(min(cv) == Inf){
    stop("All bandwidths resulted in infinite CV costs.")
  }
  if( useBW1SE ){
    # This will pick the bandwidth that is the max but it's average cost is at most
    # 1 standard error of the minimum cost /  I use means so it is more straighforward what the SE is.
    bopt = bw[max(which( 
      colMeans(cv) < min(colMeans(cv)) + apply(cv,2, stats::sd)[which.min(colMeans(cv))]/sqrt(kFolds)))]
  } else {
    bopt = bw[which.min( colMeans(cv))];
  }
  
  return(bopt)

}

CreateFolds <- function(y, k=10) {
  # Returns the test sets indices for k-fold cross-validation. Stratified sampling is used.
  # y: the response vector, for stratifying
  # k: number of folds
  # returns: a list of length k, containing the test-set indices.
  n <- length(y)
  if (n == 0)
    stop('response length is zero')
  
  uniqY <- unique(y)
  if (!is.factor(y) && length(y) / length(uniqY) >= k) {
    # Intepret the integer-valued y as class labels. Stratify if the number of class labels is <= 5.
    y <- factor(y)
  } else if (is.numeric(y)) { 
    # 5-stratum Stratified sampling
    if (n >= 5 * k) {
      breaks <- unique(stats::quantile(y, probs=seq(0, 1, length.out=5)))
      y <- as.integer(cut(y, breaks, include.lowest=TRUE))
    } else 
      y <- rep(1, length(y))
  }
  
  sampList <- tapply(seq_along(y), y, SimpleFolds, k=k, simplify=FALSE)
  list0 <- list()
  length(list0) <- k
  samp <- Reduce(function(list1, list2) {
    mapply(c, list1, list2, SIMPLIFY=FALSE)
  }, sampList, list0)
  
  return(samp)
}


SimpleFolds <- function(yy, k=10) {
  # Simple k-fold test-set samples.
  # Input a set of SAMPLES
  # Returns: a list of length k, containing the SAMPLES.
  if (length(yy) > 1)
    allSamp <- sample(yy)
  else
    allSamp <- yy
  
  n <- length(yy)
  nEach <- n %/% k
  samp <- list()
  length(samp) <- k
  for (i in seq_along(samp)) {
    if (nEach > 0)
      samp[[i]] <- allSamp[1:nEach + (i - 1) * nEach]
    else
      samp[[i]] <- numeric(0)
  }
  restSamp <- allSamp[seq(nEach * k + 1, length(allSamp), length.out=length(allSamp) - nEach * k)]
  restInd <- sample(k, length(restSamp))
  for (i in seq_along(restInd)) {
    sampInd <- restInd[i]
    samp[[sampInd]] <- c(samp[[sampInd]], restSamp[i])
  }
  
  return(samp)
}
