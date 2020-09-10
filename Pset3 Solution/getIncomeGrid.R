stdnormpdf_manual <- function(x){
  # This function gives the pdf of a standard normal
  pdf <- ((sqrt(2*pi))^-1)*(exp(-((x)^2)/2))
  return(pdf)
}

stdnormcdf_manual <- function(x){
  # this function gives the cdf of a standard normal
  approxforminusinf <- -20
  # Numerical integration
  cdf <- quad(stdnormpdf_manual,approxforminusinf,x,tol=1e-12)
  return(cdf)
}
cdfforzero <- function(a,p){
  return(stdnormcdf_manual(a)-p)
}

stdnorminv_manual <- function(p){
  # this routines gives the standard normal value associated with a particular p value
  # with truncations at -3 and 3
  if (p<stdnormcdf_manual(-3)){
    x <- -3
  }else if (p>stdnormcdf_manual(3)) {
    x <- 3
  } else{
    boundforzero <- c(-3,3)
    x <- fzero(fun=cdfforzero, p=p, x=boundforzero, tol=1e-12)$x
  }
  return(x)
}

getNormDev <- function(mu, sigma_inc, trunc, N){
  # This function returns two vectors :
  # a) Z - a vector of (N+1) normal deviates (numbers) that divide a normal
  # distribution with standard deviation sigma into N segments,
  # each of which are equiprobable
  # b) EVbetweenZ - the expected value of the random variables between those two points
  
  #-----------------
  # Initialize the output nad working arrays
  
  # output
  Z <- rep(NA,N+1)
  EVbetweenZ <- rep(NA,N)
  
  #------------------
  # Finding the points that divide the standard normal into N segments
  # the first and last of these should, if we are using an actual normal
  # distribution would be minus and plus infinity.
  # In reality we use a truncated normal distribution and use -trunc and -trunc 
  # where trunc is an arument
  Z[1]   <- -trunc*sigma_inc
  Z[N+1] <-  trunc*sigma_inc
  
  # Now recursively get the rest of the points
  for (ixi in 2:N){
    Z[ixi] <- sigma_inc * stdnorminv_manual((ixi-1)/N)+mu
  }
  stdZ <- (Z-mu*rep(1,N+1))/sigma_inc
  
  # Finding the expected value within each interval (see Adda & Cooper page 58)
  PDF <- stdnormpdf_manual(stdZ)
  for (ixi in 1:N){
    EVbetweenZ[ixi] <- N*sigma_inc*(PDF[ixi]-PDF[ixi+1]) + mu
  }
  return(list(Z,EVbetweenZ))
}

getIncomeGrid <- function(){
  # this function returns
  # an income grid
  # a Markovian transition matrix Q over income realization
  # a vector of minimum incomes in each year
  # a vector of maximum income in each year
  
  # Scenario where there is no uncertainty
  if (isUncertainty==0){
    # income is set equal to the exp of the log mean
    y <- exp(mu)
    minInc <- y
    maxInc <- y
    Q <- 1
    # transition matrix Q is simply a constant 1
    # with prob 1 each period income is 1
  } else if (isUncertainty==1){
    # first get the standard deviation of income (From sigma and rho)
    sig_inc <- sigma/((1-rho^2)^0.5)
    
    # Split the entire normal distribution into numPointsY sections that 
    # are equiprobable. The output lNormDev gives the (numPointsY + 1)
    # points that bound the sections, the output ly gives the
    # (numPointsY) expected value in each section
    
    oout <- getNormDev(mu, sig_inc, normBnd, numPointsY)
    lNormDev <- oout[[1]]
    ly <- oout[[2]]
    
    #-------------------#
    # Get transition matrix Q[i,j]. The prob of income j in t+1
    # conditional on income i in t
    
    Q <- matrix(rep(NA,numPointsY^2),ncol=numPointsY)
    for (i in 1:numPointsY){
      for (j in 1:numPointsY){
        hiDraw <- lNormDev[j+1] - (1-rho)*mu - rho*ly[i] # highest innovation that will give us income tmrw
        loDraw <- lNormDev[j]   - (1-rho)*mu - rho*ly[i] # lowest innovation that will give us income j tmrw
        Q[i,j] <- stdnormcdf_manual(hiDraw/sigma) - stdnormcdf_manual(loDraw/sigma)
      } # j
      
      # Each of the rows of Q should add up to 1. But
      # due to the truncation of the normal distribution
      # they will not. So we divide through by the sum of elements in the row
      Q[i,] <- Q[i,]/sum(Q[i,])
    } # i
    
    
    y <- exp(ly) # get y from log y
    minInc <- exp(-normBnd*sig_inc) # get the minimum income in each year
    maxInc <- exp( normBnd*sig_inc) # get the maximum income in each year
    
    if ((y[1]<1e-4) | (y[numPointsY] > 1e5)){
      print('Combination of sigma and rho give a very high income variance. Numerical instability possible')
    }
    
  } # if isUncertainty==0
  
  #----------------------------#
  # Now get a matrix, T * numPointsY that holds the grid for each income in
  # each year. Do likeiwse with minimum and maximum income
  
  Ygrid  <- t(matrix(rep(y,T),ncol=T)) # T by length(y)
  minInc <- t(matrix(rep(minInc,T),ncol=T))  
  maxInc <- t(matrix(rep(maxInc,T),ncol=T))
  
  #------------------------------#
  # replace these arrays with zeros for all years after retirements
  
  if (Tretire==0){
    # no work (retired at birth)
    Ygrid  <- matrix(rep(0,T*numPointsY),ncol=numPointsY)
    minInc <- matrix(rep(0,T*numPointsY),ncol=numPointsY) 
    maxInc <- matrix(rep(0,T*numPointsY),ncol=numPointsY)
  } else if ( (Tretire>0) & (Tretire<=T) ) { # retire at some age
    Ygrid[Tretire:T,]  <-0
    minInc[Tretire:T,] <-0
    maxInc[Tretire:T,] <-0
  }

  
  return(list(Ygrid,Q,minInc,maxInc))
}


