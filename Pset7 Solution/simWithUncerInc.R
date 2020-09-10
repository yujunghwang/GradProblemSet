checkSimExtrap <- function(lba1,y,t){
  # This is a function to check, in selecting next period's asset we haven't
  # selected a value that is less than the borrowing constraint. This could
  # occur for one of two reasons. First it could be that we are extrapolating
  # (i.e. income in the period is larger than the largest income in the grid.
  # If this is the case we set next period's assets to the lowest permissable
  # level. Otherwise it is likely that there is an error - and in this case
  # we cause the programme to stop
  
  if ((y > Ygrid[t,numPointsY]) | (y < Ygrid[t,1]) ){
    a1 <- lba1
  } else {
    stop('Next periods asset is below minimum permissible asset And we are not extrapolating. Likely there is a bug')
  }
  
  return(a1)
}

getNormalDraws <- function(mu,sigma,dim1,dim2,seed){
  
  # this function returns a dim1 * dim2 two array of pseudo random draws from
  # a normal distribution with mean mu and std sigma. 
  
  #-------------------#
  # set the seed
  set.seed(seed)
  
  #-------------------#
  # Draw standard normal draws, and transform them so they come from a 
  # distribution with our desired mean and stdev
  
  StdRandNormal <- randn(dim1,dim2)
  normalDraws   <- mu + sigma*StdRandNormal

  return(normalDraws)
}

simWithUncerInc <- function(policyA1,EV,startingA,seed1,seed2){
  # This function takes the policy functions and value functions 
  # along with starting assets and returns simulated paths of income, consumption
  # assets and value
  
  #-------------------
  # initialize arrays that will hold the paths of income consumption, value and assets
  
  # arguments for output
  y <- matrix(rep(NA,T*numSims),ncol=numSims)
  c <- matrix(rep(NA,T*numSims),ncol=numSims)
  v <- matrix(rep(NA,T*numSims),ncol=numSims)
  a <- matrix(rep(NA,(T+1)*numSims),ncol=numSims)
  
  # Other arrays that will be used below
  e          <- matrix(rep(NA,T*numSims),ncol=numSims)
  logy1      <- matrix(rep(NA,T*numSims),ncol=numSims)
  ly         <- matrix(rep(NA,T*numSims),ncol=numSims)
  ypathIndex <- matrix(rep(NA,T*numSims),ncol=numSims)
  
  #---------------------------#
  # Obtain time series of incomes for our simulated individuals
  #---------------------------#
  sig_inc <- sigma/ ((1-rho^2)^0.5)
  
  e     <- getNormalDraws(0,  sigma,  T,numSims,seed1) # normally distributed random draws
  logy1 <- getNormalDraws(mu, sig_inc,1,numSims,seed2) # a random draw for the initial income
  
  # Get all the incomes, recursively
  for (s in 1:numSims){ # loop through individuals
    ly[1,s] <- min(max(logy1[1,s],-normBnd*sig_inc),normBnd*sig_inc)
    y[ 1,s] <- exp(ly[1,s])
    for (t in 1:T){
      if (t!=T){ # Get next year's income
        ly[t+1,s] <- (1-rho)*mu + rho*ly[t,s] + e[(t+1),s]
        ly[t+1,s] <- min( max(ly[t+1,s],-normBnd*sig_inc), normBnd*sig_inc)
        y[t+1,s]  <- exp(ly[t+1,s])
      } # if (t != T)
      
      if (t>=Tretire){
        y[t,s] <- 0
      } # if (t>=Tretire)
      
    } # t
  } # s
  
  #-----------------------------#
  # Obtain consumption, asset, and value profiles
  #-----------------------------#
  
  for (s in 1:numSims){
    a[1,s] <- startingA
    for (t in 1:T){
      if (t < Tretire){
        tA1 <- policyA1[t,,]
        tV <- EV[t,,]
        
        # interp2() does not allow extrapolation. So do it manually
        if ((a[t,s]<=max(Agrid[t,])) & (a[t,s]>=min(Agrid[t,])) & (y[t,s]<=max(Ygrid[t,])) & (y[t,s]>=min(Ygrid[t,]))){
          a[t+1,s] <- interp2(Agrid[t,],Ygrid[t,],t(tA1),a[t,s],y[t,s],method=interpMethod)
          v[t  ,s] <- interp2(Agrid[t,],Ygrid[t,],t(tV ),a[t,s],y[t,s],method=interpMethod)
        } else{
          # linear interpolation
          xx   <- data.frame(x1=rep(Agrid[t,],length(Ygrid[t,])), x2=kron(Ygrid[t,],rep(1,length(Agrid[t,]))))
          yy1  <- tA1[1:length(tA1)]
          
          xx$x3 <- xx$x1*xx$x2
          xx$x4 <- (xx$x1)^2
          xx$x5 <- (xx$x2)^2
          xx$x6 <- (xx$x1)^3
          xx$x7 <- (xx$x2)^3
          xx$x8 <- (xx$x4)*xx$x2
          xx$x9 <- (xx$x5)*xx$x1
          
          oout1 <- lm(yy1 ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data=xx)
          nxx  <- data.frame(x1=a[t,s],x2=y[t,s])
          nxx$x3 <-  nxx$x1   *nxx$x2
          nxx$x4 <- (nxx$x1)^2
          nxx$x5 <- (nxx$x2)^2
          nxx$x6 <- (nxx$x1)^3
          nxx$x7 <- (nxx$x2)^3
          nxx$x8 <- (nxx$x4)  *nxx$x2
          nxx$x9 <- (nxx$x5)  *nxx$x1
          
          a[t+1,s] <- predict(oout1,newdata=nxx)
          yy2  <- tV[1:length(tV)]
          oout2 <- lm(yy2 ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9, data=xx)
          v[t,s] <- predict(oout2,newdata=nxx) 
        }
        
      } else{
        tA1 <- policyA1[t,,1]
        tV  <- EV[t,,1]
        a[t+1,s] <- interp(Agrid[t,],tA1,a[t,s],method=interpMethod)
        v[t  ,s] <- interp(Agrid[t,],tV, a[t,s],method=interpMethod)
      } # if (t<Tretire)
      
      # Check whether next period's asset is below the lowest
      # permissible
      
      if ( a[t+1,s] < Agrid[t+1,1]){
         a[t+1,s] <- checkSimExtrap( Agrid[t+1,1], y[t,s], t)
      }
      
      # Get consumption from today's assets, today's income and 
      # tomorrow's optimal assets
      c[t,s] <- a[t,s] + y[t,s] - (a[t+1,s]/(1+r))
      
    } # t
  } # s
  
  return(list(y,c,a,v))
}