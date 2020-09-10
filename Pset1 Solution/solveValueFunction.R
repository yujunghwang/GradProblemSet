utility <-function(cons){
  # this function takes consumption and returns utility
  
  if (cons <=0){
    print('error in utility. consumption is <=0')
  }
  if (gamma==1){
    utils <- log(cons)
  } else {
    utils <- ((cons)^(1-gamma))/(1-gamma)
  }
  return(utils)
}
objectivefunc <- function(A1,A0,Agrid1,V1){
 
  # this function returns
  # -(u(c) + bV(A1))
  # where c is calculated from today's assets and tomorrow's assets
  
  # get tomorrow's consumption (cons), the value of left over assets (VA1) and
  # total value (u(c) + b*VA1)
  
  cons <- A0 - A1/(1+r)
  VA1 <- interp1(Agrid1,V1,A1,method=interpMethod)
  value <- utility(cons) + beta*VA1
  
  # take minus so minimization finds maximum
  return(-value)
}

solveValueFunction <- function(){
  
  #This function obtains the value function for each time period and
  #the policy function (i.e. optimal next-period asset choice) for each time
  #period. From there we can work the optimal consumption level.
  
  #The approach taken is by backwards recursion. The optimisation each period
  #is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
  #The optimisation routine it uses is known as the 'golden search method'
  
  # ------------------------------------------------------------------------ 
  # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN
  
  # Matrices to hold the policy and value functions 
  V        <- matrix(rep(NA, (T+1)*numPointsA), ncol=numPointsA)
  policyA1 <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  policyC  <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  dU       <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  
  #------------------------------------------------------------------------ 
  #Set the terminal value function to 0
  V[(T+1),]<-0
  
  ## ------------------------------------------------------------------------ 
  # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
  # BACKWARDS TO ZERO, ONE PERIOD AT A TIME
  
  for (ixt in T:1){ # loop from time T-1 to 1
  
    V1 <- V[(ixt+1),] # tomorrow's value function
    Agrid1 <- Agrid[(ixt+1),] # tomorrow's asset
    
    for (ixA in 1:numPointsA){ # points on asset grid
      
      # solve problem at asset grid point
      # information on optimisation
      
      A <- Agrid[ixt,ixA] # asset today
      lbA1 <- Agrid1[1] # lower bound : assets tomorrow
      ubA1 <- (A-minCons)*(1+r) # upper bound : assets tomorrow
      
      # compute solution
      if (ubA1 - lbA1 < minCons){ # if liquidity constrained
        negV <- objectivefunc(lbA1, A, Agrid1,V1)
        policyA1[ixt,ixA] <- lbA1
      } else{ # if interior solution
        oout<- fminbnd(f=objectivefunc,A0=A,Agrid1=Agrid1,V1=V1,a=lbA1,b=ubA1)    
        policyA1[ixt,ixA]  <- oout$xmin
        negV <- oout$fmin
      }
      
      # store solution and its value
      policyC[ixt,ixA] <- A - policyA1[ixt,ixA]/(1+r)
      dU[ixt,ixA] <- getmargutility(policyC[ixt,ixA])
      V[ixt,ixA]  <- -negV
    }
  }
  
  return(list(policyA1, policyC, V, dU))
}

