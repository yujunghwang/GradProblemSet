utility <-function(cons,gamma){
  # this function takes consumption and returns utility
  
  if (cons <=0){
    stop('error in utility. consumption is <=0')
  }
  if (gamma==1){
    utils <- log(cons)
  } else {
    utils <- ((cons)^(1-gamma))/(1-gamma)
  }
  return(utils)
}

getmargutility <- function(cons,gamma){
  
  if (cons<=0){
    stop('Consumption is <=0')
  }
  
  if (gamma==1){
    margut <- 1/cons
  } else{
    margut <- cons^(-gamma)
  }
  
  return(margut)
}

getinversemargutility <- function(margut,gamma){
  if (gamma==1){
    invmargut <- 1/margut
  } else {
    invmargut <- margut^(-1/gamma)
  }
  return(invmargut)
}

objectivefunc <- function(A1,A0,Y,Agrid1,EV1,beta,gamma){
  
  # this function returns
  # -(u(c) + bV(A1))
  # where c is calculated from today's assets and tomorrow's assets
  
  # get tomorrow's consumption (cons), the value of left over assets (VA1) and
  # total value (u(c) + b*EV1)
  
  cons <- A0 +Y - A1/(1+r)
  VA1 <- interp1(Agrid1,EV1,A1,method=interpMethod)
  value <- utility(cons,gamma) + beta*VA1
  
  
  # take minus so minimization finds maximum
  return(-value)
}

solveValueFunctionBCIncEst <- function(beta,gamma){
  
  # This function obtains the value function for each time period and
  # the policy function (i.e. optimal next-period asset choice) for each time
  # period. From there we can work the optimal consumption level
  
  # The approach taken is by backwards recursion. The optimization each period is 
  # carried out using "fminbnd". 
  # The optimisation routine it uses is known as the "golden search method"
  
  # the following variables will be used as global
  # T, r, tol, minCons, numPointsA, numPointsY, Agrid, Ygrid, incTransitionMrx
  # Agrid1, EV1
  
  #------------------
  # Generate matrices to store numerical approximations and initiate as NA
  
  # matrices to hold the policy, value and marginal utility functions
  V        <- rep(NA, (T+1)*numPointsA*numPointsY)
  policyA1 <- rep(NA, T*numPointsA*numPointsY)
  policyC  <- rep(NA, T*numPointsA*numPointsY)
  dU       <- rep(NA, T*numPointsA*numPointsY)
  
  dim(V)        <- c((T+1),numPointsA,numPointsY)
  dim(policyA1) <- c(T,numPointsA,numPointsY)
  dim(policyC)  <- c(T,numPointsA,numPointsY)
  dim(dU)       <- c(T,numPointsA,numPointsY)
  
  # Matrices to hold expected value and marginal utility functions
  EV        <- rep(NA, (T+1)*numPointsA*numPointsY)
  EdU       <- rep(NA, T*numPointsA*numPointsY)

  dim(EV)  <- c((T+1),numPointsA,numPointsY)
  dim(EdU) <- c(T,numPointsA,numPointsY)
  
  #-----------------------------
  # Set the terminal value function and expected value function to 0
  
  EV[T+1,,] <-0
  V[T+1,,]  <-0
  
  #------------------------------
  # Solve recursively the consumer's problem, starting at time T-1
  # and moving backwrds to zero, one period at a time
  
  for (ixt in T : 1){ # Loop from time T-1 to 1
    Agrid1 <- Agrid[ixt+1,] # The grid on assets tmrw
    
    for (ixA in 1:numPointsA){ # points on asset grid
      
      # Step 1. solve problem at grid points in assets and income
      #-------------
      for (ixY in 1:numPointsY){ # points on income grid
        
        # value of income and information for optimization
        A <- Agrid[ixt,ixA] # assets today
        Y <- Ygrid[ixt,ixY] # income today
        lbA1 <- Agrid[ixt+1,1] # lower bound : assets tmrw
        ubA1 <- (A + Y - minCons)*(1+r) # upper bound : assets tmrw
        EV1 <- EV[ixt+1,,ixY] # relevant section of EV matrix (in assets tmrw)
        
        # Compute solution
        if (ubA1 - lbA1 < minCons){ # if liquidity constrained
           negV <- objectivefunc(A1=lbA1,A0=A,Y=Y,Agrid1=Agrid1,EV1=EV1,beta=beta,gamma=gamma)
           policyA1[ixt,ixA,ixY] <- lbA1
        } else {
          oout <- fminbnd(objectivefunc,A0=A,Y=Y,Agrid1=Agrid1,EV1=EV1,beta=beta,gamma=gamma,lbA1,ubA1,tol=tol)
          policyA1[ixt,ixA,ixY] <- oout$xmin
          negV                  <- oout$fmin
        } # if (ubA1 - lbA1 < minCons)
        
        # Store solution and its value
        policyC[ixt,ixA,ixY] <- A + Y - policyA1[ixt,ixA,ixY]/(1+r)
        V[ ixt,ixA,ixY] <- -negV
        dU[ixt,ixA,ixY] <- getmargutility(policyC[ixt,ixA,ixY],gamma)
        
      } # ixY
      
      # Step 2 : integrate out income today conditional on income yesterday 
      # to get EV and EdU
      
      realisedV  <- V[ ixt,ixA,]
      realiseddU <- dU[ixt,ixA,]
      
      for (ixY in 1:numPointsY){
        EV[ ixt,ixA,ixY] <- sum(incTransitionMrx[ixY,]*realisedV)
        EdU[ixt,ixA,ixY] <- sum(incTransitionMrx[ixY,]*realiseddU)
      } # ixY
      
      
    } # ixA
  } # ixt
  
  
  return(list(policyA1,policyC,V,EV,EdU))  
}
