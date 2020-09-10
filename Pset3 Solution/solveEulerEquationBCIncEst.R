eulerforzeroInc <- function(A1,A0,Y,Agrid1,dU1,beta,gamma){
  
  #----------
  # this function returns the following quantity :
  # u'(c_t) =b(1+r)u'(c_{t+1})
  # This quantity =0 when the Euler equation u'(c_t) =b(1+r)u'(c_{t+1})
  # is satisfied with equality
  
  # get marginal utility
  # get marginal utility at consumption tomorrow.
  if (linearise ==1) {
    # get linearised marginal utility tomorrow
    dU1 <- getinversemargutility(dU1,gamma)
  }
  
  du1AtA1 <- interp1(Agrid1,dU1,A1,method=interpMethod)
  if (linearise==1){
    du1AtA1 <- getmargutility(du1AtA1,gamma)
  }
  
  # check whether tomorrow's expected marginal utility negative. 
  if (du1AtA1<0){
    stop('approximated marginal utility in negative')
  }
  
  # get conusmption today and the required output
  todaycons <- A0 +Y - A1/(1+r)
  euler <- getmargutility(todaycons,gamma) - (beta*(1+r)*du1AtA1)
  
  return(euler)
}


solveEulerEquationBCIncEst <- function(beta,gamma){
  
  # ------------------------------------------------------------------------ 
  # This function obtains the value function for each time period and
  #the policy function (i.e. optimal next-period asset choice) for each time period.
  
  #The approach taken is by backwards recursion. The optimisation each period
  #is carried out using 'fzero'. This function finds the zero of an arbitrary
  #function. We use it to find the value of saving (and therefore
  #consumption) that ensures that the Euler equation holds.
  
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
  # SOLVE THE CONSUMER'S PROBLEM AT TIME T, WHEN SOLUTION IS KNOWN
  # Optimal consumption is equal to assets held as there is no value to keep them for
  # after death
  
  for (ixA in 1:numPointsA){ # points on asset grid
    
    # Step 1. solve problem for each grid point in assets and income today
    #-------------
    
    for (ixY in 1:numPointsY){ # points on income grid
      
      # Value of state variables
      Y <- Ygrid[T,ixY] # income today
      A <- Agrid[T,ixA] # assets today
      
      # Compute and store solution and its value
      policyC[T,ixA,ixY] <- A+Y # optimal consumption
      policyA1[T,ixA,ixY] <- 0  # optimal next period assets
      V[T,ixA,ixY] <- utility(policyC[T,ixA,ixY],gamma) # value of policyC
      dU[T,ixA,ixY]<- getmargutility(policyC[T,ixA,ixY],gamma) # marginal valye of policyC
      
    } # ixY
    
    # Step 2 . integrate out income today conditional on income yesterday to get 
    # EV and edU
    
    realisedV <- V[T,ixA,]
    realiseddU <- dU[T,ixA,]
    for (ixY in 1:numPointsY){ # for each point on the income grid
      EV[T,ixA,ixY]  <- incTransitionMrx[ixY,]%*%realisedV  # continuation value at T-1
      EdU[T,ixA,ixY] <- incTransitionMrx[ixY,]%*%realiseddU  # expect marginal utility at T
    } # ixY
  } # ixA
  
  #----------------------#
  # Solve recursively the consumer's problem, starting at time T-1
  
  for (ixt in (T-1):1){ # Loop from time T-1 to 1
    Agrid1 <- Agrid[ixt+1,] # The grid on assets tmrw
    for (ixA in 1:numPointsA){ # Loop through points on the asset grid
      
      # Step 1. Solve problem for each grid point in assets and income today
      #----------
      
      for (ixY in 1:numPointsY){
        # Value of income and informatino for optimisation
        A <- Agrid[ixt,ixA] # assets today
        Y <- Ygrid[ixt,ixY] # income today
        lbA1 <- Agrid[ixt+1,1] # lower bound : assets tmrw
        ubA1 <- (A+Y-minCons)*(1+r) # upper bound: assets tmrw
        bndForSol <- c(lbA1, ubA1)
        Edu1 <- EdU[ixt+1,,ixY] # relevant section of Edu matrix (in assets tmrw)

        # compute solution
        signoflowerbound <- sign(eulerforzeroInc(A1=lbA1,A0=A,Y=Y,Agrid1=Agrid1,dU1=Edu1,beta=beta,gamma=gamma))
        if ((signoflowerbound==1) | (ubA1 - lbA1 < minCons) ){ # if liquidity constrained
          policyA1[ixt,ixA,ixY] <- lbA1
        } else{ # if interior solution
          signofupperbound <- sign(eulerforzeroInc(A1=ubA1,A0=A,Y=Y,Agrid1=Agrid1,dU1=Edu1,beta=beta,gamma=gamma))
          if (signoflowerbound*signofupperbound==1){
            stop('Sign of lower bound and upper bound are the same - no solution to Euler equation. Bug likely')
          }
          policyA1[ixt,ixA,ixY] <- fzero(eulerforzeroInc,A0=A,Y=Y,Agrid1=Agrid1,dU1=Edu1,beta=beta,gamma=gamma,bndForSol,tol=tol)$x
        }
        
        # Store solution and its value
        policyC[ixt,ixA,ixY] <- A+Y-policyA1[ixt,ixA,ixY]/(1+r)
        dU[ixt,ixA,ixY] <- getmargutility(policyC[ixt,ixA,ixY],gamma)
        EV1 <- EV[ixt+1,,ixY]
        V[ixt,ixA,ixY] <- -objectivefunc(policyA1[ixt,ixA,ixY],A,Y,Agrid1,EV1,beta,gamma)
      } # ixY
      
      # STEP 2. integrate out income today conditional on income
      # yesterday to get EV and EdU
      realisedvalues <- V[ixt,ixA,]
      realisedMargUtility <- dU[ixt,ixA,]
      
      for (ixY in 1:numPointsY){
        EV[ ixt,ixA,ixY] <- incTransitionMrx[ixY,]%*%realisedvalues
        EdU[ixt,ixA,ixY] <- incTransitionMrx[ixY,]%*%realisedMargUtility
      } # ixY
        
    } # ixA
  } # ixt
  
  return(list(policyA1,policyC,V,EV,EdU))
}