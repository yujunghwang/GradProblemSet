getmargutility <- function(cons){
  
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

getinversemargutility <- function(margut){
  if (gamma==1){
    invmargut <- 1/margut
  } else {
    invmargut <- margut^(-1/gamma)
  }
  return(invmargut)
}

eulerforzero <- function(A1,A0,Agrid1,dU1){

  #----------
  # this function returns the following quantity :
  # u'(c_t) =b(1+r)u'(c_{t+1})
  # This quantity =0 when the Euler equation u'(c_t) =b(1+r)u'(c_{t+1})
  # is satisfied with equality
  
  # get marginal utility
  # get marginal utility at consumption tomorrow.
  if (linearise ==1) {
    # get linearised marginal utility tomorrow
    dU1 <- getinversemargutility(dU1)
  }
  
  du1AtA1 <- interp1(Agrid1,dU1,A1,method=interpMethod)
  
  if (linearise==1){
    du1AtA1 <- getmargutility(du1AtA1)
  }
  
  # check whether tomorrow's expected marginal utility negative. 
  if (du1AtA1<0){
    stop('approximated marginal utility in negative')
  }
  
  # get conusmption today and the required output
  todaycons <- A0- A1/(1+r)
  euler <- getmargutility(todaycons) - (beta*(1+r)*du1AtA1)
  
  return(euler)
}

solveEulerEquation <- function(){

  # ------------------------------------------------------------------------ 
  #This function obtains the value function for each time period and
  #the policy function (i.e. optimal next-period asset choice) for each time period.
  #The approach taken is by backwards recursion. The optimisation each period
  #is carried out using 'fzero'. This function finds the zero of an arbitrary
  #function. We use it to find the value of saving (and therefore
  #consumption) that ensures that the Euler equation holds.
  
  # ------------------------------------------------------------------------ 
  # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NA
  # Matrices to hold the policy, marginal utility and value functions 
  
  V        <- matrix(rep(NA, (T+1)*numPointsA), ncol=numPointsA)
  policyA1 <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  policyC  <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  dU       <- matrix(rep(NA,     T*numPointsA), ncol=numPointsA)
  
  #------------------------------------------------------------------------ 
  #Set the terminal value function to 0
  V[(T+1),]<-0
  
  # ------------------------------------------------------------------------ 
  #   SOLVE THE CONSUMER'S PROBLEM AT TIME T, WHEN SOLUTION IS KNOWN
  # Optimal consumption is equal to assets held as there is no value to keep them for
  # after death

  
  for (ixA in 1:numPointsA){
    # points on asset grid
    
    # solve problem for each grid point in assets and income today
    #------------------------
    
    # information for optimization
    
    A <- Agrid[T,ixA] # assets today
    
    # compute and store solution and its value
    policyC[T,ixA]  <-A # optimal consumption
    policyA1[T,ixA] <-0 # optimal next period assets
    V[T,ixA]  <- utility(policyC[T,ixA]) # value of policyC
    dU[T,ixA] <- getmargutility(policyC[T,ixA]) # marginal value of policyC 
  }
  
  # ------------------------------------------------------------------------ 
  #   SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
  # BACKWARDS TO ZERO, ONE PERIOD AT A TIME
  
  for (ixt in (T-1):1){
    # loop from time T-1 to 1
    Agrid1 <- Agrid[(ixt+1),] # the grid on assets tmrw
    dU1    <- dU[(ixt+1),]    # relevant section of Edu matrix (in assets tmrw)
    V1     <- V[(ixt+1),]   
    
    for (ixA in 1:numPointsA){
      # loop through points on the asset grid
      
      # Step 1 : solve problem for each grid point in assets today
      
      # value of income and information for optimisation
      A <- Agrid[ixt,ixA]  # assets today
      lbA1 <- Agrid1[1] # lower bound assets tomorrow
      ubA1 <- (A-minCons)*(1+r) # upper bound assets tomorrow
      bndForSol <- c(lbA1,ubA1)
      # if the Euler equation has a solution, it must lie within these bounds

   # compute solution
      signoflowerbound <- sign(eulerforzero(A0=A,A1=lbA1,Agrid1=Agrid1,
                                            dU1=dU1))
      if ((signoflowerbound==1) | (ubA1-lbA1<minCons)) {
        # if liquidity constrained
        policyA1[ixt, ixA] <- lbA1
      } else { # if interior solution
        signofupperbound <- sign(eulerforzero(A0=A,A1=ubA1,Agrid1=Agrid1,
                                              dU1=dU1))
        if (signoflowerbound*signofupperbound==1){
          stop('Sign of lower bound and upper bound are the same -  no solution to Euler equation. Bug likeliy.')
        }

        # find A1 making euler equation zero
        policyA1[ixt,ixA] <- fzero(eulerforzero,A0=A,Agrid1=Agrid1,
                                   dU1=dU1,bndForSol,tol=tol)$x
      }
      
      # Store solution and its value
      policyC[ixt,ixA] <- A-policyA1[ixt,ixA]/(1+r)
      dU[ixt,ixA] <- getmargutility(policyC[ixt,ixA])
      V[ixt,ixA] <- -objectivefunc(policyA1[ixt,ixA],A,Agrid1,V1)
    } # ixA
  }
  
  return(list(policyA1, policyC, V, dU))
}