utility <-function(cons,d){
  # this function takes consumption and returns utility
  
  if (cons <=0){
    stop('error in utility. consumption is <=0')
  }
  if (gamma==1){
    utils <- log(cons) - delta*as.integer(d==2)
  } else {
    utils <- ((cons)^(1-gamma))/(1-gamma) - delta*as.integer(d==2)
  }
  return(utils)
}

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

objectivefunc <- function(A1,A0,Y,Agrid1,EV1,d){
  
  # this function returns
  # -(u(c) + bV(A1))
  # where c is calculated from today's assets and tomorrow's assets
  
  # get tomorrow's consumption (cons), the value of left over assets (VA1) and
  # total value (u(c) + b*EV1)
  
  cons  <- A0 +Y - A1/(1+r)
  VA1   <- interp1(Agrid1,EV1,A1,method=interpMethod)
  value <- utility(cons,d) + beta*VA1
  
  # take minus so minimization finds maximum
  return(-value)
}

solveDCValueFunctionBCInc <- function(){

  #------------------
  # Generate matrices to store numerical approximations and initiate as NA
  
  # matrices to hold the policy, value and marginal utility functions
  V        <- rep(NA, (T+1)*numPointsA*2)
  policyA1 <- rep(NA, T*numPointsA*2)
  policyC  <- rep(NA, T*numPointsA*2)
  policyD  <- rep(NA, T*numPointsA)
  dU       <- rep(NA, T*numPointsA*2)
  optV     <- rep(NA, (T+1)*numPointsA)
  optdU    <- rep(NA, T*numPointsA)
  
  dim(V)        <- c((T+1),numPointsA,2)
  dim(policyA1) <- c(T,numPointsA,2)
  dim(policyC)  <- c(T,numPointsA,2)
  dim(policyD)  <- c(T,numPointsA)
  dim(dU)       <- c(T,numPointsA,2)
  dim(optV)     <- c((T+1),numPointsA)
  dim(optdU)    <- c(T,numPointsA)
  
  #-----------------------------
  # Set the terminal value function and expected value function to 0
  V[(T+1),,]  <-0
  optV[(T+1),] <-0
  
  #------------------------------
  # Solve recursively the consumer's problem, starting at time T-1
  # and moving backwrds to zero, one period at a time
  
  for (ixt in T : 1){ # Loop from time T-1 to 1
    
    
    Agrid1 <- Agrid[(ixt+1),] # The grid on assets tmrw
    
    for (ixA in 1:numPointsA){ # points on asset grid
      
      # Step 1. solve problem at grid points in assets and income
      #-------------
      # work decision
      for (ixd in 1:2){ 
        
        # value of income and information for optimization
        A <- Agrid[ixt,ixA] # assets today
        lbA1 <- Agrid[(ixt+1),1]
        if (ixd==2){
          ubA1 <- (A + ybar - minCons) # asset tomorrow
          Y <- ybar
        } else {
          ubA1 <- A - minCons
          Y <- 0
        }

        # Compute solution
        if (ubA1 - lbA1 < minCons){ # if liquidity constrained
          if (ixd==2){
           negV <- objectivefunc(A1=lbA1,A0=A,Y=Y,Agrid1=Agrid1,EV1=optV[(ixt+1),],d=ixd)
          } else{
           negV <- objectivefunc(A1=lbA1,A0=A,Y=Y,Agrid1=Agrid1,EV1=V[(ixt+1),,1],d=ixd)
          }
           policyA1[ixt,ixA,ixd] <- lbA1
        } else {
          if (ixd==2){
            oout <- fminbnd(objectivefunc,A0=A,Y=Y,Agrid1=Agrid1,EV1=optV[(ixt+1),],d=ixd,lbA1,ubA1,tol=tol)
          } else{
            oout <- fminbnd(objectivefunc,A0=A,Y=Y,Agrid1=Agrid1,EV1=V[(ixt+1),,1],d=ixd,lbA1,ubA1,tol=tol)
          }
          
          policyA1[ixt,ixA,ixd] <- oout$xmin
          negV                  <- oout$fmin
        } # if (ubA1 - lbA1 < minCons)
        
        # Store solution and its value
        policyC[ixt,ixA,ixd] <- A + Y - policyA1[ixt,ixA,ixd]/(1+r)
        V[ ixt,ixA,ixd] <- -negV
        dU[ixt,ixA,ixd] <- getmargutility(policyC[ixt,ixA,ixd])
        
      } # ixd
    
      # optimal decision and associated marginal utility
      optV[ixt,ixA]       <- max(V[ixt,ixA,1],V[ixt,ixA,2])
      policyD[ixt,ixA]    <- which(V[ixt,ixA,]==max(V[ixt,ixA,]))[1]
      optdU[ixt,ixA]      <- dU[ixt,ixA,policyD[ixt,ixA]]
      
    } # ixA
    
    
  } # ixt
  
  
  return(list(policyA1,policyC,policyD,V,dU,optV,optdU))
}
