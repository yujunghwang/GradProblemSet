findMomentVec <- function(paramForOpt){
  
  beta  <- paramForOpt[1]
  gamma <- paramForOpt[2]
  
  # Solve and simulate
  if (solveUsingValueFunction ==1){
    oout <- solveValueFunctionBCIncEst(beta,gamma)
  } else if (solveUsingEulerEquation ==1) {
    oout <- solveEulerEquationBCIncEst(beta,gamma)
  }
  policyA1 <- oout[[1]]
  policyC  <- oout[[2]]
  val      <- oout[[3]]
  exVal    <- oout[[4]]
  exDu     <- oout[[5]]
  
  #-----------------------#
  # SIMULATE CONSUMER'S PATH
  # start from initial level of assets and simulate optimal consumption and
  # savings profiles over lifecycle
  
  if (isUncertainty == 0){
    oout2 <- simNoUncerInc(policyA1,exVal,startA)
  } else {
    # Draw random draws for starting income and for innovations
    # It is important to set seed, otherwise estimate will keep changing whenever you run estimation
    
    seed1 <- 201827482
    seed2 <- 182993711
    
    oout2 <- simWithUncerInc(policyA1,exVal,startA,seed1,seed2)
    
  }
  ypath <- oout2[[1]]
  cpath <- oout2[[2]]
  apath <- oout2[[3]]
  vpath <- oout2[[4]]
  
  # truncate apath last period T+1
  apath <- apath[1:T,]
  
  # Get moments in the simulations 
  smomentA  <- apply(apath,1,mean) # apath dim 1 : time, 2 : num of Sims

  sbetaA    <- rep(NA,T)
  sbetaY    <- rep(NA,T)
  for (i in 1:T){
    oout <- lm(cpath[i,]~apath[i,]+ypath[i,])
    sbetaA[i] <- coef(oout)[2]
    sbetaY[i] <- coef(oout)[3]
  }
  
  simMoments    <- c(smomentA[2:T],sbetaY[1:(T-1)],sbetaA[2:(T-1)])
    
  # Get difference between sims and dat using the estimated Param
  diffEst <- simMoments - dataMoments
  
  return(diffEst)
}
