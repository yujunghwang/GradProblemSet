getPolicy_analytical <- function(){
  
  # This function gets analytical policy functions in a world where
  # a) there is no uncertainty and 
  # b) borrowing is allowed up to the NBL
  
  # --------
  # initialize
  policyC <- matrix(rep(NA,T*numPointsA),ncol=numPointsA)
  policyA1<- matrix(rep(NA,T*numPointsA),ncol=numPointsA)
  
  # the following is best read with the notes
  alpha <- (beta^(1/gamma))*((1+r)^((1-gamma)/gamma))
  
  for (ixt in 1:T){
    periodsLeft <- T-ixt+1
    for (ixA in 1:numPointsA){
        if (abs(alpha-1)<1e-5){
          policyC[ixt,ixA] <- Agrid[ixt,ixA]/periodsLeft
        } else{
          policyC[ixt,ixA] <- ((1-alpha) / (1-(alpha^periodsLeft))) * Agrid[ixt,ixA]
        }
      policyA1[ixt,ixA] <- (1+r)*(Agrid[ixt,ixA] - policyC[ixt,ixA])
    }
  }
  return(list(policyA1,policyC))
}