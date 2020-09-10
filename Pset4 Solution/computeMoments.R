computeMoments <- function(param,data){
  
  
  oout <- computeMoments2(data)
  Ymm <- oout[[1]]
  Xmm <- oout[[2]]
  
  # param 
  # sigma_zeta^2, sigma_zeta eta, sigma m^2
  sigmazeta2   <- param[1]
  sigmazetaeta <- param[2]
  sigmam2      <- param[3]
  
  diff <- c( Ymm[1] - sigmazetaeta*Xmm[1],
             Ymm[2] - sigmazeta2 - (sigmazetaeta^2)*Xmm[2],
             Ymm[3] - sigmazeta2 - (sigmazetaeta^2)*Xmm[3] - 2*sigmam2)
  
  dim(diff) <-c(3,1)
  
  # return moments
  return(diff)
}
