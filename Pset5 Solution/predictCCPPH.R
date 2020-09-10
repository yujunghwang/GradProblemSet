predictCCPPH <- function(param,p0,p1){
  
  theta1 <- param[1]
  theta2 <- param[2]
  
  # likelihood function
  # v1 - v0
  vdiff <- -theta2  + beta*(theta1*txs) + crrn
  p0hat <- 1/(1+exp(vdiff))
  p1hat <- 1-p0hat
  
  return(list(p0hat,p1hat))
  
}
