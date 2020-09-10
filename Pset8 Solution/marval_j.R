marval_j <- function(lambda){
# married j's value function
  
  # compute couple's ccp of work choices
  pr_ij <- ccp_ij(lambda)
  dim(lambda) <- c(I,J)
  
  # compute Women's value 
  V1_j <- array(0,c(4,I,J))
  for (k in 1:4){
    V1_j[k,,] <- (log( (1-lambda)*c1[k,,]) -beta*work1_j[k])
  }
  V1_j <- apply(V1_j*pr_ij, c(2,3), sum)
    
  return(V1_j)
}
