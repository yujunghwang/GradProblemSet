marval_i <- function(lambda){

  # compute couple's ccp of work choices
  pr_ij <- ccp_ij(lambda)
  
  # compute Men's value 
  dim(lambda) <- c(I,J)
  
  V1_i <- array(0,c(4,I,J))
  for (k in 1:4){
    V1_i[k,,] <- (log(lambda*c1[k,,]) -beta*work1_i[k])
  }
  V1_i <- apply(V1_i*pr_ij, c(2,3), sum)
  
  return(V1_i)
}