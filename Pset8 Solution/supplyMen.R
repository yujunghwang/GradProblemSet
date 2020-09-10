supplyMen <- function(lambda){

  # input : Pareto weight vector
  # output : measure of men in ij-type marriage
  V1_i <- marval_i(lambda)
  
  supply <- array(0, c(I,J+1))
  # comptue ccp of men
  pr_ij <- exp(cbind(V1_i,V0_i)/tsd) ### 4th column is value of single
  for (k in 1:I){
    pr_ij[k,] <- pr_ij[k,]/sum(pr_ij[k,])
    supply[k,] <- pr_ij[k,]*measure_i[k]
  }
  
  # truncate the measure of single
  return(supply[,1:J])
}
