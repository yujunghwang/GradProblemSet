demandMen <- function(lambda){
  
  # input : Pareto weight vector
  # output : measure of women in ij-type marriage
 V1_j <- marval_j(lambda)
  
 demand <- array(0,c(I+1,J))
  # comptue ccp of women
  pr_ij <- exp(rbind(V1_j,matrix(V0_j,c(1,J)))/tsd) ### 4th row is value of single
  for (k in 1:J){
    pr_ij[,k] <- pr_ij[,k]/sum(pr_ij[,k])
    demand[,k] <- pr_ij[,k]*measure_j[k]
  }
  
  # truncate the measure of single
  return(demand[1:I,])
}
