ccp_ij <- function(lambda){
  # compute couple's work choice probability
  dim(lambda) <- c(I,J) # I by J matrix
  
  # compute couple's value
  V1_ij <- array(0,c(4,I,J))
  for (k in 1:4){
    for (ixI in 1:I){
      for (ixJ in 1:J){
        V1_ij[k,ixI,ixJ] <-  lambda[ixI,ixJ] *(log(   lambda[ixI,ixJ] *c1[k,ixI,ixJ]) - beta*work1_i[k]  ) + 
                          (1-lambda[ixI,ixJ])*(log((1-lambda[ixI,ixJ])*c1[k,ixI,ixJ]) - beta*work1_j[k]  )
          
      }
    }
  }
  
  ccp_ij <- exp(V1_ij/esd) 
  for (k in 1:4){
    ccp_ij[k,,] <- ccp_ij[k,,] / apply( exp(V1_ij/esd), c(2,3), sum )
  }

  return(ccp_ij)
}