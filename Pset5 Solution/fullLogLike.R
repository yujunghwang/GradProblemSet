# likelihood based on full solution method
fullLogLike <- function(param){
  
  theta1 <- param[1]
  theta2 <- param[2]

  # value function iterations
  xgrid <- seq(0,M,1) # mileage grid from 0 to M
  v0    <- rep(1,M+1) # conditional value function for no replacement
  v1    <- rep(1,M+1) # conditional value function for replacement
  
  # convergence criteria
  tol <- 10^-5
  itermax <- 1000
  dif  <- 1
  iter <-1
  
  while (iter<=itermax & dif >tol){
    v0new <- -theta1*xgrid          + beta*log(exp(c(v0[2:(M+1)],v0[M+1])) + exp(c(v1[2:(M+1)],v1[M+1])) ) + beta*euler
    v1new <- -theta1*xgrid - theta2 + beta*log(exp(rep(v0[1],M+1))       + exp(rep(v1[1],M+1)) )           + beta*euler
    
    dif <- max(c(abs(v0new-v0),abs(v1new-v1)))
    iter <- iter+1
    v0 <- v0new
    v1 <- v1new
  }
  
  # policy function
  p0 <- exp(v0)/(exp(v0)+exp(v1))
  p1 <- 1-p0
  
  # compute likelihood
  loglike <- rep(NA,N*T)
  for (k in 1:(N*T)){
    loglike[k] <- (1-d[k])*log(p0[x[k]+1]) + d[k]*log(p1[x[k]+1])
  }

  # return
  return(loglike)  
}

