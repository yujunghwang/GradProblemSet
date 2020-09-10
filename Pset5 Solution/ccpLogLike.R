ccpLogLike <- function(param){
  
  oout  <- predictCCP(param,p0,p1)
  p0hat <- oout[[1]]
  p1hat <- oout[[2]]
  
  loglike <- (1-d)*log(p0hat) + d*log(p1hat)
  
  return(loglike)
}

