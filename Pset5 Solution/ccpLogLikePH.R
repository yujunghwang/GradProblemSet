ccpLogLikePH <- function(param){
  
  oout  <- predictCCPPH(param,p0,p1)
  p0hat <- oout[[1]]
  p1hat <- oout[[2]]
  
  loglike <- (1-td)*log(p0hat) + td*log(p1hat)
  loglike <- tqst[,1]*loglike
  
  return(loglike)
}