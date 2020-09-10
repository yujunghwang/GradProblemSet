simNoUncerInc <- function(policyA1,EV,startingA){
  # This function takes the policy functions and value functions in an environment
  # where there is no uncertainty, along with starting assets and returns
  # simulated paths of consumption assets and value
  
  #-------------------
  # initialize arrays that will hold the paths of income consumption, value and assets
  
  # arguments for output
  y <- matrix(rep(NA,T*numSims),ncol=numSims)
  c <- matrix(rep(NA,T*numSims),ncol=numSims)
  v <- matrix(rep(NA,T*numSims),ncol=numSims)
  a <- matrix(rep(NA,(T+1)*numSims),ncol=numSims)
  
  #---------------------------#
  # obtain paths using the initial condition and the policy and value functions
  
  for (s in numSims){
    a[1,s] <- startingA
    for (t in 1:T){
      y[t,s]   <- Ygrid[t]
      v[t,s]   <- interp1(Agrid[t,],EV[t,],a[t,s], method=interpMethod)
      a[t+1,s] <- interp1(Agrid[t,],policyA1[t,],a[t,s],method=interpMethod)
      c[t,s]   <- a[t,s] + y[t,s] - (a[t+1,s]/(1+r))
    } # t
  } # s
  
  return(list(y,c,a,v))
}