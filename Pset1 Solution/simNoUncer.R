simNoUncer <- function(policyA1,V,startingA){
  # This function takes the policy functions and value functions in an environment
  # where there is no uncertainty, along with starting assets and returns
  # simulated paths of consumption assets and value
  
  #-------------------
  # initialize arrays that will hold the paths of income consumption, value and assets
  
  # arguments for output
  c <- rep(NA,T)
  v <- rep(NA,T)
  a <- rep(NA,(T+1))
  
  #---------------------------#
  # obtain paths using the initial condition and the policy and value functions
    a[1] <- startingA
    for (t in 1:T){
      v[t]   <- interp1(Agrid[t,],V[t,],a[t], method=interpMethod)
      a[t+1] <- interp1(Agrid[t,],policyA1[t,],a[t],method=interpMethod)
      c[t]   <- a[t] - (a[t+1]/(1+r))
    }
  
  return(list(c,a,v))
}