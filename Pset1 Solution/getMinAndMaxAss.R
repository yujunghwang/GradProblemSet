getMinAndMaxAss <- function(startA){
  
  # A function that gives the minimum and maximum on hte asset grid in each year.
  # The minimum is the natural borrowing constraint.
  # The maximum is how much one would have if one saved everything
  
  # The following variables are global variables
  # T r minCons
  
  
  # Initialize the outputs
  BC   <- rep(NA, T+1)
  maxA <- rep(NA, T+1)
  
  # Iteratively calculate the borrowing constraints, and maximum on asset grid
  
  # borrowing constraints
  BC[T+1] <-0
  for (ixt in T:1){
    BC[ixt] <- BC[ixt+1]/(1+r) + minCons
  }
  
  # maximum asset
  maxA[1] <- startA
  for (ixt in 2:(T+1)){
    maxA[ixt] <- maxA[ixt-1] * (1+r)
  }
  
  return(list(BC, maxA))
}