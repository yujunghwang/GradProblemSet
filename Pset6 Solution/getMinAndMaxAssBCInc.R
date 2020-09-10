getMinAndMaxAssBCInc <- function(borrowingAllowed, minIncome, maxIncome, startA){
  
  # A fuction that gives the minimum and maximum on the asset grid in each
  # year. The minimum is the natural borrowing constraint. The maximum is how
  # much one would have if one saved everything
  
  # The following variables are global variables
  # T r minCons
  
  # Initialize the outputs
  BC   <- rep(NA, T+1)
  maxA <- rep(NA, T+1)
  
  # Iteratively calculate the borrowing constraints, and maximum on asset grid
  
  # borrowing constraints
  BC[T+1] <-0
  for (ixt in T:1){
    BC[ixt] <- BC[ixt+1]/(1+r) -minIncome[ixt] + minCons
  }
  
  # if borrowing is not allowed, replace negative points in the borrowing
  # constraint with zero
  if (borrowingAllowed == 0 ){
    BC[BC<0] <- 0;
  }
  
  # maximum asset

  maxA[1] <- startA 
  for (ixt in 2:(T+1)){
    maxA[ixt] <- (maxA[ixt-1] + maxIncome[ixt-1]) * (1+r)
  }
  
  #If starting assets are 0 we will have maxA<0 in period 1 and BC>0.
  #Therefore BC>maxA (minA<maxA). Ensure that this is not the case by
  #replacing any element of maxA that is less than the corresponding element
  #of BC with that corresponding element plus 1
  for (ixt in 1:(T+1)){
    if (maxA[ixt] <= BC[ixt]){
      maxA[ixt] <- BC[ixt] + 1
    }  
  }
  
  
  return(list(BC, maxA))
}