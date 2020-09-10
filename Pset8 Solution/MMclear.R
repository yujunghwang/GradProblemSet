# function
# input : Pareto weight matrix
# output : excess demand vector
MMclear <- function(lambda){
  # lambda is a Pareto weight vector 
  
  # compute demand of type i men
  supply <- supplyMen(lambda)
  
  # compute supply of type i men
  demand <- demandMen(lambda)
  
  # take difference and compute excess demand
  excessDemand <- demand - supply
  
  return(excessDemand)
}