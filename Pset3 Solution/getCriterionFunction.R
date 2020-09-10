getCriterionFunction <- function(paramForOpt,WeightMatrix){

  diff                 <- findMomentVec(paramForOpt)
  criterionFunction <- t(diff)%*%WeightMatrix%*%diff

  return(criterionFunction)
}

