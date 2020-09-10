# full commitment
solveFullCoupleValueFunctionBCInc <- function(){
  # this function solves couple's LC problem
  
  # matrices to hold the policy, value and marginal utility functions
  marV  <- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )  ## married couple's value function
  marV1 <- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )  ## spouse member 1's value within marriage
  marV2 <- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )  ## spouse member 2's value within marriage
  marA1 <- rep(NA,  T*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marC1 <- rep(NA,  T*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marC2 <- rep(NA,  T*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marD  <- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marEV <- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marEV1<- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marEV2<- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta )
  marTh1<- rep(NA, (T+1)*numPointsA*numPointsY*numPointsY*numPointsTheta ) ## tomorrow's new bargaining weight
  
  # dim : time, asset, spouse 1 income, spouse 2 income, theta
  dim(marV)  <- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marV1) <- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marV2) <- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marA1) <- c( T,numPointsA,numPointsY,numPointsY,numPointsTheta )
  dim(marC1) <- c( T,numPointsA,numPointsY,numPointsY,numPointsTheta )
  dim(marC2) <- c( T,numPointsA,numPointsY,numPointsY,numPointsTheta )
  dim(marD)  <- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marEV) <- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marEV1)<- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marEV2)<- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  dim(marTh1)<- c((T+1),numPointsA,numPointsY,numPointsY,numPointsTheta)
  
  #--------------------------#
  # final period
  marEV[(T+1),,,,]  <-0
  marEV1[(T+1),,,,] <-0
  marEV2[(T+1),,,,] <-0
  marV[ (T+1),,,,]  <-0
  marV1[(T+1),,,,]  <-0
  marV2[(T+1),,,,]  <-0
  marD[ (T+1),,,,]  <-0
  marTh1[(T+1),,,,] <-0
  
  #--------------------------#
  # backward induction
  
  for (ixt in T : 1){ # Loop from time T-1 to 1
    
    Agrid1 <- Agrid[ixt+1,] # The grid on assets tmrw
    
    for (ixA in 1:numPointsA){ # points on asset grid
      
      # Step 1. solve problem at grid points in assets and income
      #-------------
      for (ixY1 in 1:numPointsY){ # points on income grid (spouse 1)
        for (ixY2 in 1:numPointsY){ # points on income grid (spouse 2)
          for (ixTh in 1:numPointsTheta){ # points on theta grid
            
            # value of income and information for optimization
            A  <- Agrid[ixt,ixA]  # assets today
            Y1 <- Ygrid[ixt,ixY1] # income today for spouse 1
            Y2 <- Ygrid[ixt,ixY2] # income today for spouse 2
            
            lbA1 <- Agrid[ixt+1,1] # lower bound : assets tmrw
            ubA1 <- (A + Y1 + Y2 - minCons)*(1+r) # upper bound : assets tmrw
            
            marEVtm  <- marEV[ (ixt+1),,ixY1,ixY2,ixTh] # relevant section of EV matrix (in assets tmrw)
            marEV1tm <- marEV1[(ixt+1),,ixY1,ixY2,ixTh] # income and theta is in today's values (theta before updating by participation constraint)
            marEV2tm <- marEV2[(ixt+1),,ixY1,ixY2,ixTh]          
            
            oout <- solveCoupleIntraAlloc(A=A,Y1=Y1,Y2=Y2,theta=Thetagrid[ixTh],
                                          Agrid1=Agrid1,marEVtm=marEVtm,marEV1tm=marEV1tm,
                                          marEV2tm=marEV2tm,lbA1=lbA1,ubA1=ubA1)
            marA1[ixt,ixA,ixY1,ixY2,ixTh] <- oout[[1]]
            marC1[ixt,ixA,ixY1,ixY2,ixTh] <- oout[[2]]
            marC2[ixt,ixA,ixY1,ixY2,ixTh] <- oout[[3]]
            marV[ ixt,ixA,ixY1,ixY2,ixTh] <- oout[[4]]
            marV1[ixt,ixA,ixY1,ixY2,ixTh] <- oout[[5]]
            marV2[ixt,ixA,ixY1,ixY2,ixTh] <- oout[[6]]
            rm(oout)
            
            # compute outside option of 1 and 2
            outside1 <- interp1(Agrid[ixt,],divV[ixt,,ixY1],(Agrid[ixt,ixA]/2),method=interpMethod) # upon divorce, receive half asset
            outside2 <- interp1(Agrid[ixt,],divV[ixt,,ixY2],(Agrid[ixt,ixA]/2),method=interpMethod)
            
            # check participation constriants
            if ((marV1[ixt,ixA,ixY1,ixY2,ixTh]<outside1) & (marV2[ixt,ixA,ixY1,ixY2,ixTh]<outside2)){
              marD[ ixt,ixA,ixY1,ixY2,ixTh]  <-1
              marV[ ixt,ixA,ixY1,ixY2,ixTh]  <-Thetagrid[ixTh]*outside1 + (1-Thetagrid[ixTh])*outside2 ### sum of divorcee's values
              marTh1[ ixt,ixA,ixY1,ixY2,ixTh]<-Thetagrid[ixTh]
              
            } else {
              marD[ixt,ixA,ixY1,ixY2,ixTh]   <-0
              marTh1[ixt,ixA,ixY1,ixY2,ixTh] <-Thetagrid[ixTh] ## no change in bargaining weight
            }
          } # ixTh
        } # ixY2
      } # ixY1
      
      
      # Step 2 : update expected value at the beginning of period, before income transition.
      for (ixY1 in 1:numPointsY){
        for (ixY2 in 1:numPointsY){
          for (ixTh in 1:numPointsTheta){ # points on theta grid
            
            marEV[ ixt,ixA,ixY1,ixY2,ixTh] <- 0
            marEV1[ixt,ixA,ixY1,ixY2,ixTh] <- 0    
            marEV2[ixt,ixA,ixY1,ixY2,ixTh] <- 0
            
            for (ixY11 in 1:numPointsY){
              for (ixY21 in 1:numPointsY){
                marEV[ ixt,ixA,ixY1,ixY2,ixTh] <-  marEV[ ixt,ixA,ixY1,ixY2,ixTh]+
                  incTransitionMrx[ixY1,ixY11]*incTransitionMrx[ixY2,ixY21]*marV[ ixt,ixA,ixY11,ixY21,ixTh]
                marEV1[ixt,ixA,ixY1,ixY2,ixTh] <-  marEV1[ixt,ixA,ixY1,ixY2,ixTh]+ 
                  incTransitionMrx[ixY1,ixY11]*incTransitionMrx[ixY2,ixY21]*marV1[ixt,ixA,ixY11,ixY21,ixTh]
                marEV2[ixt,ixA,ixY1,ixY2,ixTh] <-  marEV2[ixt,ixA,ixY1,ixY2,ixTh]+
                  incTransitionMrx[ixY1,ixY11]*incTransitionMrx[ixY2,ixY21]*marV2[ixt,ixA,ixY11,ixY21,ixTh]
                
              }
            }
            
          }
        }
      }
      
      
    } # ixA
  } # ixt
  
  return(list(marA1,marC1,marC2,marD,marV,marV1,marV2,marEV,marEV1,marEV2,marTh1))  
}