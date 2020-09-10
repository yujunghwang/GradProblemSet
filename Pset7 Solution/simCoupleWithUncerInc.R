simCoupleWithUncerInc <- function(marA1,marC1,marC2,marD,marTh1,divA1,divC,startA,startTheta,seed1,seed2){
  # This function takes the policy functions and value functions 
  # along with starting assets and returns simulated paths of income, consumption
  # assets and value
  
  #-------------------
  # initialize arrays that will hold the paths of income consumption, value and assets
  
  # arguments for output
  y1 <- matrix(rep(NA,T*numSims),ncol=numSims)
  y2 <- matrix(rep(NA,T*numSims),ncol=numSims)
  c1 <- matrix(rep(NA,T*numSims),ncol=numSims)
  c2 <- matrix(rep(NA,T*numSims),ncol=numSims)
  d  <- matrix(rep(NA,T*numSims),ncol=numSims)
  th <- matrix(rep(NA,(T+1)*numSims),ncol=numSims)
  a  <- matrix(rep(NA,(T+1)*numSims),ncol=numSims) ## couple's asset
  
  da1 <- matrix(rep(NA,(T+1)*numSims),ncol=numSims) ## divorcee's asset
  da2 <- matrix(rep(NA,(T+1)*numSims),ncol=numSims)  
  
  # Other arrays that will be used below
  ly1         <- matrix(rep(NA,T*numSims),ncol=numSims)
  ly2         <- matrix(rep(NA,T*numSims),ncol=numSims)

  #---------------------------#
  # Obtain time series of incomes for our simulated individuals
  #---------------------------#
  sig_inc <- sigma/ ((1-rho^2)^0.5)
  
  e1    <- getNormalDraws(0,  sigma,  T,numSims,seed1) # normally distributed random draws
  logy1 <- getNormalDraws(mu, sig_inc,1,numSims,seed2) # a random draw for the initial income
  
  e2    <- getNormalDraws(0,  sigma,  T,numSims,(seed1+10)) # normally distributed random draws
  logy2 <- getNormalDraws(mu, sig_inc,1,numSims,(seed2+10)) # a random draw for the initial income
  
  # Get all the incomes, recursively
  for (s in 1:numSims){ # loop through individuals
    
    ly1[1,s] <- min(max(logy1[1,s],-normBnd*sig_inc),normBnd*sig_inc)
    y1[ 1,s] <- exp(ly1[1,s])
    
    ly2[1,s] <- min(max(logy2[1,s],-normBnd*sig_inc),normBnd*sig_inc)
    y2[ 1,s] <- exp(ly2[1,s])    
    
    for (t in 1:T){
      
      if (t!=T){ # Get next year's income
        ly1[t+1,s] <- (1-rho)*mu + rho*ly1[t,s] + e1[(t+1),s]
        ly1[t+1,s] <- min( max(ly1[t+1,s],-normBnd*sig_inc), normBnd*sig_inc)
        y1[t+1,s]  <- exp(ly1[t+1,s])
        
        ly2[t+1,s] <- (1-rho)*mu + rho*ly2[t,s] + e2[(t+1),s]
        ly2[t+1,s] <- min( max(ly2[t+1,s],-normBnd*sig_inc), normBnd*sig_inc)
        y2[t+1,s]  <- exp(ly2[t+1,s])
      } # if (t != T)
      
      if (t>=Tretire){
        y1[t,s] <- 0
        y2[t,s] <- 0
      } # if (t>=Tretire)
      
    } # t
  } # s
  
  #-----------------------------#
  # Obtain consumption, asset, and value profiles
  #-----------------------------#
  
  # initialize
  a[1,]  <- startA
  th[1,] <- startTheta 
  
  
  for (t in 1:T){
    
    # interpolate policy functions
    grid  <- list( Agrid[t,],Ygrid[t,],Ygrid[t,],Thetagrid)
    imarD <- ipol(marD[ t,,,,],grid=grid, method='multilin')
    imarA1<- ipol(marA1[t,,,,],grid=grid, method='multilin')
    imarC1<- ipol(marC1[t,,,,],grid=grid, method='multilin')
    imarC2<- ipol(marC2[t,,,,],grid=grid, method='multilin')
    imarTh<- ipol(marTh1[t,,,,],grid=grid,method='multilin')
    
    sgrid  <- list(Agrid[t,],Ygrid[t,])
    idivA1 <- ipol(divA1[t,,],grid=sgrid,method='multilin')
    idivC  <- ipol( divC[t,,],grid=sgrid,method='multilin') 

    
    # simulate agents within period t
    for (s in 1:numSims){
     
      # simulate divorce decision
      if (t==1){
        divprob <- imarD( c(a[t,s],y1[t,s],y2[t,s],th[t,s]) )
        set.seed((seed1+10*s+t))
        d[t,s]  <- which(c(1-divprob,1)>as.numeric(rand(1)) )[1] -1 
      } else {
        if (d[(t-1),s]==0){
         # simulate divorce decision only if not divorced last period 
         divprob <- imarD( c(a[t,s],y1[t,s],y2[t,s],th[t,s]) )
         set.seed((seed1+10*s+t))
         d[t,s]  <- which(c(1-divprob,1)>as.numeric(rand(1)) )[1] -1  
        } else{
         d[t,s] <- 1 
        }
      }
      
      #----------------#
      if (d[t,s]==1){
        # divorcee's simulation 
        
        if (t>1){
          # newly divorced
          if (d[(t-1),s]==0){
            da1[(t+1),s] <- idivA1( c(a[t,s]/2, y1[t,s] )) # 50% asset division upon divorce
            da2[(t+1),s] <- idivA1( c(a[t,s]/2, y2[t,s] ))
            c1[t,s]      <- idivC(  c(a[t,s]/2, y1[t,s] ))
            c2[t,s]      <- idivC(  c(a[t,s]/2, y2[t,s] ))
            th[(t+1),s]  <- th[t,s]
            
          } else{
            # old divorcee
            da1[(t+1),s] <- idivA1( c(da1[t,s], y1[t,s] )) # 50% asset division upon divorce
            da2[(t+1),s] <- idivA1( c(da2[t,s], y2[t,s] ))
            c1[t,s]      <- idivC(  c(da1[t,s], y1[t,s] ))
            c2[t,s]      <- idivC(  c(da2[t,s], y2[t,s] ))
            th[(t+1),s]  <- th[t,s]
          }

          
        } else{
          # newly divorced
          da1[(t+1),s] <- idivA1( c(a[t,s]/2, y1[t,s] )) # 50% asset division upon divorce
          da2[(t+1),s] <- idivA1( c(a[t,s]/2, y2[t,s] ))
          c1[t,s]      <- idivC(  c(a[t,s]/2, y1[t,s] ))
          c2[t,s]      <- idivC(  c(a[t,s]/2, y2[t,s] ))
          th[(t+1),s]  <- th[t,s]
        }
       
        
      } else{
        # married couple's simulation
        a[(t+1),s] <- imarA1(c(a[t,s],y1[t,s],y2[t,s],th[t,s])) 
        c1[t,s]    <- imarC1(c(a[t,s],y1[t,s],y2[t,s],th[t,s])) 
        c2[t,s]    <- imarC2(c(a[t,s],y1[t,s],y2[t,s],th[t,s])) 
        th[(t+1),s]<- imarTh(c(a[t,s],y1[t,s],y2[t,s],th[t,s])) 
      }
      
      # due to interpolatin some consumption values may fall below 0. correct it.
      if (c1[t,s]<0 & d[t,s]==0){
        a[(t+1),s] <- a[(t+1),s] + (1+r)*(c1[t,s] - minCons)
        c1[t,s]<- minCons                  
      } 
      
      if (c1[t,s]<0 & d[t,s]==1){
        da1[(t+1),s]<- da1[(t+1),s] + (1+r)*(c1[t,s] - minCons)
        c1[t,s]     <- minCons
      }
      
      if (c2[t,s]<0 & d[t,s]==0){
        a[(t+1),s] <- a[(t+1),s] + (1+r)*(c2[t,s] - minCons)
        c2[t,s]    <- minCons  
      }
      if (c2[t,s]<0 & d[t,s]==1){
        da2[(t+1),s]<- da2[(t+1),s] + (1+r)*(c2[t,s] - minCons)
        c2[t,s]     <- minCons
      }
      
    } #s
  } #t
  
  
  return(list(y1,y2,c1,c2,a,da1,da2,d,th))
}