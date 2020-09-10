computeMoments2 <- function(data){
  
  workeq <- glm(work ~ age + age2 + factor(educ) + noveliv, family = binomial(link = "probit"), data = data)
  # (a)
  #data$lambdahat <- invMillsRatio(workeq)$IMR1 ### produces same result.
  
  alphahat <- -cbind(rep(1,nobs),data$age,data$age^2,(data$educ=="College Graduate"),
                     (data$educ=="Postgraduate"),(data$educ=="Some College"),data$noveliv)%*%workeq$coefficients
  lambdahat <- dnorm(alphahat)/(1-pnorm(alphahat))
  
  age      <- data$age
  dim(age) <- c(nsim,T)
  d1age2   <- age^2 - (age-1)^2 
  dim(d1age2) <-c(nsim*T,1)
  
  # (b)
  #### work with differenced equation
  # lnwage_t - lnwage_{t-2}
  lnwage <- data$lninc
  dim(lnwage)<-c(nsim,T)
  lnwagep1 <- cbind(rep(NA,nsim),lnwage[,1:(T-1)]) ### data was already sorted by wave and id
  lnwagep2 <- cbind(rep(NA,nsim),rep(NA,nsim),lnwage[,1:(T-2)])
  
  d2lnwage <- lnwage - lnwagep2
  d1lnwage <- lnwage - lnwagep1
  
  dim(d1lnwage) <- c(nsim*T,1)
  
  reg1  <- lm(d1lnwage ~ d1age2 + lambdahat)
  # compute residual and find uhat
  duhat <- d1lnwage - reg1$coefficients[1] - reg1$coefficients[2]*d1age2
  
  dim(duhat) <- c(nsim,T)
  duhatp1 <- cbind(rep(NA,nsim),duhat[,1:(T-1)])
  duhatl1 <- cbind(duhat[,2:T],rep(NA,nsim))
  
  L      <- as.numeric(data$work)
  dim(L) <- c(nsim,T) # L_{it}
  
  Lp1 <- cbind(rep(0,nsim),L[,1:(T-1)]) #L_{it-1}
  Lp2 <- cbind(rep(0,nsim),rep(0,nsim),L[,1:(T-2)]) #L_{it-2}
  Ll1 <- cbind(L[,2:T],rep(0,nsim)) #L_{it+1}
  dim(L)   <-c(nsim*T,1)
  dim(Lp1) <-c(nsim*T,1)
  dim(Lp2) <-c(nsim*T,1)
  dim(Ll1) <-c(nsim*T,1)
  
  # indicator L_{it}=1, L_{it-1}=1
  iind1 <- rep(0,nsim*T)
  iind1[L==1 & Lp1==1] <- 1
  
  # indicator L_{it-2}=1,L_{it-1}=1,L_{it}=1,L_{it+1}=1
  iind2 <- matrix(rep(0,nsim*T),ncol=T)
  iind2[Lp2==1 & Lp1==1 & L==1 & Ll1 ==1] <- 1
  
  # transform the data size
  dim(duhat)    <-c(nsim*T,1)
  dim(duhatp1)  <-c(nsim*T,1)
  dim(duhatl1)  <-c(nsim*T,1)
  dim(iind1)    <-c(nsim*T,1)
  dim(iind2)    <-c(nsim*T,1)
  
  
  # generate data moments
  Ymm <- c(mean(duhat[iind1==1]),
           mean(duhat[iind2==1]*(duhatp1[iind2==1] + duhat[iind2==1] + duhatl1[iind2==1])),
           mean(duhat[iind1==1]^2))

  Xmm <- c(mean(lambdahat[iind1==1]),
           mean(lambdahat[iind2==1]*alphahat[iind2==1]),
           mean(lambdahat[iind1==1]*alphahat[iind1==1]) )
 
  return(list(Ymm,Xmm)) 
}
