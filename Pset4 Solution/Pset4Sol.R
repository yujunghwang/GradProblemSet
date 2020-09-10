# Pset 4
rm(list=ls()) # clean memory
# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

# load libraries
#install.packages("pracma") # install packages if you haven't installed these yet.
#install.packages("sampleSelection")
library(pracma)
library(sampleSelection)

# read data
data <- read.csv("Pset4Data.csv",header=TRUE) ### data was already sorted by wave and id
data$age2 <- data$age^2
data$educ <- relevel(data$educ,ref="High School Graduate")

nobs <- dim(data)[1]
nsim <- max(data$id)
T    <- max(data$wave)

source("computeMoments.R")
source("computeMoments2.R")

oout <- computeMoments2(data)
Ymm <- oout[[1]]
Xmm <- oout[[2]]

# moment : duhat , duhat(duhat_1 + duhat + duhat+1), duhat^2 
GMMCriterionFunction <- function(param,wghtmat,Ymm,Xmm){
  
  # param 
  # sigma_zeta^2, sigma_zeta eta, sigma m^2
  sigmazeta2   <- param[1]
  sigmazetaeta <- param[2]
  sigmam2      <- param[3]
  
  diff <- c( Ymm[1] - sigmazetaeta*Xmm[1],
             Ymm[2] - sigmazeta2 - (sigmazetaeta^2)*Xmm[2],
             Ymm[3] - sigmazeta2 - (sigmazetaeta^2)*Xmm[3] - 2*sigmam2)
  
  dim(diff) <-c(3,1)
  
  val <- t(diff)%*%wghtmat%*%diff
  
  return(val)
}


lb <- c(0.00001,-10,0.00001)
ub <- c(10,10,10)

# Nonlinear Least Squares
oout <- optim(par=c(1,0,1),GMMCriterionFunction,wghtmat=eye(3),Ymm=Ymm,Xmm=Xmm,lower=lb,upper=ub,control=list(maxit=100),
      method="L-BFGS-B")

# get the estimates
paramest <- oout$par

# bootstrap standard errors

set.seed(20200217) # set seed
nb   <- 100 # number of resample
draw <- randi(nsim,nsim,nb)
nvar <- dim(data)[2]
parlen <- length(paramest)
bmoments <- matrix(rep(NA,parlen*nb),ncol=nb)

### below part can be parallelized using foreach()
for (k in 1:nb){
  
  bdata <- data ### use same data class as data
  colnames(bdata) <- colnames(data)
  for (l in 1:nsim){
    bdata[((l-1)*T+1):(l*T),1:nvar] <- data[data$id==draw[l,k],1:nvar] 
  }
  # sort bdata by wave and id
  bdata <- bdata[order(bdata$wave,bdata$id),] 
  
  # compute moments
  bout <- computeMoments2(data=bdata)
  bYmm <- bout[[1]]
  bXmm <- bout[[2]]
  rm(bout)
  
  # Nonlinear Least Squares
  bout <- optim(par=c(1,0,1),GMMCriterionFunction,wghtmat=eye(3),Ymm=bYmm,Xmm=bXmm,lower=lb,upper=ub,control=list(maxit=100),
                method="L-BFGS-B")
  
  # get the estimates
  bmoments[,k] <- bout$par
  rm(bdata,bout)
}

dbmoments <- bmoments -matrix(rep(paramest,nb),ncol=nb)
stderr <- sqrt(apply((dbmoments)^2,1,mean))

