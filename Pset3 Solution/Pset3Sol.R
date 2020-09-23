#---------------------------#
# Pset 3 Solution
# written by Yujung Hwang
# building on codes prepared by Cormac O'Dea and Monica Costa Dias
#--------------------------#
rm(list=ls()) # clean memory
# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

# load libraries
#install.packages("pracma") # install packages if you haven't installed these yet.
#install.packages("ggplot2")
library(pracma)
library(ggplot2)

# load functions
source("getMinAndMaxAssBCInc.R")
source("getGrid.R")
source("getIncomeGrid.R")
source("solveValueFunctionBCIncEst.R")
source("solveEulerEquationBCIncEst.R")
source("simWithUncerInc.R")
source("simNoUncerInc.R")
source("getCriterionFunction.R")
source("findMomentVec.R")

# read data
data <- read.csv("Pset3Data.csv",header=TRUE)

#---------------------#
# construct data moments which will be used for estimation
N <- max(data$id)
T <- max(data$age)
C <- matrix(data$cons, ncol=N)
A <- matrix(data$asset,ncol=N)
Y <- matrix(data$inc,ncol=N)

momentA  <- apply(A,1,mean)
betaA    <- rep(NA,T)
betaY    <- rep(NA,T)


for (i in 1:T){
 oout <- lm(C[i,]~A[i,]+Y[i,])
 betaA[i] <- coef(oout)[2]
 betaY[i] <- coef(oout)[3]
}

dataMoments     <- c(momentA[2:T],betaY[1:(T-1)],betaA[2:(T-1)])
NdataMoments   <- length(dataMoments)

#---------------------#
# Set other parameters as known
# solution method : set 1 to solve using Value function
solveUsingValueFunction <-0
# solution method : set 1 to solve using Euler equation
solveUsingEulerEquation <-1
# linearize 0 or 1
linearise <- 1

# Information for simulations
numSims <-N # How many individuals to simulate

# choose interpolation method from "linear", "nearest", "spline", "cubic"
interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption

# set values of structural economic parameters
r      <- 0.03
borrowingAllowed <- 1
isUncertainty <- 1
startA <- 0

#----------------------------#
# Grids
# choose dimensions, set matrices and select methods to construct grids
# remember if you choose finer grids, it takes longer to solve and estimate.

# The grid for assets
numPointsA <-3
# method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps, or 10logsteps
gridMethod <- "logsteps"

# grid for income shocks
numPointsY <- 3
#  points in grid for income 
normBnd <- 3
# ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma 
mu     <- 0
sigma <- 1
rho <- 0.4
Tretire <- 11 ### ignore this 

#---------
# get income grid
oout  <- getIncomeGrid()
Ygrid <- oout[[1]]
incTransitionMrx <- oout[[2]]
minInc <- oout[[3]]
maxInc <- oout[[4]]

#----------------------------#
# Get Asset Grid
# populate grid for assets using 'gridMethod'
oout <- getMinAndMaxAssBCInc(borrowingAllowed,minInc,maxInc,startA)
borrowCon <- oout[[1]]
maxAss    <- oout[[2]]

Agrid <- matrix(rep(NA,(T+1)*numPointsA),ncol=numPointsA)
for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
}

#---------------------#

### compute dataMomentsVCM by bootstrap
### use Diagonally Weighted Mimimum Distance estimator (DWMD) to avoid small sample bias
set.seed(1039281) # set seed
nb   <- 100 # number of resample
draw <- randi(N,N,nb)
nvar <- dim(data)[2]
bmoments <- matrix(rep(NA,NdataMoments*nb),ncol=nb)

for (k in 1:nb){
  
  bdata <- as.data.frame(matrix(rep(NA,N*T*nvar),ncol=nvar))
  colnames(bdata) <- colnames(data)
  for (l in 1:N){
    bdata[((l-1)*T+1):(l*T),1:nvar] <- data[data$id==draw[l,k],1:nvar]
  }
  
  bC <- matrix(bdata$cons, ncol=N)
  bA <- matrix(bdata$asset,ncol=N)
  bY <- matrix(bdata$inc,ncol=N)
  bbetaA    <- rep(NA,T)
  bbetaY    <- rep(NA,T)
  for (i in 1:T){
    oout <- lm(bC[i,]~bA[i,]+bY[i,])
    bbetaA[i] <- coef(oout)[2]
    bbetaY[i] <- coef(oout)[3]
  }
  
  # compute moments
  bmomentA  <- apply(bA,1,mean)
  bmoments[,k] <- c(bmomentA[2:T],bbetaY[1:(T-1)],bbetaA[2:(T-1)])
  rm(bdata,bmomentA,bC,bA,bY,bbetaA,bbetaY)
}

dbmoments <- bmoments-matrix(rep(apply(bmoments,1,mean),nb),ncol=nb)

varMat <- (dbmoments%*%t(dbmoments))/(nb-1)
weightMat <- eye(NdataMoments)
diag(weightMat) <- diag(varMat)
weightMat <- solve(weightMat)

#---------------------#
# Estimate Parameters

# (b)
# report time
ptm <- proc.time()
getCriterionFunction(c(0.96,1.9),WeightMatrix=weightMat)
proc.time() - ptm


# (c)
lb <- c(0.90, 1)
ub <- c(0.99,2.5)

# more generous bounds
#lb <- c(0.4,1.1)
#ub <- c(0.99,9)

# set environment of getCriterionFunction to GlobalEnv
oout <- optim(par=c(0.96,1.9),fn=getCriterionFunction,WeightMatrix=weightMat,lower=lb,upper=ub,control=list(maxit=100),
                        method="L-BFGS-B")
# print out estimate
print(oout)

save.image("Pset3Sol200917.RData")

# (d) Report the standard error
dMoment <- jacobian(findMomentVec,oout$par) 
covar <- (1+1/1)*solve(t(dMoment)%*%weightMat%*%dMoment)%*%t(dMoment)%*%weightMat%*%varMat%*%weightMat%*%dMoment%*%solve(t(dMoment)%*%weightMat%*%dMoment)
standerr <- sqrt(diag(covar))
print(standerr)


# (e) sensitivity measure
sensitivity <- -solve(t(dMoment)%*%weightMat%*%dMoment)%*%t(dMoment)%*%weightMat


# visually graph sensitivity measures
sensitivity <- as.data.frame(t(sensitivity))
colnames(sensitivity) <- c("beta","gamma")
sensitivity$age <- c(c(2:T), c(1:(T-1)), c(2:(T-1)))
sensitivity$grp <- as.factor(c(rep("Mean Asset",(T-1)), rep("coef Y",(T-1)), rep("coef A", (T-2))))

png("BetaSensitivity.png")
colvec <- c("black","blue","red")
dotchart(sensitivity$beta, labels = sensitivity$age,
         groups = sensitivity$grp, gcolor = colvec,
         cex = 0.6,  pch = 19, xlab = "measure", main="Beta",
         color=colvec[sensitivity$grp], lcolor = "gray")
dev.off()


png("GammaSensitivity.png")
dotchart(sensitivity$gamma, labels = sensitivity$age,
         groups = sensitivity$grp, gcolor = colvec,
         cex = 0.6,  pch = 19, xlab = "measure", main="Gamma",
         color=colvec[sensitivity$grp], lcolor="gray")
dev.off()



# (f) model fit
simMoments <- findMomentVec(oout$par) + dataMoments
simdata <- data.frame(time=1:T,simC=simMoments[1:T],simA=simMoments[(T+1):(2*T)],
                      dataC=dataMoments[1:T],dataA=dataMoments[(T+1):(2*T)])

ggplot(simdata,aes(x=time,y=simC,color="Sim C",linetype="Sim C"),size=1.5) + 
  geom_line() + geom_line(aes(y=dataC,color="Data C",linetype="Data C"),size=1.5) +
  scale_colour_manual(values=c("Sim C"="red", "Data C"="black"),name="") + 
  scale_linetype_manual(values=c("Sim C"=1, "Data C"=2),name="") +
  ylab("Mean Consumption") + xlab("Time") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 
ggsave('modelFitCons.png')


ggplot(simdata,aes(x=time,y=simA,color="Sim A",linetype="Sim A"),size=1.5) + 
  geom_line() + geom_line(aes(y=dataA,color="Data A",linetype="Data A"),size=1.5) +
  scale_colour_manual(values=c("Sim A"="red", "Data A"="black"),name="") + 
  scale_linetype_manual(values=c("Sim A"=1, "Data A"=2),name="") +
  ylab("Mean Asset") + xlab("Time") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 
ggsave('modelFitAssets.png')