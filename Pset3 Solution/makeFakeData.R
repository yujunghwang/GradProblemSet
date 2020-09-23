#-----------------------#
# This file simulates fake data for estimation exercise
#-----------------------#

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
source("solveValueFunctionBCInc.R")
source("solveEulerEquationBCInc.R")
source("simWithUncerInc.R")
source("simNoUncerInc.R")

#------------------------#
# set parameters to generate data
#------------------------#

beta   <- 0.96
gamma  <- 1.9

#------------------------#

# choose interpolation method from "linear", "nearest", "spline", "cubic"
interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption

# set values of structural economic parameters
T      <- 10
r      <- 0.03
borrowingAllowed <- 1
isUncertainty <- 1
startA <- 0

#----------------------------#
# Grids
numPointsA <-3
# method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps, or 10logsteps
gridMethod <- "logsteps"

# where to truncate the normal distributions
numPointsY <- 3
#  points in grid for income 
normBnd <- 3
# ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma 
mu     <- 0
# increase income uncertainty
sigma <- 1
rho <- 0.4
Tretire <- 11 ### meaningless, ignore this.

#---------
# get income grid
oout  <- getIncomeGrid()
Ygrid <- oout[[1]]
incTransitionMrx <- oout[[2]]
minInc <- oout[[3]]
maxInc <- oout[[4]]

#----------------------------#
# Get Asset Grid
oout <- getMinAndMaxAssBCInc(borrowingAllowed,minInc,maxInc,startA)
borrowCon <- oout[[1]]
maxAss    <- oout[[2]]

Agrid <- matrix(rep(NA,(T+1)*numPointsA),ncol=numPointsA)
for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
}

# Graph (c) using Euler equation with linearization
solveUsingValueFunction <-0
solveUsingEulerEquation <-1
linearise  <- 1
oout_c     <- solveEulerEquationBCInc()
policyA1_c <- oout_c[[1]]
policyC_c  <- oout_c[[2]]
V_c        <- oout_c[[3]]

# simulate consumption and asset path for 100 individuals using solution (c)
numSims <-3000
# Draw random draws for starting income and for innovations
seed1 <- 10239402
seed2 <- 53927204
simout <- simWithUncerInc(policyA1_c,V_c,startA,seed1,seed2)
ypath_c <- simout[[1]]
cpath_c <- simout[[2]]
apath_c <- simout[[3]]
# truncate last a
apath_c <- apath_c[1:T,]

id  <- matrix(kron(seq(1,numSims,1),rep(1,T)),ncol=1)
age <- matrix(kron(rep(1,numSims),seq(1,T,1)),ncol=1)
dim(ypath_c) <- c(numSims*T,1)
dim(cpath_c) <- c(numSims*T,1)
dim(apath_c) <- c(numSims*T,1)

simData <- data.frame(id=id,age=age,inc=ypath_c,cons=cpath_c,asset=apath_c)
write.csv(simData,file="Pset3Data.csv")
