#---------------------------#
# Pset 2 Solution
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
source("solveValueFunctionBCInc.R")
source("solveEulerEquationBCInc.R")
source("simWithUncerInc.R")
source("simNoUncerInc.R")

# choose interpolation method from "linear", "nearest", "spline", "cubic"
interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption

# set values of structural economic parameters
T      <- 40
r      <- 0.03
beta   <- 0.95
gamma  <- 1.5
borrowingAllowed <- 1
isUncertainty <- 1
startA <- 0

#----------------------------#
# Grids
numPointsA <-20
# method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps, or 10logsteps
gridMethod <- "logsteps"

# where to truncate the normal distributions
numPointsY <- 5
#  points in grid for income (should be 2 if hard-coded)
normBnd <- 3
# ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma 
mu     <- 0
sigma  <- 0.25
rho    <- 0.75
Tretire <- 41

#---------
# get income grid
oout  <- getIncomeGrid()
Ygrid <- oout[[1]]
incTransitionMrx <- oout[[2]]
minInc <- oout[[3]]
maxInc <- oout[[4]]

# print income 
options(digits=3)
print('Income Grid')
print(Ygrid)
print('Income Transition Matrix')
print(incTransitionMrx)

#----------------------------#
# Get Asset Grid
oout <- getMinAndMaxAssBCInc(borrowingAllowed,minInc,maxInc,startA)
borrowCon <- oout[[1]]
maxAss    <- oout[[2]]

Agrid <- matrix(rep(NA,(T+1)*numPointsA),ncol=numPointsA)
for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
}


#----------------------------#
# Question 2.
#----------------------------#

# Graph (a) using value function
solveUsingValueFunction <-1
solveUsingEulerEquation <-0
linearise <- 0

oout_a     <- solveValueFunctionBCInc()
policyA1_a <- oout_a[[1]]
policyC_a  <- oout_a[[2]]
V_a        <- oout_a[[3]]

# Graph (b) using Euler equation without linearization
solveUsingValueFunction <-0
solveUsingEulerEquation <-1
linearise  <- 0
oout_b     <- solveEulerEquationBCInc()
policyA1_b <- oout_b[[1]]
policyC_b  <- oout_b[[2]]
V_b        <- oout_b[[3]]

# Graph (c) using Euler equation with linearization
solveUsingValueFunction <-0
solveUsingEulerEquation <-1
linearise  <- 1
oout_c     <- solveEulerEquationBCInc()
policyA1_c <- oout_c[[1]]
policyC_c  <- oout_c[[2]]
V_c        <- oout_c[[3]]

# plot
dat <- data.frame(asset=Agrid[20,],
                  sol_a_h=policyC_a[20,,5],
                  sol_a_l=policyC_a[20,,1],
                  sol_b_h=policyC_b[20,,5],
                  sol_b_l=policyC_b[20,,1],
                  sol_c_h=policyC_c[20,,5],
                  sol_c_l=policyC_c[20,,1])


ggplot(dat,aes(x=asset,y=sol_a_h,color="Graph (a) highest inc",linetype="Graph (a) highest inc")) + geom_line(size=1.5) + 
  geom_line(aes(y=sol_a_l,color="Graph (a) lowest inc",linetype="Graph (a) lowest inc"),size=1.5) + 
  geom_line(aes(y=sol_b_h,color="Graph (b) highest inc",linetype="Graph (b) highest inc"),size=1.5) + 
  geom_line(aes(y=sol_b_l,color="Graph (b) lowest inc",linetype="Graph (b) lowest inc"),size=1.5) + 
  geom_line(aes(y=sol_c_h,color="Graph (c) highest inc",linetype="Graph (c) highest inc"),size=1.5) +
  geom_line(aes(y=sol_c_l,color="Graph (c) lowest inc",linetype="Graph (c) lowest inc"),size=1.5) +
  scale_colour_manual(values=c("Graph (a) highest inc" ="black",
                               "Graph (a) lowest inc"="black",
                               "Graph (b) highest inc"="blue",
                               "Graph (b) lowest inc"="blue",
                               "Graph (c) highest inc" ="forestgreen",
                               "Graph (c) lowest inc" ="forestgreen"),name="") + 
  scale_linetype_manual(values=c("Graph (a) highest inc" =1,
                                 "Graph (a) lowest inc"=4,
                                 "Graph (b) highest inc"=1, 
                                 "Graph (b) lowest inc"=4,
                                 "Graph (c) highest inc" =1,
                                 "Graph (c) lowest inc"=4),name="") +
  ylab("Consumption") + xlab("Asset") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 

ggsave('Question2.png')

#----------------------------#
# Question 3.
#----------------------------#
# simulate consumption and asset path for 100 individuals using solution (c)
numSims <-100
# Draw random draws for starting income and for innovations
seed1 <- 10239402
seed2 <- 53927204
simout <- simWithUncerInc(policyA1_c,V_c,startA,seed1,seed2)
cpath_c <- simout[[2]]
apath_c <- simout[[3]]
# truncate last a
apath_c <- apath_c[1:T,]


# compute consumption asset moments by age
# plot
ddat <- data.frame(time=1:T,cons=apply(cpath_c,1,mean),asset=apply(apath_c,1,mean))

ggplot(ddat,aes(x=time,y=cons,color="Consumption",linetype="Consumption"),size=1.5) + 
  geom_line() + geom_line(aes(y=asset,color="Asset",linetype="Asset"),size=1.5) +
  scale_colour_manual(values=c("Consumption"="red", "Asset"="black"),name="") + 
  scale_linetype_manual(values=c("Consumption"=1, "Asset"=2),name="") +
  ylab("Consumption/Asset") + xlab("Time") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 
ggsave('Question3.png')
