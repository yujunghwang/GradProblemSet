#---------------------------#
# Pset 1 Solution
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
source("getMinAndMaxAss.R")
source("getGrid.R")
source("solveValueFunction.R")
source("solveEulerEquation.R")
source("getPolicy_analytical.R")
source("simNoUncer.R")

# choose interpolation method from "linear", "nearest", "spline", "cubic"
interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption

# set values of structural economic parameters
T      <- 40
r      <- 0.03
beta   <- 0.95
gamma  <- 1.5
startA <- 1

#----------------------------#
# Grids
numPointsA <-20
# method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps, or 10logsteps
gridMethod <- "logsteps"

#----------------------------#
# Get Asset Grid
# populate grid for assets using 'gridMethod'
oout <- getMinAndMaxAss(startA)
MinAss <- oout[[1]]
MaxAss <- oout[[2]]

Agrid <- matrix(rep(NA,(T+1)*numPointsA),ncol=numPointsA)
for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
}

#----------------------------#
# Question 2.
#----------------------------#
# get analytical policy functions if we have no uncertainty and borrowing is allowed
oout <- getPolicy_analytical()
policyA1_analytic <- oout[[1]]
policyC_analytic  <- oout[[2]]

# Graph (a) using value function
solveUsingValueFunction <-1
solveUsingEulerEquation <-0
linearise <- 0

oout_a <- solveValueFunction()
policyA1_a <- oout_a[[1]]
policyC_a  <- oout_a[[2]]
V_a        <- oout_a[[3]]

# Graph (b) using Euler equation without linearization
solveUsingValueFunction <-0
solveUsingEulerEquation <-1
linearise <- 0
oout_b <- solveEulerEquation()
policyA1_b <- oout_b[[1]]
policyC_b  <- oout_b[[2]]
V_b        <- oout_b[[3]]

# Graph (c) using Euler equation with linearization
solveUsingValueFunction <-0
solveUsingEulerEquation <-1
linearise <- 1
oout_c <- solveEulerEquation()
policyA1_c <- oout_c[[1]]
policyC_c  <- oout_c[[2]]
V_c        <- oout_c[[3]]


# plot
dat <- data.frame(asset=Agrid[20,],true=policyC_analytic[20,],sol_a=policyC_a[20,],sol_b=policyC_b[20,],sol_c=policyC_c[20,])
ggplot(dat,aes(x=asset,y=true,color="Analytical",linetype="Analytical"),size=1.5) + geom_line() + 
  geom_line(aes(y=sol_a,color="Graph (a)",linetype="Graph (a)"),size=1.5) +
  geom_line(aes(y=sol_b,color="Graph (b)",linetype="Graph (b)"),size=1.5) + 
  geom_line(aes(y=sol_c,color="Graph (c)",linetype="Graph (c)"),size=1.5) +
  scale_colour_manual(values=c("Analytical"="red", "Graph (a)" ="black",
                               "Graph (b)"="blue", "Graph (c)" ="forestgreen"),name="") + 
  scale_linetype_manual(values=c("Analytical"=1, "Graph (a)" =2,
                                 "Graph (b)"=3, "Graph (c)" =4),name="") +
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
# simulate consumption and asset path for 10 individuals using solution (c)
nsim <-10
simout <- simNoUncer(policyA1_c,V_c,startA)
cpath_c <- simout[[1]]
apath_c <- simout[[2]]
# truncate last a
apath_c <- apath_c[1:T]

# plot
ddat <- data.frame(time=1:T,cons=cpath_c,asset=apath_c)
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


