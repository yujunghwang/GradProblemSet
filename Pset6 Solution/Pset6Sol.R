# pset 6
# written by Yujung Hwang

rm(list=ls()) # clean memory
# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

library(pracma)
library(ggplot2)

# load functions
source("getMinAndMaxAssBCInc.R")
source("getGrid.R")
source("SolveDCValueFunctionBCInc.R")


interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption

# Grids
numPointsA <-1000
# method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps, or 10logsteps
gridMethod <- "logsteps"

# set values of structural economic parameters
T      <- 20
r      <- 0
beta   <- 0.98
gamma  <- 1 # log utility
delta  <- 1
  
borrowingAllowed <- 1
startA <- 0
ybar   <- 20

#----------------------------#
# Get Asset Grid differ by work status

Agrid <- array(rep(NA,(T+1)*numPointsA),c((T+1),numPointsA))

oout <- getMinAndMaxAssBCInc(borrowingAllowed,minInc=rep(0,T),maxInc=rep(ybar,T),startA)
borrowCon <- oout[[1]]
maxAss    <- oout[[2]]

for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
}


#--------------------#

# Graph using Value Function Iteration
solveUsingValueFunction <-1
solveUsingEulerEquation <-0
linearise  <- 1
oout_a     <- solveDCValueFunctionBCInc()
policyA1_a <- oout_a[[1]]
policyC_a  <- oout_a[[2]]
policyD_a  <- oout_a[[3]]
V_a        <- oout_a[[4]]
dU_a       <- oout_a[[5]]
optV_a     <- oout_a[[6]]
optdU_a    <- oout_a[[7]]

#--------------------------------------------------#
# graph only VI since Euler solution is not valid
#--------------------------------------------------#


# Graph policy function by time
poldat <- data.frame(c18_a=policyC_a[cbind(rep(18,numPointsA),c(1:numPointsA),policyD_a[18,])],a18=Agrid[18,],
                     c15_a=policyC_a[cbind(rep(15,numPointsA),c(1:numPointsA),policyD_a[15,])],a15=Agrid[15,],
                     c10_a=policyC_a[cbind(rep(10,numPointsA),c(1:numPointsA),policyD_a[10,])],a10=Agrid[10,])

# report kinks at t=18, t=15, t=10, t=1
ggplot(aes(x=a18,y=c18_a,color="t=18",linetype="t=18"),dat=poldat) + geom_line() + 
  geom_line(aes(x=a15,y=c15_a,color="t=15",linetype="t=15")) +
  geom_line(aes(x=a10,y=c10_a,color="t=10",linetype="t=10")) +
  scale_colour_manual(values=c("t=18"="black","t=15"="red","t=10"="blue"),name="") +
  scale_linetype_manual( values=c("t=18"=1,"t=15"=2,"t=10"=1) ,name="") +
  ylab("Consumption") + xlab("Asset") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 
ggsave('Consumption.png')


# graph nonmonotonicity of optdU
dUdat <- data.frame(dU18_a=optdU_a[cbind(rep(18,numPointsA),c(1:numPointsA))],a18=Agrid[18,],
                    dU15_a=optdU_a[cbind(rep(15,numPointsA),c(1:numPointsA))],a15=Agrid[15,],
                    dU10_a=optdU_a[cbind(rep(10,numPointsA),c(1:numPointsA))],a10=Agrid[10,])

ggplot(aes(x=a18,y=dU18_a,color="t=18",linetype="t=18"),dat=dUdat) + geom_line() + 
  geom_line(aes(x=a15,y=dU15_a,color="t=15",linetype="t=15")) +
  geom_line(aes(x=a10,y=dU10_a,color="t=10",linetype="t=10")) +
  scale_colour_manual(values=c("t=18"="black","t=15"="red","t=10"="blue"),name="") +
  scale_linetype_manual( values=c("t=18"=1,"t=15"=2,"t=10"=1) ,name="") +
  ylab("MU_c") + xlab("Asset") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30)) 
ggsave('MUc.png')


# conditional value function and value function
# find the primary kink point
prkink <- which(optV_a[20,]>V_a[20,,2])[1]
print(paste0("primary kink point (of asset) is ",Agrid[20,prkink]))

valdat <- data.frame(v20_a=optV_a[cbind(rep(20,numPointsA),c(1:numPointsA))],a20=Agrid[20,],
                     wv20_a=V_a[cbind(rep(20,numPointsA),c(1:numPointsA),rep(2,numPointsA))],
                     nv20_a=V_a[cbind(rep(20,numPointsA),c(1:numPointsA),rep(1,numPointsA))],
                     prkink=rep(Agrid[20,prkink],numPointsA))

ggplot(aes(x=a20,y=v20_a,color="Uncond.",linetype="Uncond."),dat=valdat,size=1) + geom_line() + 
  geom_line(aes(x=a20,y=wv20_a,color="Cond. Working",linetype="Cond. Working")) +
  geom_line(aes(x=a20,y=nv20_a,color="Cond. Not Working",linetype="Cond. Not Working")) +
  geom_vline(aes(xintercept=prkink),color="purple",size=1) +
  scale_colour_manual(values=c("Uncond."="black","Cond. Working"="red", "Cond. Not Working"="blue"),name="") +
  scale_linetype_manual( values=c("Uncond."=1,"Cond. Working"=2,"Cond. Not Working"=2) ,name="") +
  ylab("Value") + xlab("Asset") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=20),
        legend.position="bottom") 
ggsave('ValueFun.png')

