#---------------------------#
# Pset 7 Solution
# written by Yujung Hwang
#--------------------------#

rm(list=ls()) # clean memory
# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

# load libraries
#install.packages("pracma") # install packages if you haven't installed these yet.
#install.packages("ggplot2")
#install.packages("chebpol")
library(pracma)
library(ggplot2)
library(chebpol)

# load functions
source("getMinAndMaxAssBCInc.R")
source("getGrid.R")
source("getIncomeGrid.R")
source("solveValueFunctionBCInc.R")
source("solveLimitedCoupleValueFunctionBCInc.R")
source("solveFullCoupleValueFunctionBCInc.R")
source("simWithUncerInc.R")
source("simCoupleWithUncerInc.R")

# set parameters
source("setParams.R")

# get income grid
oout  <- getIncomeGrid()
Ygrid <- oout[[1]]
incTransitionMrx <- oout[[2]]
minInc <- oout[[3]]
maxInc <- oout[[4]]

#### Get Asset Grid ####
#### construct an asset grid which can be used for both couple and single problem
#### upon divorce, the single is endowed with half of couple's asset by assumption
#### therefore, couple's maximum income affects single's asset grid range.
oout <- getMinAndMaxAssBCInc(borrowingAllowed,(minInc),(2*maxInc),startA)
borrowCon <- oout[[1]]
maxAss    <- oout[[2]]
rm(oout)

Agrid <- matrix(rep(NA,(T+1)*numPointsA),ncol=numPointsA)
for (ixt in 1:(T+1)){
  Agrid[ixt,] <- getGrid(borrowCon[ixt], maxAss[ixt], numPointsA, gridMethod)
}

#### Get theta (Pareto weight) grid ####
#### theta = spouse 1's bargaining weight
Thetagrid <- seq(from=0.1,to=0.9,by=0.2)
numPointsTheta <- length(Thetagrid)

##### divorcee's value function and policy functions (solution of pset 2) ######
solveUsingValueFunction <-1
solveUsingEulerEquation <-0
linearise <- 0
oout  <- solveValueFunctionBCInc()
divA1 <- oout[[1]]
divC  <- oout[[2]]
divV  <- oout[[3]]
rm(oout)

##### solve limited commitment couple's problem #####
ptm <- proc.time()
oout  <- solveLimitedCoupleValueFunctionBCInc()
proc.time() - ptm
LmarA1 <- oout[[1]]
LmarC1 <- oout[[2]]
LmarC2 <- oout[[3]]
LmarD  <- oout[[4]]
LmarV  <- oout[[5]]
LmarV1 <- oout[[6]]
LmarV2 <- oout[[7]]
LmarEV <- oout[[8]]
LmarEV1<- oout[[9]]
LmarEV2<- oout[[10]]
LmarTh1<- oout[[11]]
rm(oout)

##### solve full commitment couple's problem #####
ptm <- proc.time()
oout  <- solveFullCoupleValueFunctionBCInc()
proc.time() - ptm
FmarA1 <- oout[[1]]
FmarC1 <- oout[[2]]
FmarC2 <- oout[[3]]
FmarD  <- oout[[4]]
FmarV  <- oout[[5]]
FmarV1 <- oout[[6]]
FmarV2 <- oout[[7]]
FmarEV <- oout[[8]]
FmarEV1<- oout[[9]]
FmarEV2<- oout[[10]]
FmarTh1<- oout[[11]]
rm(oout)


#### simulated married couple's LC ####
numSims <- 1000
# Draw random draws for starting income and for innovations
seed1   <- 20200302
seed2   <- 20200303

startTheta <- 0.3

save.image(file="pset7Sol.RData")
  
simout  <- simCoupleWithUncerInc(marA1=LmarA1,marC1=LmarC1,marC2=LmarC2,marD=LmarD,marTh1=LmarTh1,
                                 divA1=divA1,divC=divC,startA=startA,startTheta=startTheta,seed1=seed1,seed2=seed2)
lcy1 <- simout[[1]]
lcy2 <- simout[[2]]
lcc1 <- simout[[3]]
lcc2 <- simout[[4]]
lca  <- simout[[5]]
lca1 <- simout[[6]]
lca2 <- simout[[7]]
lcd  <- simout[[8]]
lcth <- simout[[9]]

simout  <- simCoupleWithUncerInc(marA1=FmarA1,marC1=FmarC1,marC2=FmarC2,marD=FmarD,marTh1=FmarTh1,
                                 divA1=divA1,divC=divC,startA=startA,startTheta=startTheta,seed1=seed1,seed2=seed2)
fcy1 <- simout[[1]]
fcy2 <- simout[[2]]
fcc1 <- simout[[3]]
fcc2 <- simout[[4]]
fca  <- simout[[5]]
fca1 <- simout[[6]]
fca2 <- simout[[7]]
fcd  <- simout[[8]]
fcth <- simout[[9]]

#### report the Var(log(c1)), Var(log(c2)) under LC and FC
# under limited commitment
var(log(lcc1[1:(numSims*T)]))
var(log(lcc2[1:(numSims*T)]))

# under full commitment
var(log(fcc1[1:(numSims*T)]))
var(log(fcc2[1:(numSims*T)]))

# graph average log consumption over time
ddat <- data.frame(time=1:T,lc1=apply(log(lcc1),1,mean),  lc2=apply(log(lcc2),1,mean), 
                            fc1=apply(log(fcc1),1,mean),  fc2=apply(log(fcc2),1,mean) )

ggplot(ddat,aes(x=time,y=lc1,color="LC Spouse1",linetype="LC Spouse1")) + 
  geom_line(size=1) + geom_line(aes(y=lc2,color="LC Spouse2",linetype="LC Spouse2"),size=1) +
  geom_line(aes(y=fc1,color="FC Spouse1",linetype="FC Spouse1"),size=1) +
  geom_line(aes(y=fc2,color="FC Spouse2",linetype="FC Spouse2"),size=1) +
  scale_linetype_manual(values=c("LC Spouse1"=1, "LC Spouse2"=2, "FC Spouse1"=1, "FC Spouse2"=2),name="") +
  scale_colour_manual(values=c("LC Spouse1"="red", "LC Spouse2"="red",
                               "FC Spouse1"="blue","FC Spouse2"="blue"),name="") + 
  ylab("Average of Log Consumption") + xlab("Time") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30),
        legend.key.size = unit(1, 'cm')) 

ggsave('AvLogCons.png')

# graph variance of log consumption over time
ddat <- data.frame(time=1:T,lc1=apply(log(lcc1),1,var),  lc2=apply(log(lcc2),1,var), 
                   fc1=apply(log(fcc1),1,var),  fc2=apply(log(fcc2),1,var) )

ggplot(ddat,aes(x=time,y=lc1,color="LC Spouse1",linetype="LC Spouse1")) + 
  geom_line(size=1) + geom_line(aes(y=lc2,color="LC Spouse2",linetype="LC Spouse2"),size=1) +
  geom_line(aes(y=fc1,color="FC Spouse1",linetype="FC Spouse1"),size=1) +
  geom_line(aes(y=fc2,color="FC Spouse2",linetype="FC Spouse2"),size=1) +
  scale_linetype_manual(values=c("LC Spouse1"=1, "LC Spouse2"=2, "FC Spouse1"=1, "FC Spouse2"=2),name="") +
  scale_colour_manual(values=c("LC Spouse1"="red", "LC Spouse2"="red",
                               "FC Spouse1"="blue","FC Spouse2"="blue"),name="") + 
  ylab("Variance of Log Consumption") + xlab("Time") + 
  theme(plot.title=element_text(size=30),
        axis.title=element_text(size=30),
        axis.text=element_text(size=30),
        strip.text =element_text(size=30),
        legend.text=element_text(size=30),
        legend.key.size = unit(1, 'cm')) 
ggsave('VarLogCons.png')
