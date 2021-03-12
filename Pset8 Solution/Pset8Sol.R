# Problem Set 8
# written by Yujung Hwang
# building on a code provided by Andrew Shephard
rm(list=ls())

# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

library(devtools)

source("MMclear.R")
source("supplyMen.R")
source("demandMen.R")
source("ccp_ij.R")
source("marval_i.R")
source("marval_j.R")

###### set parameters
I <- 3 # num of type of men
J <- 3
wage_i <- c(1,2,3)
wage_j <- c(1,2,3)
measure_i <- c(1/3,1/3,1/3)
measure_j <- c(1/3,1/3,1/3)
beta <- 0.8 # leisure value
esd  <- 0.4 # time allocation state specific sd
tsd  <- 0.7 # choo-siow marriage sd
y    <- 1 # non labor income

#### (c)
# tax schedule
tax0_i <- 0 # tax rate for single men
tax0_j <- 0 # tax rate for single women
tax1_i <- 0 # tax rate for married men
tax1_j <- 0 # tax rate for married women

# solve single problem
source("prepareMatrices.R")

objfun <- function(lambda){
  return(sum(MMclear(lambda)^2))
}

lambda <- optim(par=rep(0.5,I*J),fn=objfun,control=list(abstol=10^-25))$par
dim(lambda) <- c(I,J) # I by J matrix
print(lambda) ### pareto weight
print(supplyMen(lambda)) ### marriage matching function
print((measure_i - apply(supplyMen(lambda),1,sum))/measure_i ) ## prob of single men of each type


#### (d)
tax0_i <- 0.5 # tax rate for single men
tax0_j <- 0 # tax rate for single women
tax1_i <- 0 # tax rate for married men
tax1_j <- 0 # tax rate for married women

# solve single problem
source("prepareMatrices.R")
lambda <- optim(par=rep(0.5,I*J),fn=objfun)$par
dim(lambda) <- c(I,J) # I by J matrix
print(lambda) ### pareto weight
print(supplyMen(lambda)) ### marriage matching function
print((measure_i - apply(supplyMen(lambda),1,sum))/measure_i ) ## prob of single men of each type
