# make fake data for pset 4

rm(list=ls()) # clean memory
# set the right directory
# setwd("PATH")

# load libraries
#install.packages("pracma") # install packages if you haven't installed these yet.
#install.packages("ggplot2")
library(pracma)
library(ggplot2)
library(MASS)

# set seed
set.seed(191121)

# set structural parameters
varp  <- 0.03
varm  <- 0.02
covpe <- 0.01
delta <- 0.13
lam0  <- 0.02
lam1  <- -0.005

nsim <- 3000
T    <- 8

# generate shocks
shocks <- as.data.frame(mvrnorm(n=nsim*T,mu=c(0,0),Sigma=rbind(c(varp,covpe),c(covpe,1))))
colnames(shocks) <- c("income","work")

# generate id
id <- kron(rep(1,T),seq(1,nsim,1))

# generate wave
wave <- kron(seq(1,T,1), rep(1,nsim))

# generate age
age0 <- randi(imax=30,n=nsim,m=1) + 20
age  <-c()
for (k in 1:T){
  age <- c(age, age0+(k-1))
}
rm(age0)
  
# generate education
educ <- kron(rep(1,T),randi(4,n=nsim,m=1))

# instrument
noveliv <- matrix(randn(n=nsim,m=T),ncol=1)

work  <- ( 1.3+ 0.44*age -0.01*age^2 + 0.2*(educ==2) + 0.3*(educ==3) + 0.5*(educ==4) + 0.8*noveliv + shocks$work > 0 )

# stochastic income
pcomp0 <- matrix(rep(randn(n=nsim,m=1),T),ncol=T)
temp <- matrix(shocks$income,ncol=T)
temp <- t(apply(temp,1,cumsum))
pcomp <- pcomp0 + temp
dim(pcomp) <- c(nsim*T,1)
rm(temp,pcomp0)
mcomp  <- randn(n=nsim*T,m=1)*sqrt(varm)

#h0 <- rep(0.2*as.numeric(educ[1:nsim]) + randn(nsim,1),T)
#h0 <- rep(randn(nsim,1),T)
#dhc <- -delta*matrix(1-work,ncol=T) + (lam0 + lam1*matrix(age,ncol=T))*matrix(work,ncol=T)
#dhc <- t(apply(dhc,1,cumsum))
#dim(dhc) <- c(nsim*T,1)

#lninc <- 0.2 + 0.3*(educ==2) + 0.4*(educ==3) + 0.6*(educ==4) + h0 + dhc + pcomp + mcomp
lninc <- 0.2 + 0.4*age -0.02*age^2  + 0.3*(educ==2) + 0.4*(educ==3) + 0.6*(educ==4) + pcomp + mcomp

lninc[work==0] <- NA

educ <- factor(educ,labels=c("High School Graduate","Some College","College Graduate","Postgraduate"))

simData <- data.frame(id=id,wave=wave,age=age,educ=educ,noveliv=noveliv,lninc=lninc,work=work)

write.csv(simData,file="Pset4Data.csv")

