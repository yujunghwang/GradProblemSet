# Pset 5 solution
rm(list=ls()) # clean memory
# set the right directory
library(rstudioapi)
rootdir <- dirname(getSourceEditorContext()$path)
setwd(rootdir)

# load libraries
#install.packages("pracma") # install packages if you haven't installed these yet.
#install.packages("maxLik")
#install.packages("plm")
library(pracma)
library(maxLik)
library(stats)
library(plm)

source("fullLogLike.R")
source("predictCCP.R")
source("ccpLogLike.R")

#--------------------#
# Question 1.

# read data
data <- read.csv("Pset5DataQ1.csv",header=TRUE)

x <- data$x
d <- data$d

# set known parameters
N <- 500 # number of buses
T <- 20  # number of periods per bus
M <- 9   # maximum mileage
beta <- 0.95 # discount factor
euler <- 0.5772 # Euler constant

# Question 1.
# (a) estimate parameter theta1, theta2, using Rust full solution method
oout <- maxLik(fullLogLike,start=c(0,0))
summary(oout)
coef(oout)  # estimate
stdEr(oout) # standard error

# (b) estimate parameter theta1, theta2, using AM's 1-step NPL.
# 1st stage : compute CCP
p1 <- rep(N,(M+1))
for (k in 1:(M+1)){
  p1[k] <- mean(d[x==(k-1)])
}
p0 <- 1-p1

# 2nd stage : ML estimation
oout <- maxLik(ccpLogLike,start=c(0,0))
summary(oout)
coef(oout)  # estimate
stdEr(oout) # standard error

Est1st <- coef(oout)

# (c) estimate parameter theta1, theta2 using AM's 2-step NPL.
# update CCP using 1st stage estimate

oout <- predictCCP(param=Est1st,p0=p0,p1=p1)
p0New <- oout[[1]] # NT by 1 
p1New <- oout[[2]] # NT by 1

p1 <- rep(N,(M+1))
for (k in 1:(M+1)){
  p1[k] <- mean(p1New[x==(k-1)])
}
p0 <- 1-p1 ### compute grid of length M+1

# 2nd stage
oout <- maxLik(ccpLogLike,start=Est1st)
summary(oout)
coef(oout)  # estimate
stdEr(oout) # standard error
Est2nd <- coef(oout)


#-----------------#
# Question 3.
rm(list=ls())
data <- read.csv("Pset5DataQ3.csv",header=TRUE)

source("ccpLogLikePH.R")
source("predictCCPPH.R")


set.seed(11215)

# set known parameters
N <- max(data$id)    # number of buses
T <- max(data$time)  # number of periods per bus
M <- max(data$x)     # maximum mileage
beta <- 0.95 # discount factor
euler <- 0.5772 # Euler constant

S <- 2 # number of unobs type

# copy observations
tx <- kron(rep(1,S),data$x)
td <- kron(rep(1,S),data$d)
ttype <- kron(seq(1,S,1),rep(1,N*T))

# estimate fixed effect regression for initialization
fereg <- plm(d~factor(x),data=data,index=c("id","time"),model="within")
feest <- fixef(fereg)[1:N]
oout <- kmeans(feest,S,nstart=30)
label <- order(oout$centers)
inigp <- rep(NA,N)
for (k in 1:S){
  inigp[oout$cluster==k] <- label[k] ### ordering
}
inigp <- rep(inigp, T*S)

# set qst 
tqst      <- rep(NA,N*T*S)
for (k in 1:S){
 tqst[ttype==k] <- 100*as.integer(inigp[ttype==k]==k) +   runif(sum(ttype==k),1)*(1-as.integer(inigp[ttype==k]==k))
}
dim(tqst) <- c(N*T,S)
tqst <- tqst/matrix(rep(apply(tqst,1,sum),S),ncol=S)
iniprob <- apply(tqst,2,mean)
dim(tqst) <- c(N*T*S,1)

# prepare regressors
xx <- c(as.integer(ttype==2))
for (k in 1:M){
  xx <- cbind(xx,as.integer(tx==k))
}
for (k in 1:M){
  xx <- cbind(xx,as.integer(tx==k)*as.integer(ttype==2))
}

# start EM algorithm
iter <-1
tol  <-10^-3
maxiter <-200
dif  <-1

# estimate in two steps
# first stage EM
while (iter<=maxiter & dif>tol){
  
  otqst <- tqst
  
  # compute reduced form CCP using weighted logit
  aa   <-glm(td~xx, family=quasibinomial(link="logit"),weights=tqst,control=glm.control(trace=FALSE,maxit=2000))
  bb   <-predict(aa,type="response")  
  like  <- bb*(td==1) + (1-bb)*(td==0)
  rm(aa,bb)  
  
  # update tqst
  dim(like)<-c(N,T,S)
  like <- apply(like,c(1,3),prod)
  for (k in 1:S){
    like[,k] <- iniprob[k]*like[,k]
  }
  den <- apply(like,1,sum)
  den <- matrix(rep(den,S),ncol=S)
  qq <- like/den # N by S
  
  for (k in 1:S){
    tqst[((k-1)*N*T+1):(k*N*T)] <- rep(qq[,k],T)
  }

  # update initial type prob
  iniprob <- apply(qq,2,mean)
    
  iter <- iter+1
  dif <- max(abs(tqst-otqst))
  print(dif)
  
}

# estimate the flow payoff parameters
# estimated ccp
p0 <- array(rep(NA,(M+1)*S),c(M+1,S))
for (k1 in 1:(M+1)){
  for (k2 in 1:2){
    iind <- tx==(k1-1) & ttype==k2 
    p0[k1,k2] <- sum(tqst[iind==1]*(1-td[iind==1]))/sum(tqst[iind==1])
  } 
}
p1 <- 1-p0

# correction term 
txs  <- pmin(tx+1,M)*ttype
# compute once before ML to reduce computational time
crrn <- rep(NA,N*T*S)
for (k in 1:(N*T*S)){
  crrn[k] <-  beta*(log(1+p0[1,ttype[k]]/p1[1,ttype[k]]) - log(1+p0[pmin(tx[k]+2,M),ttype[k]]/p1[pmin(tx[k]+2,M),ttype[k]]) )
}

oout <- maxLik(ccpLogLikePH,start=c(0,0))
summary(oout)
coef(oout)  # estimate
stdEr(oout) # standard error
