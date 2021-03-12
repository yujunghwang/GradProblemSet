# prepare matrices

# indicator for working
work0   <- matrix(c(0,1))
work1_i <- matrix(c(0,1,0,1))
work1_j <- matrix(c(0,0,1,1))

# Pre-calculate household budget constraint
c0_i <- matrix(0, nrow = 2, ncol = I)
c0_j <- matrix(0, nrow = 2, ncol = J)
c1   <- array(0, dim = c(4, I, J))

for(ixI in 1:I) {
  c0_i[ , ixI] <- y + wage_i[ixI] * work0 * (1 - tax0_i)
}

for(ixJ in 1:J) {
  c0_j[ , ixJ] <- y + wage_j[ixJ] * work0 * (1 - tax0_j)
}

for(ixJ in 1:J) {
  for(ixI in 1:I) {
    c1[ , ixI, ixJ] <- 2*y + wage_i[ixI] * work1_i * (1 - tax1_i) +
      wage_j[ixJ] * work1_j * (1 - tax1_j)
  }
}

logc0_i <- log(c0_i)
logc0_j <- log(c0_j)
logc1   <- log(c1)

#### Calculate expected values when single
work0_pr_i <- exp((logc0_i - beta*matrix(rep(work0,I),ncol=I))/esd)
for (ixI in 1:I){
  work0_pr_i[,ixI] <- work0_pr_i[,ixI]/sum(work0_pr_i[,ixI])
}

work0_pr_j <- exp((logc0_j - beta*matrix(rep(work0,J),ncol=J))/esd)
for (ixJ in 1:J){
  work0_pr_j[,ixJ] <- work0_pr_j[,ixJ]/sum(work0_pr_j[,ixJ])
}

V0_i <- apply((logc0_i - beta*matrix(rep(work0,I),ncol=I))*work0_pr_i, 2, sum)
V0_j <- apply((logc0_j - beta*matrix(rep(work0,J),ncol=J))*work0_pr_j, 2, sum)


