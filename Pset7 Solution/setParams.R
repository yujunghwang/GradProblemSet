
#### Set parameters ####
interpMethod <- "linear"
tol     <- 10^-10    # maximum allowed error
minCons <- 10^-5     # minimum allowe consumption
T      <- 10
r      <- 0.03
beta   <- 0.95
gamma  <- 1.5
borrowingAllowed <- 1
isUncertainty <- 1
startA <- 0

# Grids
numPointsA <-10
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


