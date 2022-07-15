# Function info:
# This function creates the required Population Projection Matrix sometimes 
# referred to as A : see Caswell 2018 page 10

# Based on the required components:
# Ps - probabilities of remaining in stage (i.e. survive and not changing stage)
# Gs - probabilities of moving to next stage
# Bs - production of newborns

# Note that the function is partly hardwired to work with 16 stages
# when static numbers are used to alocate the Gs
# this might require updates later for other species not the sperm whale


buildPPM <- function(Ps, Gs, Bs) {
  
  # define number of stages
  nstages <- length(Ps)
  A <- matrix(0, nstages, nstages)
  
  # fill in the main diagonal with probabilities of remaining in stage
  diag(A) <- Ps
  
  #
  b <- Bs[1]
  be <- Bs[2]
  
  # generating girl babies
  A[1, 3] <- b  # from unexposed females
  A[1, 8] <- be # from exposed females
  
  # generating boy babies
  A[11, 3] <- b      # from unexposed females
  A[11, 8] <- be     # from exposed females
  A[3, 5]  <- Gs[5]  # unexposed post calving becomes mature, again
  A[8, 10] <- Gs[10] # exposed post calving becomes mature, again
  
  # unexposed female calf, juvenile, mature move on to next stage
  for (i in 2:4) {
    A[i, i - 1] <- Gs[i - 1]
  }
  # exposed female calf, juvenile, mature move on to next stage
  for (i in 7:10) {
    A[i, i - 1] <- Gs[i - 1]
  }
  # unexposed male calf, juvenile move on to next stage
  for (i in 12:13) {
    A[i, i - 1] <- Gs[i - 1]
  }
  # exposed male calf, juvenile move on to next stage
  for (i in 15:16) {
    A[i, i - 1] <- Gs[i - 1]
  }
  
  # exposed mother becomes post calving
  A[5, 4] <- Gs[4]
  # unexposed mother becomes post calving
  A[10, 9] <- Gs[9]
  
  # return Population Projection Matrix
  return(A)
}
