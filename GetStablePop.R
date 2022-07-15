# Function info:
# This function uses a species current population size and population parameters
# and gets the correct age structure for a stable population size.

# Unlike before, the notion of exposed animals was removed from 
# the stable structure population as the stable age structure 
# for exposed cohort would otherwise be all zeros-because there
# is no production of new exposed animals.
# Nothing generates new either exposed females or males
# So later we also use the stable age structure of the unexposed as the
# starting point for exposed.

# Inputs:
#  - sp:          the species details
#  - ls:          aka lambda.start - the growth rate to start aiming for
#  - convergence: the criterium for convergence
#                 (defines how close to 1 the observed population
#                 growth needs to be to achieve convergence)

GetStablePop <- function(sp.list, ls = 1, converge = 0.0000001) {
  
  # get number of states, from sp object
  n.stage <- length(sp.list$survs)
  # create object to hold the probability of moving out of stage
  gt <- rep(n.stage / 2, 0)
  # value for testing if convergence was achieved
  test <- 1
  # iteration counter
  itcount <- 1
  # b is birth rate for calves from unexposed moms
  b <- 0.5 * 1 / sp.list$ibnom
  # index for non-exposed animals
  ine <- c(1:5, 11:13)
  
  repeat {
    # Given lambda start and sigma, solve for gammas except for 3 and 8
    # i.e. except for mature females
    surs <- sp.list$survs[ine]
    T.inits <- sp.list$tinit[ine]
    gt <- (((surs / ls) ^ T.inits) - ((surs / ls) ^ (T.inits - 1))) / (((surs / ls) ^ T.inits) - 1)
    
    # get gamma for mature females
    gt[3] <- b / (0.5 * surs[3] * sqrt(surs[4]))
    
    # build A.test
    # note this means, of those that survive, sigma, some persist in the same stage, P, some grow out of it, G
    P.test <- surs * (1 - gt)
    G.test <- surs * gt
    # create a matrix to hold values
    A.test <- mat.or.vec(n.stage / 2, n.stage / 2)
    # the probability of surviving and staing in same stage
    diag(A.test) <- P.test
    # new girl babies
    A.test[1, 3] <- b
    # new boy babies
    A.test[6, 3] <- b
    # post-calving females that go back to mature
    A.test[3, 5] <- G.test[5]
    
    # for females
    # calves becoming juveniles, juveniles becoming mature, mature becoming moms - unexposed
    for (i in 2:4) {
      A.test[i, i - 1] <- G.test[i - 1]
    }
    # for males
    # calves becoming juveniles, juveniles becoming mature - unexposed
    for (i in 7:8) {
      A.test[i, i - 1] <- G.test[i - 1]
    }
    # the probability of a mother becoming post calving
    # unexposed
    A.test[5, 4] <- G.test[4]
    
    # Calculates the population growth rate of the current projection matrix
    lambda.test <- lambda(A.test)
    
    # calculate the diference between estimated growth rate of the current population's 
    # age structure if close enough to 1 that means, by definition, a stable population
    test <- abs(lambda.test - ls)
    # DO NOT PRINT CONVERGENCE STATS
    if (test < converge) {
      #  print("converge")
      #  print(paste("absolute diference between start and end population growth = ", test))
      #  print(paste("iterations required = ", itcount))
      #  print(paste("population growth = ", lambda.test))
      break
    }
    # gets here only if you haven't converged yet
    # print(lambda.test)
    ls <- lambda.test
    itcount <- itcount + 1
  }
  # once converged, i.e. age structure is such that the population is stable
  # get the stable age
  stable.age <- eigen.analysis(A.test)$stable.stage
  # get the population size by stage
  n.init <- sp.list$totn * stable.age
  # check the final N is indeed the starting value
  N <- sum(n.init)
  # returns a list containing
  # stable.dist:  the number per stage
  # N:       the population size, should equal sp$N.pop
  # trans.mat:  the transition matrix used
  # pop.growth: the growth rate at convergence, should be close to 1 (that is a stable population, by definition)
  return(list(
    stable.dist = n.init,
    N           = N,
    trans.mat   = A.test,
    pop.growth  = lambda.test
  ))
}
# Test run with sperm whale
# PmacStablePop=GetStablePop(Pmac)
# KospStablePop=GetStablePop(Kosp)

