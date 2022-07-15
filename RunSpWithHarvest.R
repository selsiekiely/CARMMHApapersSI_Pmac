#LT: Based on Tiago's RunSp, I have created RunSpWithHarvest that simply
# subtracts off a fixed harvest after projecting forward each year

# Function info:
# this function runs a simulation forward for a given species and number of years
# assuming a stage structured model
# considering, or not, there was an oil spill

# Arguments:
# sp    - a list with all the relevant species details
# years - a scalar representing the number of years to run the model for
# pe    - a scalar that represents the proportion of the population exposed to oil.
#         By default pe is 0, i.e., project forward under a no oil scenario
# harvest - a matrix with years rows and nstage columns containing the number of animals
#           in each year and stage to subtract after projecting forward.

runSpWithHarvest <- function(sp.list, years, pe = 0, harvest = NULL) {
  
  # create objects to hold results
  # p holds the population size per stage per year
  nstage <- length(sp.list$survs)
  p <- mat.or.vec(years, nstage)

  #create a harvest vector if null passed in
  if(is.null(harvest)) harvest <- matrix(0, years, nstage)
    
  # These were in Lance's code, but since they were not looked at I've commented them out
  # Interbirth interval, per year
  # ib.store <- rep(0, years)
  # ib.store[1] <- sp.list$ibnom
  # Lambda - growth rate for population, per year
  # lambda.years <- rep(0, years)
  
  # get survival probabilities
  sigma      <- sp.list$survs
  sigma.base <- sp.list$survs
  
  # get maximum and nominal fecundity
  f.max <- 1 / sp.list$ibmin
  f.nom <- 1 / sp.list$ibnom
  
  #Density dependent model to influence interbirth intervals
  beta <- (1 / sp.list$totn) * ((f.max - f.nom) / f.nom) ^ (1 / sp.list$rho)
  
  
  # b is birth rate for calves from unexposed moms
  b <- 0.5 * 1 / sp.list$ibnom
  # get initial population size per stage
  stapop <- GetStablePop(sp.list)$stable.dist
  # pop is females unexposed, females exposed, males unexposed, males exposed
  p[1, ] <- c(stapop[1:5] * (1 - pe), stapop[1:5] * pe, stapop[6:8] * (1 - pe), stapop[6:8] * pe)
  #harvest at time 1
  p[1, ] <- p[1, ] - harvest[1, ]
  
  # get reduced values for survival and increased mortalities, if oil spill ocurred
  # yearly baseline values for added mortality  - 0 means no added mortality
  mort.vector <- rep(0, years)
  # yearly base line values for reduced reproduction - 1 means no reduced reproduction
  repro.vector <- rep(1, years)
  
  if (pe != 0) {
    
    # if there was a non-0 proportion pe of animals exposed to oil
    # ---
    # both durations of effect of oil (in reproduction and survival) were
    # unstated elsewhere - unless I missed it - so these values are only
    # available directly from the codee, with no details about how they
    # were calculated/obtained
    # ---
    # For the perturbed model, create vectors for annual mortality and reproductive 
    # rates for the exposed cohort
    
    # how long are effects constant
    mc <- sp.list$mc # mortality
    rc <- sp.list$rc # reproduction
    
    # how long they take to return to baseline
    mc2b <- sp.list$mc2b # mortality
    rc2b <- sp.list$rc2b # reproduction
    
    # increased mortality
    mort.effect       <- sp.list$me
    mort.vector[1:mc] <- mort.effect
    
    # note this implies that the effect of the oil lasts only 15 years
    for (i in 1:mc2b) {
      mort.vector[mc + i] <- mort.effect - ((1 / mc2b) * i * mort.effect)
    }
    
    # ---
    
    # reduced fertility
    repro.effect <- sp.list$re
    repro.vector[1:rc] <- 1 - repro.effect
    # note this implies that the effect of the oil lasts 20 years
    
    for (i in 1:rc2b) {
      repro.vector[rc + i] <- 1 - (repro.effect - ((1 / rc2b) * i * repro.effect))
    }
  }
  
  # before the oil, the ib of exposed animals is the same as unexposed
  ib.yr <- sp.list$ibnom
  
  # for each year one if propagating the population forward for
  for (y in 2:years) {
    
    # apply year specific added mortality perturbation to exposed
    # note if pe=0 then mort.vector is all 0's, so nothing happens
    sigma[6:10]  <- sigma.base[6:10] * (1 - mort.vector[y])
    sigma[14:16] <- sigma.base[14:16] * (1 - mort.vector[y])
    
    # get the birth induced changes by oil
    b <- 0.5 * 1 / ib.yr
    # note if pe=0 then repro.vector is all 1's, so nothing happens
    be <- b * repro.vector[y]
    
    # get the gamma's
    # this bit I do not understand - how do you get gammas
    lambda.start <- 1
    gamma.test <- (((sigma / lambda.start) ^ sp.list$tinit) - ((sigma / lambda.start) ^ (sp.list$tinit - 1))) / (((sigma / lambda.start) ^ sp.list$tinit) - 1)
    # get gamma for mature females
    gamma.test[3] <- b / (0.5 * sigma[3] * sqrt(sigma[4]))
    gamma.test[8] <- be / (0.5 * sigma[8] * sqrt(sigma[9]))
    
    # get the Ps and the Gs required for the Population Projection Matrix
    # i.e. conditional on survival, who stays in state (P) and who grows out of it (G)
    P.test <- sigma * (1 - gamma.test)
    G.test <- sigma * gamma.test
    # build Population Projection Matrix
    A.test <- buildPPM(P.test, G.test, c(b, be))
    # set A to A.test
    A.yr <- A.test
    # project the population forward
    p[y,] <- A.yr %*% p[y - 1,]
    #implement harvest
    p[y, ] <- p[y, ] - harvest[y, ]
    # now calculate ib for next iteration
    # get population size
    N.yr <- sum(p[y,])
    # corresponding fecundity
    # this is not in Lance's code... but should it be?
    # beta <- (1 / N.yr) * ((f.max - f.nom) / f.nom) ^ (1 / rho)
    f.yr <- f.max / (1 + (beta * N.yr) ^ sp.list$rho)
    # corresponding ib
    ib.yr <- 1 / f.yr
    # ib.store[y] <- ib
  }
  
  # get the final results
  return(p)
}
