# Function info:
# This function reads in the information in a file to create a species
# object which can then be used by runSp.R
# The file must have a specific fixed format

getSpData <- function(sp.code, file = "SpeciesDefinitionFile.xlsx") {
  
  # the argument for this function is the species code
  # sp: a 4 letter character, e.g. Pmac for sperm whale
  require(readxl)
  SDF    <- read_excel(file)
  Spinfo <- SDF[SDF$Species == sp.code,]
  Sp     <- list(
    
    # Total population size
    totn <- Spinfo$N,
    # minimum inter-birth interval
    ibmin <- Spinfo$ibmin,
    # inter-birth interval resulting in a stable population
    ibnom <- Spinfo$ibnom,
    # f.max <- 1 / ib.min
    # f.nom <- 1 / ib.nom
    # proportion of the population exposed to oil
    pe <- Spinfo$pe,
    
    # Density dependence rho parameter
    # model to influence interbrith intervals
    rho <- Spinfo$rho,
    # stage specific survival rates
    survs <- as.numeric(Spinfo[, (3 + 15 + 1):(3 + 15 + 1 + 15)]),
    # expected years in stage
    tinit <- as.numeric(Spinfo[, 3:(3 + 15)]),
    
    #----------------------------------------------------
    
    # added mortality due to oil
    me <- Spinfo$me,
    # years it remains constant
    mc <- Spinfo$mc,
    # years it takes to get back to baseline after constant
    mc2b <- Spinfo$mc2b,
    
    #----------------------------------------------------
    
    # lowered reproduction due to oil
    re <- Spinfo$re,
    # years it remains constant
    rc <- Spinfo$rc,
    # years it takes to get back to baseline after constant
    rc2b <- Spinfo$rc2b
  )
  names(Sp) <- c("totn", "ibmin", "ibnom", "pe", "rho", "survs", "tinit", "me", "mc", "mc2b", "re", "rc", "rc2b")
  return(Sp)
}

# Pmac <- getSpData("Pmac")
