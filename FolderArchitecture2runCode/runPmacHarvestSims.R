#This runs sperm whale simulations, with harvest

source("Functions/runPopSims.R")
Sp <- "Pmac"

#Define number of iterations
# Warning: starting with a small number of iterations is highly recommended.
nsims <- 10

#Number of years to run model forward
# 1788 was the first year of documented harvest in the Reeves paper
# Their effort plot is decadal, starting in 1785, so we will start from then, 
#   and assume whaling started that year
# DWH oil spill was in 2010, so run through to 2010.
# Total run, therefore is from 1785 to 2010 - 226 years
nyears <- 226 

#Whaling effort by decade, from figure 3 of Reeves paper - sum to 214
whaling.effort.bydecade <- c(18,20,2,2,20,25,77,20,20,10)
#Assume equally spread over the decade
whaling.effort.byyear <- rep(whaling.effort.bydecade/10, each = 10)
#Add in the other years
whaling.effort.byyear <- c(whaling.effort.byyear, numeric(nyears - length(whaling.effort.byyear)))
#Standardize
whaling.effort.byyear <- whaling.effort.byyear / sum(whaling.effort.byyear)
#Total of 1070 whales estimated to have been taken
harvest <- 1070
#Divide into years
harvest <- harvest * whaling.effort.byyear

hScenario <- 2

#Run simulations
runPopSims(Sp = Sp, nsims = nsims, nyears = nyears, seed = 7134672,
           harvest = harvest, hScenario = hScenario)

#Load results
SpInfo <- getSpData(Sp)
load(file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, "Sim.RData"))

poptraj = matrix(NA, ncol = nyears, nrow = nsims)
for(i in 1:nsims) {
  poptraj[i, ] = colSums(simres[, , i, 2])
}
ylims <- range(poptraj)
plot(poptraj[i, ], type = "n", ylim = ylims, xlab = "Year", ylab = "Predicted population size", las = 1)
for(i in 1:nsims){
  lines(poptraj[i, ], type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
}
lines(colMeans(poptraj), type = "l", lwd = 3, col = "#1b2ac8")  

