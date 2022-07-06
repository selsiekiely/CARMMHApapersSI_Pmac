#This runs sperm whale simulations, with harvest


source("Functions/runPopSims.R")
source("Functions/reqfuns.R")
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
#Load up species info into list
SpInfo <- getSpData(Sp)

#Scenario 1 - best guess scenario
hScenario <- 1
runPopSims(Sp = Sp, nsims = nsims, nyears = nyears, seed = 7134672,
           harvest = harvest, hScenario = hScenario)
load(file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, "Sim.RData"))
s1res = matrix(NA, ncol = nyears, nrow = nsims)
for(i in 1:nsims) {
  s1res[i, ] = colSums(simres[, , i, 2])
}

#Scenario 2 - best case (smallest take, biggest pop size, proportional harvest)
hScenario <- 1
#Create larger population size, temporarily
# Note - best to back up the N_boot.csv file!
Nfile <- paste0("InOutBySp/", SpInfo$folder, "/N_boot.csv")
Nstart <- read.csv(Nfile, header = TRUE)
Nstart.larger <- Nstart * (1 / 0.35)
#Replace N_boot file, temporarily
write.csv(Nstart.larger$x, Nfile)
#Run sim
runPopSims(Sp = Sp, nsims = nsims, nyears = nyears, seed = 7134672,
           harvest = harvest, hScenario = hScenario)
#Put back old N_boot numbers
write.csv(Nstart$x, Nfile)
load(file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, "Sim.RData"))
s2res = matrix(NA, ncol = nyears, nrow = nsims)
for(i in 1:nsims) {
  s2res[i, ] = colSums(simres[, , i, 2])
}

#Scenario 3 - worst case (largest take, smallest pop size, harvest more on females)
hScenario <- 2
runPopSims(Sp = Sp, nsims = nsims, nyears = nyears, seed = 7134672,
           harvest = harvest * 2, hScenario = hScenario)
load(file=paste0("InOutBySp/", SpInfo$folder, "/", Sp, "simres", nsims, "Sim.RData"))
s3res = matrix(NA, ncol = nyears, nrow = nsims)
for(i in 1:nsims) {
  s3res[i, ] = colSums(simres[, , i, 2])
}

#plot the results
par(mfrow = c(3, 1))

ylims <- range(c(s1res))
plot(s1res[i, ], type = "n", ylim = ylims, xlab = "Year", ylab = "Predicted population size", las = 1, main = "Best guess")
for(i in 1:nsims){
  lines(s1res[i, ], type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
}
lines(colMeans(s1res), type = "l", lwd = 3, col = "#1b2ac8")  

ylims <- range(c(s2res))
plot(s2res[i, ], type = "n", ylim = ylims, xlab = "Year", ylab = "Predicted population size", las = 1, main = "Best case")
for(i in 1:nsims){
  lines(s2res[i, ], type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
}
lines(colMeans(s2res), type = "l", lwd = 3, col = "#1b2ac8")  

ylims <- range(c(s3res))
plot(s3res[i, ], type = "n", ylim = ylims, xlab = "Year", ylab = "Predicted population size", las = 1, main = "Worst case")
for(i in 1:nsims){
  lines(s3res[i, ], type = "l", lwd = 0.7, col = rgb(0, 0, 0, 0.15))
}
lines(colMeans(s3res), type = "l", lwd = 3, col = "#1b2ac8")  

par(mfrow = c(1, 1))
