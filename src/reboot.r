# Cargar framework
source("./src/configuration.R")
source("./src/plot_functions.R")

fly.rnd <- read.csv(
  "./data/stem15/fly.random.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.bound <- read.csv(
  "./data/stem15/fly_bound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.nonbound <- read.csv(
  "./data/stem15/fly_nonbound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

# First random file
fly.repressed <- read.csv(
  "./data/stem15/fly_repressed.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.bound$Tipo <- "unidos"
fly.nonbound$Tipo <- "no unidos"
fly.repressed$Tipo <- "reprimidos"
fly.rnd$Tipo <- "Random"

fly.bound$AccessionID <- as.character(fly.bound$AccessionID)
fly.nonbound$AccessionID <- as.character(fly.nonbound$AccessionID)
fly.repressed$AccessionID <- as.character(fly.repressed$AccessionID)


aux <-
  fly.repressed %>% select (AccessionID) %>% unique() %>% inner_join(fly.bound %>% select(AccessionID) %>% unique(),
                                                                     by =
                                                                       c("AccessionID"))

fly.desest <- fly.bound %>% inner_join(aux, by = c("AccessionID"))
fly.desest$Tipo <- "desestabilizado"

fly.bound <- fly.bound %>% anti_join(aux, by = c("AccessionID"))
fly.repressed <-
  fly.repressed %>% anti_join(aux, by = c("AccessionID"))

aux <- NULL
modo <- 2 # desest vs random / todos vs todos / bound vs random ?

if (modo == 1) {
  fly <- rbind(fly.rnd, fly.desest)
} else {
  fly <- rbind(fly.bound, fly.desest)
  fly <- rbind(fly, fly.repressed)
  fly <- rbind(fly, fly.nonbound)
  fly <- rbind(fly, fly.rnd)
  
}
fly$Tipo <- as.factor(fly$Tipo)

# revisar las proporciones informadas en cada set
featuresOfThePlot(fly) # ok 1
geneHitsByPattern(fly) # ok 1 
locationAndPairs(fly) # ok 1
plotLoopPatterns(fly) # ok 1
plotNBases(fly) # ok 1
wpbPlots(fly) # ok 1
basePercentsPlot(fly) # no se
pairsMFEPlot(fly) #boo
purinePlot(fly)
relativePositionPlot(fly)
locationPlot(fly)
