# Cargar framework
source("./src/configuration.R")
source("./src/plot_functions.R")

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
fly.bound$AccessionID <- as.character(fly.bound$AccessionID)
fly.nonbound$AccessionID <- as.character(fly.nonbound$AccessionID)
fly.repressed$AccessionID <- as.character(fly.repressed$AccessionID)

glimpse(aux)

aux <-
  fly.repressed %>% select (AccessionID) %>% unique() %>% inner_join(fly.bound %>% select(AccessionID) %>% unique(), by =
                                                               c("AccessionID"))

fly.desest <- fly.bound %>% inner_join(aux, by=c("AccessionID")) 
fly.desest$Tipo <- "desestabilizado"

fly.bound <- fly.bound %>% anti_join(aux, by=c("AccessionID"))
fly.repressed <- fly.repressed %>% anti_join(aux, by=c("AccessionID"))

fly <- rbind(fly.bound, fly.desest)
fly <- rbind(fly, fly.repressed)

# revisar las proporciones informadas en cada set
featuresOfThePlot(fly) #boo
geneHitsByPattern(fly) #boo
chromoPattAndPairs(fly) #boo
plotNBases(fly) #ok?
wpbPlots(fly) #boo
basePercentsPlot(fly) #shit
pairsMFEPlot(fly) #boo
purinePlot(fly)
relativePositionPlot(fly)
locationPlot(fly)
