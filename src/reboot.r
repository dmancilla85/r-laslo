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
  "./data/stem15/fly.bound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

#nounidos <- read.table("./data/fly_non_bound.txt",header = T,stringsAsFactors = F)
#reprimidos <- read.table("./data/fly_repressed.txt",header = T,stringsAsFactors = F)
nounidos <- read.table("./data/fly_non_bound.lst",header = F,stringsAsFactors = F)
reprimidos <- read.table("./data/fly_repressed.lst",header = F,stringsAsFactors = F)
nounidos <- nounidos %>% anti_join(reprimidos)

write.table(nounidos,"nounidos.lst")

fly.nonbound <- read.csv(
  "./data/stem15/fly.non_bound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

# First random file
fly.repressed <- read.csv(
  "./data/stem15/fly.repressed.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.bound$Tipo <- "Unidos y no silenciados"
fly.nonbound$Tipo <- "No unidos y no silenciados"
fly.repressed$Tipo <- "No unidos y silenciados"
fly.rnd$Tipo <- "Random"


fly.bound$AccessionID <- as.character(fly.bound$AccessionID)
fly.nonbound$AccessionID <- as.character(fly.nonbound$AccessionID)
fly.repressed$AccessionID <- as.character(fly.repressed$AccessionID)

# Tomar intersección entre unidos y silenciados
aux <-
  fly.repressed %>% select (AccessionID) %>% unique() %>% 
  inner_join(fly.bound %>% select(AccessionID) %>% unique(), by = c("AccessionID"))

fly.desest <- fly.bound %>% inner_join(aux, by = c("AccessionID"))
fly.desest$Tipo <- "Unidos y silenciados"

fly.bound <- fly.bound %>% anti_join(aux, by = c("AccessionID"))
fly.repressed <-
  fly.repressed %>% anti_join(aux, by = c("AccessionID"))

# Tomar intersección entre no unidos y silenciados

aux <-
  fly.repressed %>% select (AccessionID) %>% unique() %>% 
  inner_join(fly.nonbound %>% select(AccessionID) %>% unique(), by = c("AccessionID"))

fly.nonbound <- fly.nonbound %>% anti_join(aux, by = c("AccessionID"))

aux <- NULL

modo <- 2 # 1: Unidos y silenciados vs random / 2: Todos vs todos

if (modo == 1) {
  fly <- rbind(fly.rnd, fly.desest)
} else {
  fly <- rbind(fly.bound, fly.desest)
  fly <- rbind(fly, fly.repressed)
  fly <- rbind(fly, fly.nonbound)
  fly <- rbind(fly, fly.rnd)
  
}
fly$Tipo <- as.factor(fly$Tipo)

fly.bound     %>% distinct(Gen) %>% nrow()  # 67
fly.nonbound  %>% distinct(Gen) %>% nrow()  # 161
fly.repressed %>% distinct(Gen) %>% nrow()  # 1171
fly.rnd       %>% distinct(Gen) %>% nrow()  # 233
fly.desest    %>% distinct(Gen) %>% nrow()  # 183

fly.bound     %>% distinct(AccessionID) %>% nrow()  # 174
fly.nonbound  %>% distinct(AccessionID) %>% nrow()  # 393
fly.repressed %>% distinct(AccessionID) %>% nrow()  # 3035
fly.rnd       %>% distinct(AccessionID) %>% nrow()  # 479
fly.desest    %>% distinct(AccessionID) %>% nrow()  # 469

nrow(fly.bound)     # 293
nrow(fly.nonbound)  # 588
nrow(fly.repressed) # 6327
nrow(fly.rnd)       # 762
nrow(fly.desest)    # 1090

summary(fly.bound$LoopPattern)     # 293
summary(fly.nonbound$LoopPattern)  # 588
summary(fly.repressed$LoopPattern) # 6327
summary(fly.rnd$LoopPattern)       # 762
summary(fly.desest$LoopPattern)    # 1090

fly$Id <- 1:nrow(fly)
fly$N10 <- fly$N9
fly$N9 <- fly$N8
fly$N8 <- fly$N7
fly$N7 <- fly$N6
fly$N6 <- fly$N5
fly$N5 <- fly$N4
fly$N4 <- fly$N3
fly$N3 <- fly$N2
fly$N2 <- fly$N1
fly$N1 <- fly$N0
fly$N0 <- NULL

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
