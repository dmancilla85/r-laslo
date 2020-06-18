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

fly.bound     %>% distinct(Gen) %>% nrow()  # 66
fly.nonbound  %>% distinct(Gen) %>% nrow()  # 161
fly.repressed %>% distinct(Gen) %>% nrow()  # 1170
fly.rnd       %>% distinct(Gen) %>% nrow()  # 237
fly.desest    %>% distinct(Gen) %>% nrow()  # 184

fly.bound     %>% distinct(AccessionID) %>% nrow()  # 174
fly.nonbound  %>% distinct(AccessionID) %>% nrow()  # 393
fly.repressed %>% distinct(AccessionID) %>% nrow()  # 3035
fly.rnd       %>% distinct(AccessionID) %>% nrow()  # 479
fly.desest    %>% distinct(AccessionID) %>% nrow()  # 469

nrow(fly.bound)     # 292
nrow(fly.nonbound)  # 586
nrow(fly.repressed) # 6326
nrow(fly.rnd)       # 758
nrow(fly.desest)    # 1083

summary(fly.bound$LoopPattern)     # 292
summary(fly.nonbound$LoopPattern)  # 586
summary(fly.repressed$LoopPattern) # 6326
summary(fly.rnd$LoopPattern)       # 758
summary(fly.desest$LoopPattern)    # 1083

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

####################################################
## ESTADISTICOS
####################################################


myTable <- fly %>%
  select (Tipo, N.2, N.1, N2, N5, N6, N7, N8) %>%
  gather(metric, value, -Tipo) %>%
  filter(value %in% c("A", "G", "C", "U")) %>%
  mutate(metric = str_replace_all(metric, "N.1", "N(-1)")) %>%
  mutate(metric = str_replace_all(metric, "N.2", "N(-2)")) 

#install.packages("Cairo")
#install.packages("cairoDevice")
#install.packages("vcd")
library(corrplot)
library(Cairo)
library(vcd)
library(cairoDevice)

n_2 <- myTable %>% filter(metric=="N(-2)")  
n_2 <- table(n_2$Tipo, n_2$value)
cq <- chisq.test(n_2) # X-squared = 88.471, df = 12, p-value = 9.756e-14

corrplot(cq$residuals, is.cor = FALSE)

n_1 <- myTable %>% filter(metric=="N(-1)")  
n_1 <- table(n_1$Tipo, n_1$value)
cq1 <- chisq.test(n_1) # X-squared = 80.001, df = 12, p-value = 4.126e-12
corrplot(cq1$residuals, is.cor = FALSE,method="square")
vcd::assocstats(n_1)
contrib1 <- 100*cq1$residuals^2/cq1$statistic
round(contrib, 3)
corrplot(contrib1, is.cor = FALSE)

n2 <- myTable %>% filter(metric=="N2")  
n2 <- table(n2$Tipo, n2$value)
cq2 <- chisq.test(n2) # X-squared = 164.23, df = 12, p-value < 2.2e-16
corrplot(cq1$residuals, is.cor = FALSE,method="square")
contrib2 <- 100*cq2$residuals^2/cq2$statistic
corrplot(contrib2, is.cor = FALSE)

n5 <- myTable %>% filter(metric=="N5")  
n5 <- table(n5$Tipo, n5$value)
cq <- chisq.test(n5) # X-squared = 51.813, df = 12, p-value = 6.692e-07
corrplot(cq$residuals, is.cor = FALSE)

n6 <- myTable %>% filter(metric=="N6")  
n6 <- table(n6$Tipo, n6$value)
cq <- chisq.test(n6) # X-squared = 42.266, df = 12, p-value = 3.004e-05
corrplot(cq$residuals, is.cor = FALSE)

n7 <- myTable %>% filter(metric=="N7")  
n7 <- table(n7$Tipo, n7$value)
cq <- chisq.test(n7) # X-squared = 48.854, df = 12, p-value = 2.218e-06
corrplot(cq$residuals, is.cor = FALSE)

n8 <- myTable %>% filter(metric=="N8")  
n8 <- table(n8$Tipo, n8$value)
cq <- chisq.test(n8) # X-squared = 39.747, df = 12, p-value = 7.92e-05 *may be incorrect
corrplot(cq$residuals, is.cor = FALSE)


