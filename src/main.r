# Cargar framework
source("./src/configuration.R")
source("./src/plot_functions.R")

fly.rnd <- read.csv(
  "./data/stem15/fly.random.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.uns <- read.csv(
  "./data/stem15/fly.bound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)


unidos <-
  read.table("./data/fly_bound.lst",
    header = F,
    stringsAsFactors = F
  )


nounidos <-
  read.table("./data/fly_non_bound.lst",
    header = F,
    stringsAsFactors = F
  )
reprimidos <-
  read.table("./data/fly_repressed.lst",
    header = F,
    stringsAsFactors = F
  )

uns.total <- unidos %>%
  anti_join(reprimidos) %>%
  nrow()
us.total <- unidos %>%
  inner_join(reprimidos) %>%
  nrow()
nus.total <- reprimidos %>%
  anti_join(unidos) %>%
  nrow()

nounidos %>% nrow()
nuns.total <- nounidos %>%
  anti_join(reprimidos) %>%
  anti_join(unidos) %>%
  nrow()

fly.nuns <- read.csv(
  "./data/stem15/fly.non_bound.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

# First random file
fly.nus <- read.csv(
  "./data/stem15/fly.repressed.local.4.15.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.uns$Tipo <- "Unidos y no silenciados"
fly.nuns$Tipo <- "No unidos y no silenciados"
fly.nus$Tipo <- "No unidos y silenciados"
fly.rnd$Tipo <- "Random"


fly.uns$AccessionID <- as.character(fly.uns$AccessionID)
fly.nuns$AccessionID <- as.character(fly.nuns$AccessionID)
fly.nus$AccessionID <- as.character(fly.nus$AccessionID)

# Tomar intersección entre unidos y silenciados
aux <-
  fly.nus %>%
  select(AccessionID) %>%
  unique() %>%
  inner_join(fly.uns %>% select(AccessionID) %>% unique(),
    by = c("AccessionID")
  )

fly.us <- fly.uns %>% inner_join(aux, by = c("AccessionID"))
fly.us$Tipo <- "Unidos y silenciados"

fly.uns <- fly.uns %>% anti_join(aux, by = c("AccessionID"))
fly.nus <-
  fly.nus %>% anti_join(aux, by = c("AccessionID"))

# Tomar intersección entre no unidos y silenciados

aux <-
  fly.nus %>%
  select(AccessionID) %>%
  unique() %>%
  inner_join(fly.nuns %>% select(AccessionID) %>% unique(),
    by = c("AccessionID")
  )

fly.nuns <-
  fly.nuns %>% anti_join(aux, by = c("AccessionID"))

aux <- NULL

modo <- 2 # 1: Unidos y silenciados vs random / 2: Todos vs todos

if (modo == 1) {
  fly <- rbind(fly.rnd, fly.us)
} else {
  fly <- rbind(fly.uns, fly.us)
  fly <- rbind(fly, fly.nus)
  fly <- rbind(fly, fly.nuns)
  fly <- rbind(fly, fly.rnd)
}
fly$Tipo <- as.factor(fly$Tipo)

fly.uns %>%
  distinct(Gen) %>%
  nrow() # 66
fly.nuns %>%
  distinct(Gen) %>%
  nrow() # 161
fly.nus %>%
  distinct(Gen) %>%
  nrow() # 1170
fly.rnd %>%
  distinct(Gen) %>%
  nrow() # 237
fly.us %>%
  distinct(Gen) %>%
  nrow() # 184

fly.uns %>%
  distinct(AccessionID) %>%
  nrow() # 174
fly.nuns %>%
  distinct(AccessionID) %>%
  nrow() # 393
fly.nus %>%
  distinct(AccessionID) %>%
  nrow() # 3035
fly.rnd %>%
  distinct(AccessionID) %>%
  nrow() # 479
fly.us %>%
  distinct(AccessionID) %>%
  nrow() # 469

total <- c("US" = 636, "UNS" = 263, "NUS" = 4086, "NUNS" = 1028, "RND"=899)
aciertos <- c("US" = 468, "UNS" = 179, "NUS" = 3036, "NUNS" = 395,"RND"=474)

######################################
# prueba de grupos
######################################
ztest <- prop.test(aciertos, n = total, correct = F)
zscore <- c(4,2,3,0,1)
ztest
zprop <- c(0.7358,0.6806,0.7430,0.3842,0.5272)
ztrend <- prop.trend.test(aciertos, n = total,zscore)
ztrend
ztb <- cbind("Puntaje"=zscore,"Proporciones estimadas"=zprop)
ztb <- as.data.frame(ztb)
ztb <- cbind(ztb, Grupo= c("Unidos y silenciados","Unidos y no silenciados","No unidos y silenciados","No unidos ni silenciados","Random"))

zplot <- ztb %>% 
  ggplot(aes(x=Puntaje,y=`Proporciones estimadas`,color=Grupo)) + 
  geom_point() + 
  geom_text(aes(label=round(`Proporciones estimadas`,3)),hjust=0.5, vjust=1.5,size=3)+
  scale_y_continuous(labels = comma)+
  ggtitle("Prueba X-cuadrado de Tendencia en Proporciones", subtitle = "Estadístico X-cuadrado = 398,15 - df = 1 - p-valor < 0,01")

zplot

# USxUNS # p=0.09342 x=2.8144
prop.test(c(468,179), c(636,263), correct = F) 

# USxNUS  p =0.7004 x=0.14806
prop.test(c(468,3036),c(636,4086), correct = F) # 

# USxNUNS  p < 0.01 x=194.57
prop.test(c(468,395),c(636,1028), correct = F) 

# USxRND  # p < 0.01 x= 68,365
prop.test(c(468,474),c(636,899), correct = F) 

# UNSxNUS # p = 0.02544 x=4.994
prop.test(c(179,3036),c(263,4086), correct = F) 

# UNSxNUNS # p < 0.01 x= 74.491
prop.test(c(179,395),c(263,1028), correct = F) 

# UNSxRND # p < 0.01 x = 19.44
prop.test(c(179,474),c(263,899), correct = F) 

# NUSxNUNS # p < 0.01 x = 478,86
prop.test(c(3036,395),c(4086,1028), correct = F) 

# NUSxRND # p < 0.01 x = 164.67
prop.test(c(3036,474),c(4086,899), correct = F) 

# NUNSxRND # p < 0.01 x = 39.616
prop.test(c(395,474),c(1028,899), correct = F) 

######################################
# prueba de grupos
######################################

nrow(fly.uns) # 292
nrow(fly.nuns) # 586
nrow(fly.nus) # 6326
nrow(fly.rnd) # 758
nrow(fly.us) # 1083

summary(fly.uns$LoopPattern) # 292
summary(fly.nuns$LoopPattern) # 586
summary(fly.nus$LoopPattern) # 6326
summary(fly.rnd$LoopPattern) # 758
summary(fly.us$LoopPattern) # 1083

fly$Id <- 1:nrow(fly)
fly$N10 <- as.factor(fly$N9)
fly$N9 <- as.factor(fly$N8)
fly$N8 <- as.factor(fly$N7)
fly$N7 <- as.factor(fly$N6)
fly$N6 <- as.factor(fly$N5)
fly$N5 <- as.factor(fly$N4)
fly$N4 <- as.factor(fly$N3)
fly$N3 <- as.factor(fly$N2)
fly$N2 <- as.factor(fly$N1)
fly$N1 <- as.factor(fly$N0)
fly$N0 <- NULL

levels(fly$N1) <- c("A", "C", "G", "U")
levels(fly$N2) <- c("A", "C", "G", "U")
levels(fly$N3) <- c("A", "C", "G", "U")
levels(fly$N4) <- c("A", "C", "G", "U")
levels(fly$N5) <- c("A", "C", "G", "U", NA)
levels(fly$N6) <- c("A", "C", "G", "U", NA)
levels(fly$N7) <- c("A", "C", "G", "U", NA)
levels(fly$N8) <- c("A", "C", "G", "U", NA)
levels(fly$N9) <- c("A", "C", "G", "U", NA)
levels(fly$N10) <- c("A", "C", "G", "U", NA)


# revisar las proporciones informadas en cada set
featuresOfThePlot(fly) # ok 1
geneHitsByPattern(fly) # ok 1
locationAndPairs(fly) # ok 1
plotLoopPatterns(fly) # ok 1
plotNBases(fly) # ok 1
wpbPlots(fly) # ok 1
basePercentsPlot(fly) # no se
pairsMFEPlot(fly) # boo
purinePlot(fly)
relativePositionPlot(fly)
locationPlot(fly)

# ANOVA
tb <- table(fly$Tipo, fly$LoopPattern)
cq <- chisq.test(tb)


# HEATMAP
library(reshape2)
un <- filter(fly, Tipo == "Unidos y no silenciados")
un <- scale(table(un$TerminalPair, un$Location))
un <- melt(un)
un$Tipo <- "Unidos y no silenciados"
us <- filter(fly, Tipo == "Unidos y silenciados")
us <- scale(table(us$TerminalPair, us$Location))
us <- melt(us)
us$Tipo <- "Unidos y silenciados"
nn <- filter(fly, Tipo == "No unidos y no silenciados")
nn <- scale(table(nn$TerminalPair, nn$Location))
nn <- melt(nn)
nn$Tipo <- "No unidos y no silenciados"
ns <- filter(fly, Tipo == "No unidos y silenciados")
ns <- scale(table(ns$TerminalPair, ns$Location))
ns <- melt(ns)
ns$Tipo <- "No unidos y silenciados"
ra <- filter(fly, Tipo == "Random")
ra <- scale(table(ra$TerminalPair, ra$Location))
ra <- melt(ra)
ra$Tipo <- "Random"
data <- rbind(un, us, nn, ns, ra)

pal <- wes_palette("Zissou1", 100, type = "continuous")

ggplot(data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "black", high =
      "blue"
  ) +
  xlab("Cierre terminal") +
  ylab("Sitio") +
  facet_wrap(~Tipo) +
  guides(
    fill =
      F
  ) +
  ggtitle("Cierres terminales por cada sitio del transcrito")


####################################################
## ESTADISTICOS
####################################################

myTable <- fly %>%
  select(Tipo, N.2, N.1, N2, N5, N6, N7, N8) %>%
  gather(metric, value, -Tipo) %>%
  filter(value %in% c("A", "G", "C", "U")) %>%
  mutate(metric = str_replace_all(metric, "N.1", "N(-1)")) %>%
  mutate(metric = str_replace_all(metric, "N.2", "N(-2)"))

# install.packages("Cairo")
# install.packages("cairoDevice")
# install.packages("vcd")
library(corrplot)
library(Cairo)
library(vcd)
library(cairoDevice)

n_2 <- myTable %>% filter(metric == "N(-2)")
n_2 <- table(n_2$Tipo, n_2$value)
cq <-
  chisq.test(n_2) # X-squared = 88.471, df = 12, p-value = 9.756e-14

corrplot(cq$residuals, is.cor = FALSE, tl.srt = 0, tl.offset = 0.8)

corrplot(cq$residuals, is.corr = FALSE)
corrplot(cq$stdres, is.corr = FALSE)

vcd::assocstats(n_2)
contrib <- 100 * cq$residuals^2 / cq$statistic
round(contrib, 3)
corrplot(contrib, is.cor = FALSE)

n_1 <- myTable %>% filter(metric == "N(-1)")
n_1 <- table(n_1$Tipo, n_1$value)
cq1 <-
  chisq.test(n_1) # X-squared = 80.001, df = 12, p-value = 4.126e-12
corrplot(cq1$stdres,
  is.cor = FALSE, method = "square", tl.srt = 0,
  tl.offset = 0.8
)
vcd::assocstats(n_1)
contrib1 <- 100 * cq1$residuals^2 / cq1$statistic
round(contrib, 3)
corrplot(contrib1, is.cor = FALSE)

n2 <- myTable %>% filter(metric == "N2")
n2 <- table(n2$Tipo, n2$value)
cq2 <-
  chisq.test(n2) # X-squared = 164.23, df = 12, p-value < 2.2e-16
corrplot(cq2$residuals, is.cor = FALSE, method = "square", tl.srt = 0, tl.offset = 0.8)
contrib2 <- 100 * cq2$residuals^2 / cq2$statistic
corrplot(contrib2, is.cor = FALSE)

n5 <- myTable %>% filter(metric == "N5")
n5 <- table(n5$Tipo, n5$value)
cq <-
  chisq.test(n5) # X-squared = 51.813, df = 12, p-value = 6.692e-07
corrplot(cq$residuals, is.cor = FALSE, tl.srt = 0, tl.offset = 0.8)

n6 <- myTable %>% filter(metric == "N6")
n6 <- table(n6$Tipo, n6$value)
cq <-
  chisq.test(n6) # X-squared = 42.266, df = 12, p-value = 3.004e-05
corrplot(cq$residuals, is.cor = FALSE, tl.srt = 0, tl.offset = 0.8)

n7 <- myTable %>% filter(metric == "N7")
n7 <- table(n7$Tipo, n7$value)
cq <-
  chisq.test(n7) # X-squared = 48.854, df = 12, p-value = 2.218e-06
corrplot(cq$residuals, is.cor = FALSE)

n8 <- myTable %>% filter(metric == "N8")
n8 <- table(n8$Tipo, n8$value)
cq <-
  chisq.test(n8) # X-squared = 39.747, df = 12, p-value = 7.92e-05 *may be incorrect
corrplot(cq$residuals, is.cor = FALSE, tl.srt = 0, tl.offset = 0.8)


###############################################################################
# ANALISIS DE COMPONENTES PRINCIPALES
###############################################################################
# install.packages("Gifi")
library(Gifi)
prueba <-
  princals(
    select(
      fly,
      Location,
      LoopPattern,
      TerminalPair,
      N.2,
      N.1,
      N2,
      N5,
      Pairments,
      WooblePairs,
      Bulges,
      InternalLoops,
      RnaFoldMFE,
      CG_PercentPairs,
      AU_PercentPairs
    )
  )

plot(prueba, plot.type = "transplot")
plot(prueba, plot.type = "loadplot")
plot(prueba, plot.type = "biplot")


###############################################################################
# k-prototipos
###############################################################################
install.packages("clustMixType")
library(clustMixType)

dt <-
  select(
    fly,
    Location,
    LoopPattern,
    TerminalPair,
    N.2,
    N.1,
    N2,
    N5,
    Pairments,
    WooblePairs,
    Bulges,
    InternalLoops,
    RnaFoldMFE,
    CG_PercentPairs,
    AU_PercentPairs
  )
clst <- kproto(
  dt,
  k = 4,
  iter.max = 1000,
  nstart = 1,
  na.rm = TRUE
)

clprofiles(clst$cluster, dt, col = wes_palette("Royal1", 4, type = "continuous")) # figure 1

Es <- numeric(10)
for (i in 1:10) {
  kpres <- kproto(dt, k = i, nstart = 5)
  Es[i] <- kpres$tot.withinss
}
plot(
  1:10,
  Es,
  type = "b",
  ylab = "Objective Function",
  xlab = "# Clusters",
  main = "Scree Plot"
) # figure 2


###############################################################################
# SRE por cada GEN
###############################################################################


glimpse(fly.uns)
fly$Gen <- as.character(fly$Gen)

df <- fly %>%
  filter(Tipo != "Random") %>%
  distinct(Gen, Loop, TerminalPair, Pairments, Tipo) %>%
  group_by(Gen) %>%
  select(Gen, Tipo) %>%
  count()
View(df)

df %>%
  filter(Tipo != "Random") %>%
  ggplot(aes(x = Tipo, y = freq)) +
  geom_boxplot() +
  ylab("SREs únicos") +
  xlab("Grupo") +
  ggtitle("Cantidad de SRE por cada gen")
View(df %>% filter(Tipo == "Random"))
df %>%
  filter(freq > 3) %>%
  arrange(desc(freq))

summary(df)

df$Gen <- as.factor(df$Gen)
View(df)
df %>%
  filter(freq > 4) %>%
  arrange(desc(freq)) %>%
  ggplot(aes(
    x = freq,
    y = reorder(Gen, freq),
    color = Tipo
  )) +
  geom_point(size = 4, alpha = 0.7) +
  geom_segment(aes(
    x = 0,
    y = Gen,
    xend = freq,
    yend = Gen
  ), color = "grey50") +
  theme(axis.text.y = element_text(size = 7)) +
  scale_x_continuous(
    "Cantidad de SRE detectados por la aplicación",
    limits = c(0, 14),
    breaks = seq(1, 14, 2)
  ) +
  ylab("Gen") +
  ggtitle("Genes con mayor cantidad de SREs") +
  scale_color_manual(
    name = "Grupo",
    values = c("green", "red", "yellow", "orange")
  )

#######################

library(corrplot)
library(Cairo)
library(vcd)
library(cairoDevice)



# Si la prueba de Chi-cuadrado es significativa necesitamos mirar a los residuos.
# Si el valor de los residuos estandarizados es menor que -2 significa que la celda contiene menos
# observaciones de lo esperado (la variable sería independiente). Si el valor de los residuos estandarizados
# es mayor que 2 significa que la celda contiene mas observaciones de lo esperado.
# Residuos positivos (en azul) en celdas especifican una atracción (asociación positiva)
# Residuos negativos (en rojo) implican una repulsión (asociación negativa) entre las variables de la fila y columna correspondiente.

# If the chi squared test is significant we need to take a look at residuals.
# If the value of standardized residual is lower than -2 it means that the cell contains fewer
# observations that it was expected (the case of variables independence). If the value of standardized
# residual is higher than 2 it means that the cell contains more observations that it was expected.
# Positive residuals are in blue. Positive values in cells specify an attraction (positive association)
# between the corresponding row and column variables.
# Negative residuals are in red. This implies a repulsion (negative association)between the corresponding
# row and column variables.
# We can plot the results using assocplot() function

myTable <- fly %>%
  select(Tipo, N.2, N.1, N2, N5, N6, N7, N8) %>%
  gather(metric, value, -Tipo) %>%
  filter(value %in% c("A", "G", "C", "U")) %>%
  mutate(metric = str_replace_all(metric, "N.1", "N(-1)")) %>%
  mutate(metric = str_replace_all(metric, "N.2", "N(-2)"))

n_2 <- myTable %>% filter(metric == "N(-2)")
n_2 <- table(n_2$Tipo, n_2$value)
cq <-
  chisq.test(n_2) # X-squared = 88.471, df = 12, p-value = 9.756e-14

corrplot(cq$stdres, is.cor = FALSE, tl.srt = 0, tl.offset = 0.8, method = "square")

n_1 <- myTable %>% filter(metric == "N(-1)")
n_1 <- table(n_1$Tipo, n_1$value)
cq1 <-
  chisq.test(n_1) # X-squared = 80.001, df = 12, p-value = 4.126e-12
corrplot(cq1$stdres,
  is.cor = FALSE, method = "square", tl.srt = 0,
  tl.offset = 0.8
)
vcd::assocstats(n_1)
contrib1 <- 100 * cq1$residuals^2 / cq1$statistic
round(contrib, 3)
corrplot(contrib1, is.cor = FALSE)

n2 <- myTable %>% filter(metric == "N2")
n2 <- table(n2$Tipo, n2$value)
cq2 <-
  chisq.test(n2) # X-squared = 164.23, df = 12, p-value < 2.2e-16
corrplot(cq2$stdres, is.cor = FALSE, method = "square", tl.srt = 0, tl.offset = 0.8)
contrib2 <- 100 * cq2$residuals^2 / cq2$statistic
corrplot(contrib2, is.cor = FALSE)

n5 <- myTable %>% filter(metric == "N5")
n5 <- table(n5$Tipo, n5$value)
cq <-
  chisq.test(n5) # X-squared = 51.813, df = 12, p-value = 6.692e-07
corrplot(cq$stdres, is.cor = FALSE, method = "square", tl.srt = 0, tl.offset = 0.8)

n6 <- myTable %>% filter(metric == "N6")
n6 <- table(n6$Tipo, n6$value)
cq <-
  chisq.test(n6) # X-squared = 42.266, df = 12, p-value = 3.004e-05
corrplot(cq$stdres, is.cor = FALSE, method = "square", tl.srt = 0, tl.offset = 0.8)

n7 <- myTable %>% filter(metric == "N7")
n7 <- table(n7$Tipo, n7$value)
cq <-
  chisq.test(n7) # X-squared = 48.854, df = 12, p-value = 2.218e-06
corrplot(cq$stdres, is.cor = FALSE, method = "square", tl.srt = 0, tl.offset = 0.8)

n8 <- myTable %>% filter(metric == "N8")
n8 <- table(n8$Tipo, n8$value)
cq <-
  chisq.test(n8) # X-squared = 39.747, df = 12, p-value = 7.92e-05 *may be incorrect
corrplot(cq$stdres, is.cor = FALSE, method = "square")


# pares terminales
myTable <- fly %>%
  select(Tipo, TerminalPair) %>%
  gather(metric, value, -Tipo)


pairs <- myTable %>% filter(metric == "TerminalPair")
pairs <- table(pairs$Tipo, pairs$value)
vcd::assocstats(pairs)
cq <-
  chisq.test(pairs) # X-squared = 88.471, df = 12, p-value = 9.756e-14

corrplot(cq$stdres, is.cor = FALSE, tl.srt = 0, tl.offset = 1, method = "square")

# FISHER
locs <- table(fly$Tipo, fly$Location)
locs.prop <- prop.table(locs, 1)
locs.cq <- chisq.test(locs)
corrplot::corrplot(locs.cq$stdres, is.corr = F, method = "square")
vcdd <- vcd::assocstats(locs)
