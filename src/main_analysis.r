# Cargar framework
source("./src/configuration.R")
source("./src/plot_functions.R")

###############################################################################
# Carga y limpieza de datos
###############################################################################

# Cargar datos
fly.bound <- read.csv(
  "./data/stem15/fly_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.non_bound <- read.csv(
  "./data/stem15/fly_non_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

# First random file
fly.random <- read.csv(
  "./data/stem15/fly_bound_rnd.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)
# Second random file
# fly.random <- rbind(fly.random,
#                     read.csv("./data/stem15/fly_bound_rnd1.csv",
#                               sep =";",
#                               dec =",", stringsAsFactors = TRUE))


# Eliminar variables innecesarias
fly.non_bound$Tipo <- "No unidos"
#fly.non_bound <- formatEnsembl(fly.non_bound, FALSE)
fly.bound$Tipo <- "Unidos"
#fly.bound <- formatEnsembl(fly.bound, FALSE)
fly.non_bound$Gen <- as.character(fly.non_bound$Gen)
fly.random$Column3 <- as.character(fly.random$Column3)
fly.bound$Gen <- as.character(fly.bound$Gen)


# add random
fly.random$Tipo <- "Random"
#fly.random <- formatEnsembl(fly.random, FALSE)

# check intersections
fly.bound %>% inner_join(fly.non_bound, by = c("Gen"))

# igualar datasets
fly.bound$GeneSynonym <- NULL
fly.bound$Note <- NULL
fly.bound$AccessionID <- NULL
fly.non_bound$GeneSynonym <- NULL
fly.non_bound$Note <- NULL
fly.non_bound$AccessionID <- NULL

fly_gb <- rbind(fly.bound, fly.non_bound)

fly.random$Gen <- fly.random$Column3
fly.random$Column1 <- NULL
fly.random$Column2 <- NULL
fly.random$Column3 <- NULL
fly.random$Column4 <- NULL
fly.random$Column5 <- NULL
fly.random$Column6 <- NULL

fly <-
  rbind(select(fly_gb, -CDS_Start, -CDS_End, -Location), fly.random)

# Determinando niveles
fly$Tipo <-
  factor(fly$Tipo, levels = c("Unidos", "No unidos", "Random"))

glimpse(fly.bound)
glimpse(fly.non_bound)
glimpse(fly.random)
glimpse(fly_gb)


###############################################################################
# Análisis de las variables
###############################################################################

# Hay correlaciones??
cfly <- fly.bound %>%
  select(
    -Tipo,
    -Gen,
    -Location,
    -Loop,
    -SequenceLength,
    -TerminalPair,
    -AdditionalSeqMatches,
    -EndsAt,
    -CDS_Start,
    -CDS_End,
    -StartsAt,
    -StemLoopSequence,
    -GU_PercentPairs,
    -AU_PercentPairs,
    -A_PercentSequence,
    -Sense,
    -AdditionalSeqPositions,
    -PredictedStructure,
    -ViennaBracketStr,
    -LoopPattern
  )

cfly %>% getCorPearson()
cfly %>% getCorTau()

fly$id <- fly$Gen #paste(fly$Gen, fly$N.2,fly$N.1, fly$Loop)


#######################################
# Plot 1 - Features of the stem-loop
featuresOfThePlot(fly)

# Plot  2 - Gene hits by loop pattern
geneHitsByPattern(fly)

# Plot  3 - Chromosome, patterns and pairs
chromoPattAndPairs(fly)

# Plot  4 - N Bases
plotNBases(fly)

# Plot  5 - Pairments, WooblePairs, Bulges InternalLoops
wpbPlots(fly)

# Plot  6 - Base Percent sequences
basePercentsPlot(fly)

# Plot  7 - Percent bases in pairs, MFE
pairsMFEPlot(fly)

# Plot  8 - Purine percent in pairs
purinePlot(fly)

# Plot  9 - Relative position
relativePositionPlot(fly)


#######################################
# Plot  x - Sin desarrollar
#######################################

# # "GenSymbol"
# a <- fly %>% distinct(Tipo,GenSymbol, id) %>%
#   filter(Tipo=="Unidos") %>%
#   group_by(GenSymbol) %>%
#   dplyr::summarise(Unidos=n())
#
# b <- fly %>% distinct(Tipo,GenSymbol, id) %>%
#   filter(Tipo=="Random") %>%
#   group_by(GenSymbol) %>%
#   dplyr::summarise(Random=n())
#
# nrow(a)
# nrow(b)
#
# a %>% inner_join(b) %>% filter(Unidos>2) %>%
#   mutate(y1=Unidos/sum(Unidos),
#          y2=Random/sum(Random)) %>%
#   arrange(Unidos) %>%
#   head(100) %>%
#   ggplot() +
#   geom_segment( aes(x=GenSymbol, xend=GenSymbol,
#                     y=y1, yend=y2), color="grey") +
#   geom_point( aes(x=GenSymbol, y=y1), color=rgb(0.2,0.7,0.1,0.5), size=2.5 ) +
#   geom_point( aes(x=GenSymbol, y=y2), color=rgb(0.7,0.2,0.1,0.5), size=1.5 ) +
#   coord_flip() + ylab("Frecuencia") +
#   theme_light() +
#   theme(
#     legend.position = "none",
#     panel.border = element_blank(),
#     axis.text.y = element_text(size=8)
#   ) + xlab("")


###############################################################################
# Validación estadística
###############################################################################


mfly <- fly %>% filter(Tipo != "No unidos")
mfly$Tipo <- factor(mfly$Tipo)
tb <- table(mfly$Tipo, mfly$LoopPattern)
# balloonplot(t(tb), main ="patterns", xlab ="", ylab="",
#             label = FALSE, show.margins = FALSE)
#
# mosaicplot(tb, shade = TRUE, las=2,
#            main = "housetasks")
print(tb)
chisq <- chisq.test(tb)
chisq$p.value < 0.05
chisq
