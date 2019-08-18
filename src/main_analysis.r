# Cargar framework
source("./src/configuration.R")
library(dplyr)

###############################################################################
# Carga y limpieza de datos
###############################################################################

# Cargar datos
fly.bound <- read.csv("./data/stem15/fly_bound.csv",
                      sep =";", 
                      dec =",", stringsAsFactors = TRUE)

fly.non_bound <- read.csv("./data/stem15/fly_non_bound.csv",
                          sep =";", 
                          dec =",", stringsAsFactors = TRUE)

# First random file
fly.random <- read.csv("./data/stem15/fly_bound_rnd.csv",
                          sep =";",
                          dec =",", stringsAsFactors = TRUE)
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
fly.bound %>% inner_join(fly.non_bound, by=c("Gen"))

glimpse(fly.bound)
glimpse(fly.non_bound)
glimpse(fly.random)

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

fly <- rbind(select(fly_gb, -CDS_Start, -CDS_End, -Location), fly.random)

# Determinando niveles
fly$Tipo <- factor(fly$Tipo, levels=c("Unidos", "No unidos", "Random"))


###############################################################################
# Análisis de las variables
###############################################################################

if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridEXtra)
}

# Hay correlaciones??
cfly <- fly.bound %>% 
  select(-Tipo,-Gen, -Location, -Loop, -SequenceLength,
         -TerminalPair, -AdditionalSeqMatches, -EndsAt,-CDS_Start,-CDS_End,-StartsAt,-StemLoopSequence,
         -GU_PercentPairs,-AU_PercentPairs, -A_PercentSequence, -Sense, -AdditionalSeqPositions,
         -PredictedStructure,-ViennaBracketStr,-LoopPattern) 

cfly %>% getCorPearson()
cfly %>% getCorTau()

fly$id <- fly$Gen #paste(fly$Gen, fly$N.2,fly$N.1, fly$Loop) 


#######################################
# Plot 1 - Features of the stem-loop
#######################################

# Plot by pairs
pairs <- fly %>%  filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, Pairments, id) %>% 
  group_by(Tipo,LoopPattern, Pairments) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=Pairments, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed") + 
  xlab("Longitud del stem") +
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  facet_grid(~LoopPattern) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

# Plot by wooble
wooble <- fly %>% filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, WooblePairs, id) %>% 
  group_by(Tipo,LoopPattern, WooblePairs) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=WooblePairs, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed", show.legend = F) + 
  xlab("Apareamientos GU") +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  facet_grid(~LoopPattern) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

# Plot by bulges
bulges <- fly %>% filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, Bulges, id) %>% 
  group_by(Tipo,LoopPattern, Bulges) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=Bulges, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed", show.legend = F) + 
  xlab("Bulges") +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  facet_grid(~LoopPattern)+
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

# Plot by internal loops
internals <- fly %>% filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, InternalLoops, id) %>% 
  group_by(Tipo,LoopPattern, InternalLoops) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=InternalLoops, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed", show.legend = F) + 
  xlab("Loops internos") +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  facet_grid(~LoopPattern)+
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

# Arrange grid with all plots
grid.arrange(pairs, wooble, bulges, internals, ncol=1,
             top = text_grob("% Genes acertados por patrón del loop", 
                             face="bold",size = 13.5), bottom = "")


#######################################
# Plot  2 - Gene hits by loop pattern
#######################################

n1 <- fly %>%  filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, N.1, id) %>% 
  group_by(Tipo,LoopPattern, N.1) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=N.1, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed") + 
  xlab("Base predecesora (N-1)") +
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  facet_grid(~LoopPattern) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

n2 <- fly %>% filter(str_length(LoopPattern) <= 8) %>%
  distinct(Tipo,LoopPattern, N2, id) %>% 
  group_by(Tipo,LoopPattern, N2) %>% 
  dplyr::summarise(n=n()) %>%
  mutate(freq=n/sum(n)) %>%
  ggplot(aes(x=N2, y=freq*100, group=Tipo,color=Tipo )) +
  geom_line(linetype="dashed", show.legend = F) + 
  xlab("Base N2") +
  geom_point(aes(shape=Tipo), show.legend = F) + 
  ylab("%") + labs(color="Listado") +
  facet_grid(~LoopPattern)+
  scale_y_continuous(breaks=seq(0, 60, by=10), limits=c(0,60)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="floralwhite"))

grid.arrange(n1, n2, ncol=1,
             top = text_grob("% Genes acertados por patrón del loop", 
                             face="bold",size = 13.5),
             bottom = "")


#######################################
# Plot  3 - Chromosome, patterns and pairs
#######################################

if(!require(wesanderson)){
  install.packages("wesanderson")
  library(wesanderson)
  }

# Chromosome
chromo <- fly %>% filter(Tipo!="Random") %>%
  ggplot(aes(x=Chromosome, fill=Tipo)) +
  geom_bar(color = "darkgray", position="fill", width = 0.6) + theme_gray() +
  ggtitle("Distribución de los dos conjuntos en cromosomas",
          subtitle = "") + labs(fill="Grupo") +
  scale_fill_manual(values=c('blue2','brown3')) + xlab("Cromosoma") + ylab("Recuento") #480x480px


# LoopPattern 
#source("https://install-github.me/larmarange/JLutils")

patt <- fly %>%
  ggplot(aes(fill=LoopPattern, x=1)) +
  geom_bar(position="fill", color="darkgrey", ) + coord_polar(theta="y") +
  geom_text(stat = "fill_labels", fontface="bold",
            position = position_dodge(width = .8), size=3.5, 
            check_overlap = TRUE, ) + labs(fill="Patrón del loop") +
  facet_wrap(~ Tipo) + ylab("") +
  ggtitle("Secuencias consenso buscadas en los bucles", 
          subtitle="Sobre todos los stem-loops hallados en los transcriptos.") +
  theme(axis.text.x=element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold")) + 
  ylab("") + xlab("") +
  scale_fill_brewer(palette="Spectral") 

# TerminalPair
pairs <- fly %>% ggplot(aes(Tipo,fill=TerminalPair)) + 
  geom_bar(stat="count") + labs(fill="Par de cierre") +
  ggtitle("Apareamientos de cierre del loop", 
          subtitle="Sobre todos los stem-loops hallados en los transcriptos.") +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="darkgrey", fill=NA)) + xlab("Grupo") +
  stat_stack_labels() +
  scale_fill_manual(values = wes_palette(6, name = "Darjeeling1", type = "continuous"))

grid.arrange(pairs, patt)


#######################################
# Plot  4 - N Bases
#######################################

fly %>% 
  select (Tipo,N.2,N.1,N2,N5,N6,N7,N8) %>% 
  gather(metric,value,-Tipo ) %>%
  filter(value %in% c("A","G","C","U")) %>%
  mutate(metric = str_replace_all(metric, "N.1", "N(-1)")) %>%
  mutate(metric = str_replace_all(metric, "N.2", "N(-2)")) %>%
  ggplot(aes(fill=value, x=1)) +
  geom_bar(position="fill", color="darkgrey", ) + coord_polar(theta="y") +
  labs(fill="Base") +
  ggtitle("Frecuencia de bases en cada posición variable (N) del loop",
          subtitle="Sobre todos los stem-loops hallados en los transcriptos.") +
  geom_text(stat = "fill_labels", #fontface="bold",
            position = position_dodge(width = .8), size=3.2, 
            check_overlap = TRUE, color="red") +
  facet_grid(Tipo ~ metric) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill="beige"),
        panel.border = element_rect(colour="darkgrey", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.background = element_rect(fill="seashell")) +
  scale_fill_manual(values=wes_palette("IsleofDogs2"))
  

#######################################
# Plot  5 - Pairments, WooblePairs, Bulges InternalLoops
#######################################

fly %>%  
  select(Tipo,WooblePairs,Bulges,InternalLoops) %>%
  gather(metric, value, -Tipo) %>%
  mutate(metric=str_replace_all(metric, "WooblePairs", "Pares GU")) %>%
  mutate(metric=str_replace_all(metric, "InternalLoops", "Loops internos")) %>%
  ggplot(aes(x=value,fill=metric)) +
  geom_histogram(show.legend = F, bins = 8, color="black", lwd=0.8) +
  labs(fill="") +
  theme(panel.background = element_rect(fill="azure"),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        strip.background = element_blank(),
        strip.text = element_text(face="bold"),
        plot.background = element_rect(fill="whitesmoke")) +
  ggtitle("Propiedades de los stem-loops", subtitle="") + 
  xlab("Cantidad") + ylab("Stem-loops") +
  facet_grid(Tipo ~ metric) +
  scale_fill_manual(values=wes_palette("Darjeeling2"))


#######################################
# Plot  6 - Base Percent sequences
#######################################

fly %>% filter(Tipo!="Random") %>%
  select(Tipo, A_PercentSequence, C_PercentSequence, 
         G_PercentSequence, U_PercentSequence) %>% 
  gather(metric, value, -Tipo) %>%
  mutate(metric=paste ("", str_replace(metric, "_PercentSequence", "")), value = value * 100) %>%
  ggplot(aes(x=value, fill=Tipo)) + ylab("Densidad") +
  geom_histogram(aes(y=..density..), position="identity", bins=40, alpha=0.5, color="black") + 
  xlab("Ratio de cada nucleótido en la secuencia entera") +
  facet_wrap(~ metric,ncol=2, nrow=2 ) + theme_pubr() +
  ggtitle("Composición de las secuencias cDNA", 
          subtitle = "Sobre el total de bases en transcriptos") +
  scale_fill_manual(values=wes_palette("Cavalcanti1"))


#######################################
# Plot  7 - Percent bases in pairs, MFE
#######################################

# Percent Pairs in Stem 
pairs2 <- fly %>% filter(Tipo!="Random") %>%
  select(Tipo, AU_PercentPairs, CG_PercentPairs, GU_PercentPairs) %>% 
  gather(metric, value, -Tipo) %>%
  mutate(metric=paste ("", str_replace(metric, "_PercentPairs", "")), value = value * 100) %>%
  ggplot(aes(y=value, x=Tipo)) + ylab("Densidad") +
  geom_boxplot(aes(colour=metric),
               outlier.colour = "red", outlier.shape = 1,
               lwd = 1.3) + 
  labs(colour="Apareamiento") +
  xlab("Ratio de cada apareamiento en los stem-loop") +
  ggtitle("Composición de apareamientos en los stem-loops", 
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  theme(panel.border = element_rect(colour="darkgrey", fill=NA)) +
  scale_colour_manual(values=wes_palette("FantasticFox1")) 

# Minimum free energy
mfe <- fly %>% 
  select(Tipo, RnaFoldMFE) %>% 
  ggplot(aes(x=RnaFoldMFE, group=Tipo)) + ylab("Densidad") +
  geom_density(alpha=1, size=0.7) +
  facet_wrap(~ Tipo) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold")) +
  scale_x_continuous(breaks = round(seq(min(fly$RnaFoldMFE), 
                                        max(fly$RnaFoldMFE), by = 4),1))  +
  scale_y_continuous(breaks=seq(0, 0.2, by=0.02)) +
  xlab("Mínima energía libre (kcal/mol)") +
  ggtitle("Energía libre en los stem-loops", 
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  geom_vline(data=filter(fly,Tipo=="Unidos"), aes(xintercept=mean(RnaFoldMFE), group=Tipo),
             color="red", linetype="dashed", show.legend = F, size=1) +
  geom_vline(data=filter(fly,Tipo=="No unidos"), aes(xintercept=mean(RnaFoldMFE), group=Tipo),
             color="blue", linetype="dashed", show.legend = F, size=1) +
  geom_vline(data=filter(fly,Tipo=="Random"), aes(xintercept=mean(RnaFoldMFE), group=Tipo),
             color="darkgreen", linetype="dashed", show.legend = F, size=1) 

grid.arrange(pairs2, mfe)

#######################################
# Plot  8 - Purine percent in pairs
#######################################

# PurinePercentPairs
fly %>%
  select(Tipo, PurinePercentPairs) %>%
  ggplot(aes(x=PurinePercentPairs, group=Tipo)) + ylab("Densidad") +
  geom_density(alpha=1, size=1) +
  facet_wrap(~ Tipo) +
  xlab("Ratio de A-G en apareamientos de los stems") +
  theme_pubr() +
  theme(strip.background = element_rect(colour="white", fill="lightgrey"),
        strip.text = element_text(face="bold")) +
  ggtitle("Contenido de bases A-G (purinas) en los stem-loops",
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  geom_vline(data=fly, aes(xintercept=mean(PurinePercentPairs), color="blue"),
             linetype="dashed", show.legend = F, size=1)


#######################################
# Plot  9 - Relative position
#######################################

# Cambiar por posicion UTR/CDS
fly %>% ggplot(aes(color=Tipo, x=RelativePosition, y=(-1)*RnaFoldMFE)) +
  geom_smooth(method="loess") +
  ggtitle("Variación de la estabilidad con respecto a la posición relativa en la secuencia") +
  theme_pubr()


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

if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
}


# if(!require(graphics)){
#   install.packages("graphics")
#   library(graphics)
# }

fly <- fly %>% filter(Tipo != "No unidos")
fly$Tipo <- factor(fly$Tipo)
tb <- table(fly$LoopPattern)
balloonplot(t(tb), main ="patterns", xlab ="", ylab="",
            label = FALSE, show.margins = FALSE)

mosaicplot(tb, shade = TRUE, las=2,
           main = "housetasks")
chisq <- chisq.test(tb)
chisq$p.value < 0.05


