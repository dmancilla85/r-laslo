
setwd(paste(getwd(), "/data/mitocondriales", sep=""))




# ****************************************************************

#fly_bound <- fly_bound %>% distinct(Column3, Loop, TerminalPair, .keep_all = TRUE)


# ****************************************************************

file <- "yeast_she"
yeast <- read.csv(paste(file, ".csv", sep=""), header=TRUE, sep=";", dec=",", 
                  stringsAsFactors = FALSE)
yeast$SerieID <- 1
yeast$Serie = "Secuencia cDNA"

print("Cargando randoms...")
for(i in 1:150){
  df <- read.csv(paste(file,"_rnd", i, ".csv", sep=""), header=TRUE, sep=";", 
                 dec=",", stringsAsFactors = TRUE)
  df$Serie <- "Secuencia permutada" #paste("Control", i)
  df$SerieID <- 2
  yeast <- rbind(yeast, df)
}


fly_bound <- filter(fly_bound, RnaFoldMFE < 0)
fly_repressed <- filter(fly_repressed, RnaFoldMFE < 0)
yeast <- filter(yeast, RnaFoldMFE < 0)
mouse <- filter(mouse, RnaFoldMFE < 0)

head(yeast)
fly_bound$GenID <- fly_bound$Column1
fly_bound$Transcripto <- fly_bound$Column2
fly_bound$Gen <- fly_bound$Column3
mouse$GenID <- mouse$Column1
mouse$Transcripto <- mouse$Column2
mouse$Gen <- mouse$Column5
fly_repressed$GenID <- fly_repressed$Column1
fly_repressed$Transcripto <- fly_repressed$Column2
fly_repressed$Gen <- fly_repressed$Column3
yeast$GenID <- yeast$Column1
yeast$Transcripto <- yeast$Column2
yeast$Gen <- yeast$Column4
fly_bound$Specie <- "D. melanogaster"
fly_repressed$Specie <- "D. melanogaster"
mouse$Specie <- "M. musculus"
yeast$Specie <- "S. cerevisieae"

all <- rbind(fly_bound, yeast, mouse)

all$Column1 <- NULL
all$Column2 <- NULL
all$Column3 <- NULL
all$Column4 <- NULL
all$Column5 <- NULL
all$Column6 <- NULL
summary(all)


fileConn<-file("fruitfly_bound.txt")
asd <- fly_bound[ fly_bound$SerieID ==1,]
writeLines(paste(paste(">", paste(asd$Gen,asd$StartsAt,sep=".")),
                 paste(asd$N.2,
                       asd$N.1,
                       asd$Loop,
                 sep=""),sep='\n'),
                 fileConn)
close(fileConn)


fileConn<-file("fruitfly_bound_rnd.txt")
asd <- fly_bound[ fly_bound$SerieID ==2,]
asd <- asd[1:4000,]
asd$Gen <- 1:4000
writeLines(paste(paste(">", paste(asd$Gen,asd$StartsAt,sep=".")),
                 paste(asd$N.2,
                       asd$N.1,
                       asd$Loop,
                       sep=""),sep='\n'),
           fileConn)
close(fileConn)

summary(asd$N.2)

############################
fly_repressed$Specie <- "Reprimidos o silenciados por Smaug"
fly_bound$Specie <- "Unidos a Smaug"

fly_repressed <- anti_join(fly_repressed, fly_bound, by=c("Gen","Transcripto"))

fly <- rbind(fly_bound, fly_repressed)

fly.bases <- fly %>% dplyr::select(Specie, Serie, "N(-1)"=N.1, "N(-2)"=N.2, 
                            "N2"=N2, "N5"=N5, "N6"=N6,"N7"=N7,"N8"=N8) %>%
  melt(id=c("Specie","Serie")) %>% dplyr::select(Specie, Serie, 
                                                 "Posicion"=variable, Base=value)

fly.bases <- fly.bases[fly.bases$Base!=" " & !is.na(fly.bases$Base),] 


# Plot generales
# PLOT 1
ggplot(fly.bases, aes(x=Posicion, fill= Base)) +
  geom_bar(position="fill", na.rm = TRUE, color="black") +
  facet_grid(Serie  ~ Specie , scales="free", space="free") + 
  xlab("Posición en el bucle") +
  ylab("Frecuencia") +
  ggtitle("Frecuencia de aparición de bases en cada sitio del bucle") +
  scale_fill_brewer(palette="Spectral")

# Plot 2
ggplot(fly, aes(x=TerminalPair, y=RnaFoldMFE)) +
  geom_boxplot(outlier.colour="red",
               outlier.size=1, notch=FALSE, stat="boxplot") +
  facet_grid(Serie  ~ Specie) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Mínima energía libre para cada secuencia consenso") + 
  xlab("Par terminal") + 
  ylab("MFE")

# Plot 3
ggplot(fly, aes(x=RnaFoldMFE, stat(density), color=Specie)) +
  geom_freqpoly() +
  facet_grid(TerminalPair ~ Serie)

ggplot(fly, aes(x=Pairments, fill=TerminalPair)) +
  geom_histogram(aes(y=..density..), binwidth = 1, color="black") +
  facet_grid(Serie ~ Specie )

fly <- fly[fly$SerieID ==1,]

ggplot(fly, aes(x=RelativePosition, y=RnaFoldMFE, color=Specie))+ 
  geom_point(alpha=0.5) + geom_smooth() +
  ggtitle("Mínima energía libre vs. posición relativa en la secuencia - Ajuste GAM") + 
  xlab("Posición relativa del stem-loop en la secuencia") + 
  ylab("Mínima energía libre") +
  scale_fill_discrete(name="Listado")

ggplot(fly, aes(x=Serie, fill=TerminalPair)) +
geom_bar(position="fill") +
  scale_y_continuous(labels = scales::percent) + facet_grid(. ~ Specie ) +
  scale_fill_brewer(palette="Spectral")

ggplot(fly, aes(x=Pairments)) +
  geom_histogram(aes(y=..density..), binwidth = 1, color="black") 
+ facet_grid(Serie ~ Specie )
