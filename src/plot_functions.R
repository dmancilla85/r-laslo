###############################################################################
# Funciones generadoras de gráficas
###############################################################################


# Cargar libraries
source("./src/configuration.R")
print("Cargando funciones de gráficos...")

# Pearson correlation
###########################################################
getCorPearson <- function(df) {
  df <- dplyr::select_if(df, is.numeric)
  cor <- rcorr(as.matrix(df, type = "pearson"))
  
  print(summary(cor))
  print(cor)
  
  pl1 <- ggcorrplot(
    cor(df, method = "pearson"),
    lab = TRUE,
    lab_size = 3,
    method = "square",
    tl.cex = 10,
    show.diag = TRUE,
    show.legend = FALSE,
    outline.col = "white",
    colors = c("tomato2", "white", "springgreen3"),
    ggtheme = theme_bw()
  ) +
    ggtitle("Correlación de Pearson",
            subtitle = "Lista de genes unidos a Smaug de Chen")
  
  
  print(pl1)
  return(cor)
  
}

# Goodman-Kruskal Tau correlation
###########################################################
getCorTau <- function(df) {
  corTauDF <- cor.GK.tau(select_if(df, negate(is.numeric)))
  
  pl2 <- ggcorrplot(
    corTauDF,
    lab = TRUE,
    lab_size = 3,
    method = "circle",
    tl.cex = 8,
    show.diag = FALSE,
    show.legend = FALSE,
    outline.col = "white",
    colors = c("tomato2", "white", "springgreen3"),
    ggtheme = theme_bw()
  ) +
    ggtitle("Estadístico Tau de Goodman - Kruskal",
            subtitle = "Lista de genes unidos a Smaug de Chen")
  
  print(pl2)
}

## fix numbers
fixFreq <- function(filtered, myFactor) {
  df <- fixFreq2(filtered, myFactor, "CNGG")
  df <- fixFreq2(df, myFactor, "CNGGN")
  df <- fixFreq2(df, myFactor, "CNGGNN")
  df <- fixFreq2(df, myFactor, "CNGGNNN")
  df <- fixFreq2(df, myFactor, "CNGGNNNN")
  return(df)
}

fixFreq2 <- function(filtered, myFactor, myPattern) {
  sum1 <-
    filtered[filtered$Tipo == myFactor &
               filtered$LoopPattern == myPattern,]$n %>% sum()
  
  filtered[filtered$Tipo == myFactor &
             filtered$LoopPattern == myPattern, ]$freq <-
    filtered[filtered$Tipo == myFactor &
               filtered$LoopPattern == myPattern, ]$n / sum1
  
  return(filtered)
}

###########################################################
featuresOfThePlot <- function(df) {
  filtered <- df %>%  filter(str_length(LoopPattern) <= 8) %>%
    distinct(Tipo, LoopPattern, Pairments, StemLoopSequence,Gen) %>%
    group_by(Tipo, LoopPattern, Pairments) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  limit <- 85
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  filtered$Tipo <- as.factor(filtered$Tipo)
  #write.table(filtered, "test1.txt",sep="\t")
  #View(filtered)
  
  # Plot by pairs
  pairs <- filtered  %>%
    ggplot(aes(
      x = Pairments,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    geom_point(aes(color = Tipo), show.legend = T, alpha = 0.7) +
    ylab("%") + #labs(color = "Listado") +
    geom_line(linetype = "dashed",
              show.legend = F,
              alpha = 0.3) +
    xlab("Longitud del stem") +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    facet_grid(~ LoopPattern) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  # Plot by wooble
  filtered <- df %>% filter(str_length(LoopPattern) <= 8) %>%
    distinct(Tipo, LoopPattern, WooblePairs, Gen, StartsAt) %>%
    group_by(Tipo, LoopPattern, WooblePairs) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  
  
  wooble <- filtered %>%
    ggplot(aes(
      x = WooblePairs,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    geom_line(linetype = "dashed",
              show.legend = F,
              alpha = 0.3) +
    xlab("Apareamientos GU") +
    geom_point(aes(color = Tipo), show.legend = F, alpha = 0.7) +
    ylab("%") + labs(color = "Listado") +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    facet_grid(~ LoopPattern) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  # Plot by bulges
  filtered <- df %>% filter(str_length(LoopPattern) <= 8) %>%
    distinct(Tipo, LoopPattern, Bulges, Gen, StartsAt) %>%
    group_by(Tipo, LoopPattern, Bulges) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  
  bulges <- filtered %>%
    ggplot(aes(
      x = Bulges,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    geom_line(linetype = "dashed",
              show.legend = F,
              alpha = 0.3) +
    xlab("Bulges") +
    geom_point(aes(color = Tipo), show.legend = F, alpha = 0.7) +
    ylab("%") + labs(color = "Listado") +
    facet_grid(~ LoopPattern) +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  # Plot by internal loops
  filtered <- df %>% filter(str_length(LoopPattern) <= 8) %>%
    distinct(Tipo, LoopPattern, InternalLoops, Gen, StartsAt) %>%
    group_by(Tipo, LoopPattern, InternalLoops) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  
  internals <- filtered %>%
    ggplot(aes(
      x = InternalLoops,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    geom_line(linetype = "dashed",
              show.legend = F,
              alpha = 0.3) +
    xlab("Loops internos") +
    geom_point(aes(color = Tipo), show.legend = F, alpha = 0.7) +
    ylab("%") + labs(color = "Listado") +
    facet_grid(~ LoopPattern) +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  # Arrange grid with all plots
  print(grid.arrange(
    pairs,
    wooble,
    bulges,
    internals,
    ncol = 1,
    top = text_grob(
      "% SRE detectados por cada gen",
      face = "bold",
      size = 13.5
    ),
    bottom = ""
  ))
}


###########################################################
geneHitsByPattern <- function(df) {
  limit <- 70
  
  filtered <- df %>%  filter(str_length(LoopPattern) <= 8) %>%  
    distinct(Tipo, LoopPattern, N.1,StemLoopSequence,Gen) %>%
    group_by(Tipo, LoopPattern, N.1) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  
  write.table(filtered, "test1.txt", sep = "\t")
  
  N2 <- filtered %>%
    ggplot(aes(
      x = N.1,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    #geom_line(linetype = "dashed") +
    xlab("Base predecesora (N-1)") +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    geom_point(aes(size=Tipo),
               alpha = 0.5,
               show.legend = T) +
    ylab("%") + #labs(color = "Listado") +
    geom_segment(aes(
      x = N.1,
      xend = N.1,
      y = 0,
      yend = freq * 100
    ), alpha = 0.5) +
    facet_grid( ~ LoopPattern) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  filtered <- df %>% filter(str_length(LoopPattern) <= 8) %>%
    distinct(Tipo, LoopPattern, N2, Gen) %>%
    group_by(Tipo, LoopPattern, N2) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  filtered <- fixFreq(filtered, "Unidos y silenciados")
  filtered <- fixFreq(filtered, "No unidos y silenciados")
  filtered <- fixFreq(filtered, "Unidos y no silenciados")
  filtered <- fixFreq(filtered, "No unidos y no silenciados")
  filtered <- fixFreq(filtered, "Random")
  
  n2 <- filtered %>%
    ggplot(aes(
      x = N2,
      y = freq * 100,
      group = Tipo,
      color = Tipo
    )) +
    #geom_line(linetype = "dashed", show.legend = F) +
    xlab("Base N2") +
    geom_point(aes(size = Tipo), alpha = 0.5, show.legend = T) +
    geom_segment(aes(
      x = N2,
      xend = N2,
      y = 0,
      yend = freq * 100
    ), alpha = 0.5) +
    ylab("%") + #labs(color = "Listado") +
    facet_grid( ~ LoopPattern) +
    scale_y_continuous(breaks = seq(0, limit, by = 10),
                       limits = c(0, limit)) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "floralwhite")
    )
  
  print(grid.arrange(
    N2,
    n2,
    ncol = 1,
    top = text_grob(
      "% Genes acertados por patrón del loop",
      face = "bold",
      size = 13.5
    ),
    bottom = ""
  ))
}


###########################################################
locationAndPairs <- function(df) {
  # Chromosome
  chromo <- df %>% #filter(Tipo != "Random") %>%
    ggplot(aes(x = Tipo, fill = Location)) +
    geom_bar(color = "darkgray",
             position = "fill") + theme_gray() +
    ggtitle("Localización de los stems encontrados",
            subtitle = "") + labs(fill = "Location") +
    geom_text(
      stat = "fill_labels",
      size = 3,
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_manual(values = c('blue2', 'brown3', 'orange', "yellow", "green")) +
    xlab("Posición") + ylab("Recuento") #480x480px
  
  
  # TerminalPair
  pairs <- df %>%
    ggplot(aes(x = Tipo, fill = TerminalPair)) +
    geom_bar(position = "fill",
             color = "darkgray") + labs(fill = "Par de cierre") +
    ggtitle("Tipo de apareamiento de cierre 3' a 5'",
            subtitle = "") +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#DCDCDC"),
      panel.border = element_rect(colour = "darkgrey", fill = NA)
    ) + xlab("Grupo") +
    geom_text(
      stat = "fill_labels",
      size = 3,
      color = "white",
      fontface = "bold"
    ) + #stat_stack_labels() +
    scale_fill_manual(values = wes_palette(6, name = "Rushmore1", type = "continuous"))
  
  print(grid.arrange(chromo, pairs))
}

plotLoopPatterns <- function(df) {
  # LoopPattern
  
  patt <- df %>%
    ggplot(aes(fill = LoopPattern, x = 1)) +
    geom_bar(position = "fill",
             color = "darkgrey",
             alpha = 0.7) + coord_polar(theta = "y") +
    geom_text(
      stat = "fill_labels",
      fontface = "bold",
      position = position_dodge(width = .8),
      size = 3.5,
      check_overlap = TRUE,
      
    ) + labs(fill = "Secuencia") +
    facet_wrap( ~ Tipo, nrow=1) + ylab("") +
    ggtitle("Motivos buscados en el Loop",
            subtitle = "Sobre todos los stem-loops hallados en los transcriptos.") +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#DCDCDC"),
      panel.border = element_rect(colour = "darkgrey", fill = NA),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    ) +
    ylab("") + xlab("") #+
    #scale_fill_manual(values = wes_palette(7, name = "BottleRocket1", type = "continuous"))
  
  # "#3B9AB2" "#78B7C5" "#EBCC2A" "#E1AF00" "#F21A00"
  
  print(patt)
}


###########################################################
plotNBases <- function(df) {
  windowsFonts(Calibri = windowsFont("TT Calibri Light"))
  
  myPlot <- df %>%
    select (Tipo, N.2, N.1, N2, N5, N6, N7, N8) %>%
    gather(metric, value, -Tipo) %>%
    filter(value %in% c("A", "G", "C", "U")) %>%
    mutate(metric = str_replace_all(metric, "N.1", "N(-1)")) %>%
    mutate(metric = str_replace_all(metric, "N.2", "N(-2)")) %>%
    ggplot(aes(fill = value, x = 1)) +
    geom_bar(position = "fill",
             color = "darkgrey",
             alpha = 0.8) + coord_polar(theta = "y") +
    labs(fill = "Base") +
    ggtitle("Frecuencia de bases en cada posición variable (N) del loop",
            subtitle = "Sobre todos los stem-loops hallados en los transcriptos.") +
    geom_text(
      stat = "fill_labels",
      fontface = "bold",
      position = position_dodge(width = .8),
      size = 3.1,
      check_overlap = TRUE,
      color = "black",
      family = "Calibri"
    ) +
    facet_grid(Tipo ~ metric) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.background = element_rect(fill = "beige"),
      panel.border = element_rect(colour = "darkgrey", fill = NA),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "seashell")
    ) +
    scale_fill_manual(values = wes_palette("IsleofDogs2"))
  print(myPlot)
}


###########################################################
wpbPlots <- function(df) {
  
  filtered <- df %>%
    select(Tipo, WooblePairs, Bulges, InternalLoops) %>%
    gather(metric, value, -Tipo) %>%
    mutate(metric = str_replace_all(metric, "WooblePairs", "Pares GU")) %>%
    mutate(metric = str_replace_all(metric, "InternalLoops", "Loops internos"))
  
  myPlot <- filtered %>%
    ggplot(aes(x = value, fill = metric)) +
    geom_histogram(
      show.legend = F,
      bins = 9,
      color = "black",
      lwd = 0.8,
      aes(y = stat(width*density),group=Tipo)) +
    labs(fill = "") +
    theme(
      panel.background = element_rect(fill = "azure"),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.background = element_rect(fill = "whitesmoke")
    ) + 
    scale_y_continuous(labels=percent_format()) +
    scale_x_continuous(breaks = 0:9) +
    ggtitle("Propiedades de los stem-loops", subtitle = "") +
    xlab("Cantidad") + ylab("Stem-loops") +
    facet_grid(Tipo ~ metric) +
    scale_fill_manual(values = wes_palette("Darjeeling2"))
  
  print(myPlot)
}


###########################################################
basePercentsPlot <- function(df) {
  
  filtered <- df %>% filter(Tipo != "Random") %>%
    select(Tipo,
           A_PercentSequence,
           C_PercentSequence,
           G_PercentSequence,
           U_PercentSequence) %>%
    gather(metric, value, -Tipo) %>%
    mutate(metric = paste ("", str_replace(metric, "_PercentSequence", "")), value = value * 100)
  
  #View(filtered)
  
  myPlot <- filtered %>%
    ggplot(aes(x = value, fill = metric)) + ylab("Densidad") +
    geom_density(
      aes(y = ..density..),
      position = "identity",
      #bins = 40,
      alpha = 0.3,
      color = "black"
    ) +
    xlab("Ratio de cada nucleótido en la secuencia entera") +
    facet_wrap( ~ Tipo,nrow=1) + theme_pubr() +
    ggtitle("Composición de las secuencias cDNA",
            subtitle = "Sobre el total de bases en transcriptos") +
    scale_fill_manual(values = wes_palette("Cavalcanti1")) +
    scale_x_continuous(breaks=seq(1,50,by=10))
  
  print(myPlot)
}


###########################################################
pairsMFEPlot <- function(df) {
  # Percent Pairs in Stem
  pairs2 <- df %>% #filter(Tipo != "Random") %>%
    select(Tipo, AU_PercentPairs, CG_PercentPairs, GU_PercentPairs) %>%
    gather(metric, value, -Tipo) %>%
    mutate(metric = paste ("", str_replace(metric, "_PercentPairs", "")), value = value * 100) %>%
    ggplot(aes(y = value, x = Tipo)) + ylab("Densidad") +
    geom_boxplot(
      aes(colour = metric),
      outlier.colour = "red",
      outlier.shape = 1,
      lwd = 1.3
    ) +
    labs(colour = "Apareamiento") +
    xlab("Ratio de cada apareamiento en los stem-loop") +
    ggtitle("Composición de apareamientos en los stem-loops",
            subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
    theme(panel.border = element_rect(colour = "darkgrey", fill = NA)) +
    scale_colour_manual(values = wes_palette("FantasticFox1"))
  
  # Minimum free energy
  mfe <- df %>%
    select(Tipo, RnaFoldMFE) %>%
    ggplot(aes(x = RnaFoldMFE, group = Tipo)) + ylab("Densidad") +
    geom_density(alpha = 1, size = 0.7) +
    facet_wrap(~ Tipo, nrow=1) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    scale_x_continuous(breaks = round(seq(
      min(fly$RnaFoldMFE),
      max(fly$RnaFoldMFE), by = 4
    ), 1))  +
    scale_y_continuous(breaks = seq(0, 0.2, by = 0.02)) +
    xlab("Mínima energía libre (kcal/mol)") +
    ggtitle("Energía libre en los stem-loops",
            subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
    geom_vline(
      data = filter(fly, Tipo == "Unidos y no silenciados"),
      aes(xintercept = median(RnaFoldMFE), group = Tipo),
      color = "red",
      linetype = "dashed",
      show.legend = F,
      size = 1
    ) +
    geom_vline(
      data = filter(fly, Tipo == "Unidos y silenciados"),
      aes(xintercept = median(RnaFoldMFE), group = Tipo),
      color = "blue",
      linetype = "dashed",
      show.legend = F,
      size = 1
    ) +
    geom_vline(
      data = filter(fly, Tipo == "No unidos y silenciados"),
      aes(xintercept = median(RnaFoldMFE), group = Tipo),
      color = "darkgreen",
      linetype = "dashed",
      show.legend = F,
      size = 1
    )+
    geom_vline(
      data = filter(fly, Tipo == "No unidos y no silenciados"),
      aes(xintercept = median(RnaFoldMFE), group = Tipo),
      color = "darkgreen",
      linetype = "dashed",
      show.legend = F,
      size = 1
    )+
    geom_vline(
      data = filter(fly, Tipo == "Random"),
      aes(xintercept = median(RnaFoldMFE), group = Tipo),
      color = "darkgreen",
      linetype = "dashed",
      show.legend = F,
      size = 1
    )
  
  print(grid.arrange(pairs2, mfe))
}


###########################################################
purinePlot <- function(df) {
  # PurinePercentPairs
  myPlot <- df %>%
    select(Tipo, PurinePercentPairs) %>%
    ggplot(aes(x = PurinePercentPairs, group = Tipo)) + ylab("Densidad") +
    geom_density(alpha = 1, size = 1) +
    facet_wrap(~ Tipo) +
    xlab("Ratio de A-G en apareamientos de los stems") +
    theme_pubr() +
    theme(
      strip.background = element_rect(colour = "white", fill = "lightgrey"),
      strip.text = element_text(face = "bold")
    ) +
    ggtitle("Contenido de bases A-G (purinas) en los stem-loops",
            subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
    geom_vline(
      data = fly,
      aes(
        xintercept = mean(PurinePercentPairs),
        color = "blue"
      ),
      linetype = "dashed",
      show.legend = F,
      size = 1
    )
  print(myPlot)
}


###########################################################
relativePositionPlot <- function(df) {
  # Cambiar por posicion UTR/CDS
  myPlot <-
    df %>% ggplot(aes(
      color = Tipo,
      x = RelativePosition,
      y = (-1) * RnaFoldMFE
    )) + 
    geom_smooth(method = "auto", se=F) +
    ggtitle("Variación de la estabilidad con respecto a la posición relativa en la secuencia") +
    theme_pubr()
  
  print(myPlot)
}

###########################################################
locationPlot <- function(df) {
  # Cambiar por posicion UTR/CDS
  myPlot <-
    df %>% ggplot(aes(
      fill = Location,
      x = Tipo,
      y = (-1) * RnaFoldMFE
    )) +
    geom_boxplot()
  ggtitle("Variación de la estabilidad con respecto a la posición relativa en la secuencia") +
    theme_pubr()
  
  print(myPlot)
}

print("Listo.")