print("Cargando paquetes...")


################################
# Statistical functions
################################

# if (!require(psych)){
#   install.packages('psych')
#   library(psych)
# }

# if (!require(e1071)){
#   install.packages('e1071')
#   library(e1071)
# }

# if (!require(rsample)){
#   install.packages('rsample')
#   library(rsample)
# }

if (!require(Hmisc)) {
  install.packages('Hmisc')
  library(Hmisc)
}

# if (!require(RVAideMemoire)){
#   install.packages('RVAideMemoire')
#   library(RVAideMemoire)
# }

################################
# String manipulations and data
################################

if (!require(stringr)) {
  install.packages('stringr')
  library(stringr)
}

if (!require(readxl)) {
  install.packages('readxl')
  library(readxl)
}


################################
# Cleaning and tidying data
################################

if (!require(dplyr)) {
  install.packages('dplyr')
  library(dplyr)
}

if (!require(tidyverse)) {
  install.packages('tidyverse')
  library(tidyverse)
}

# Read Excel
if (!require(tidyr)) {
  install.packages('tidyr')
  library(tidyr)
}

# if (!require(reshape2)){
#   install.packages('reshape2')
#   library(reshape2)
# }

################################
# Graphics
################################

if (!require(ggplot2)) {
  install.packages('ggplot2')
  library(ggplot2)
}

if (!require(ggcorrplot)) {
  install.packages('ggcorrplot')
  library(ggcorrplot)
}

if (!require(ggalt)) {
  install.packages('ggalt')
  library(ggalt)
}

if (!require(gridExtra)) {
  install.packages("gridExtra")
  library(gridEXtra)
}

if (!require(scales)) {
  install.packages("scales")
  library(scales)
}

if (!require(wesanderson)) {
  install.packages("wesanderson")
  library(wesanderson)
}

# if (!require(beeswarm)){
#   install.packages('beeswarm')
#   library(beeswarm)
# }

if (!require(ggpubr)) {
  install.packages('ggpubr')
  library(ggpubr)
}


################################
# Miscelaneous functions
################################

if (!require(devtools)) {
  install.packages('devtools')
  library(devtools)
}

if (!require(JLutils)) {
  #source("https://install-github.me/larmarange/JLutils")
  install_github("larmarange/JLutils",force=TRUE)
  library(JLutils)
}

if (!require(DescTools)) {
  install.packages('DescTools')
  library(DescTools)
}


#################################
# Funciones
#################################

print("Cargando funciones...")


# Carga de archivos en dos tandas (Original vs Random)
#######################################

loadFiles <- function(path,
                      pattern = "*.csv",
                      takeRandoms = TRUE) {
  temp <- list.files(path = path, pattern = pattern)
  
  for (i in 1:length(temp)) {
    fileName <- stringr::str_replace(temp[i], ".csv", "")
    print(fileName)
    
    if (grepl("_rnd", fileName)) {
      # Is Random
      if (takeRandoms) {
        if (i == 1) {
          df <- read.csv(
            paste(path, temp[i], sep = ""),
            sep = ";",
            dec = ",",
            stringsAsFactors = FALSE
          )
          df$Serie = 0
        } else {
          aux <- read.csv(
            paste(path, temp[i], sep = ""),
            sep = ";",
            dec = ",",
            stringsAsFactors = FALSE
          )
          aux$Serie = 0
          df <- rbind(df, aux)
        }
      }
      
    } else {
      if (i == 1) {
        df <- read.csv(
          paste(path, temp[i], sep = ""),
          sep = ";",
          dec = ",",
          stringsAsFactors = FALSE
        )
        df$Serie = 1
      } else {
        aux <- read.csv(
          paste(path, temp[i], sep = ""),
          sep = ";",
          dec = ",",
          stringsAsFactors = FALSE
        )
        aux$Serie = 1
        df <- rbind(df, aux)
      }
    }
  }
  
  remove(aux)
  return(df)
}


# Formato de los datos para Ensembl
#######################################
formatEnsembl <- function(df, factorizeAll = TRUE) {
  df$Column4 <- NULL
  df$Column6 <- NULL
  df$LoopID <- NULL
  df$Sense <- NULL
  df$StemLoopSequence <- NULL
  df$PredictedStructure <- NULL
  df$Additional3Seq <- NULL
  df$Additional5Seq <- NULL
  df$ViennaBracketStr <- NULL
  df$EndsAt <- NULL
  df$StartsAt <- NULL
  df$AdditionalSeqMatches <- NULL
  df$AdditionalSeqPositions <- NULL
  df$StartsAt <- NULL
  names(df)[names(df) == 'Column1'] <- 'GenID'
  names(df)[names(df) == 'Column2'] <- 'TranscriptoID'
  names(df)[names(df) == 'Column3'] <- 'GenSymbol'
  names(df)[names(df) == 'Column5'] <- 'Chromosome'
  df$Chromosome <- factor(df$Chromosome)
  
  if (factorizeAll) {
    df$U_PercentSequence <-
      cut(df$U_PercentSequence, seq(0, 1, .25), include.lowest = TRUE)
    df$RelativePosition <-
      cut(df$RelativePosition, seq(0, 1, .20), include.lowest = TRUE)
    df$A_PercentSequence <-
      cut(df$A_PercentSequence, seq(0, 1, .25), include.lowest = TRUE)
    df$C_PercentSequence <-
      cut(df$C_PercentSequence, seq(0, 1, .25), include.lowest = TRUE)
    df$G_PercentSequence <-
      cut(df$G_PercentSequence, seq(0, 1, .25), include.lowest = TRUE)
    df$AU_PercentPairs <-
      cut(df$AU_PercentPairs, seq(0, 1, .25), include.lowest = TRUE)
    df$GU_PercentPairs <-
      cut(df$GU_PercentPairs, seq(0, 1, .25), include.lowest = TRUE)
    df$CG_PercentPairs <-
      cut(df$CG_PercentPairs, seq(0, 1, .25), include.lowest = TRUE)
    df$PurinePercentStem <-
      cut(df$PurinePercentPairs, seq(0, 1, .25), include.lowest = TRUE)
    df$RnaFoldMFE <-
      cut(df$RnaFoldMFE, seq(-30, 0, 5), include.lowest = TRUE)
    df$InternalLoops <-
      cut(df$InternalLoops, seq(0, 5, 1), include.lowest = TRUE)
    df$Bulges <- cut(df$Bulges, seq(0, 5, 1), include.lowest = TRUE)
    df$Pairments <-
      cut(df$Pairments, seq(4, 15, 2), include.lowest = TRUE)
    df$WooblePairs <-
      cut(df$WooblePairs, seq(0, 6, 1), include.lowest = TRUE)
  }
  return(df)
}


# Goodman and Kruskal Tau Measure
#######################################

# (https://www.r-bloggers.com/measuring-associations-between-non-numeric-variables/)
# La tau de Goodman-Kruskal mide la asociación para las tabulaciones cruzadas de
# las variables de niveles nominales.
#
# La tau de Goodman-Kruskal se basa en la asignación aleatoria de las categorías.
# Mide la mejora porcentual en la capacidad de predicción de la variable
# dependiente (variable de columna o fila) dado el valor de otras variables
# (variables de fila o columna). La tau de Goodman-Kruskal es igual a la lambda de
# Goodman-Kruskal salvo que los cálculos del estadístico tau se basan en las
# probabilidades de asignación especificadas por las proporciones marginales o
# condicionales.
#
# Las probabilidades de clasificación incorrecta se basan en la asignación
# aleatoria de las categorías con probabilidades especificadas por la proporción
# marginal o condicional.
# La característica mas importante, y que otorga una asimetría al comparar X-Y con
# Y-X, es que cuantifica la medida en que la variable X es útil para predecir Y.

GKtau <- function(x, y, debug = FALSE) {
  #  First, compute the IxJ contingency table between x and y
  Nij = table(x, y, useNA = "ifany")
  
  if (debug) {
    print("IxJ contingency table:")
    print(Nij)
  }
  
  #  Next, convert this table into a joint probability estimate
  PIij = Nij / sum(Nij)
  
  if (debug) {
    print("IxJ joint probability estimate:")
    print(PIij)
  }
  
  #  Compute the marginal probability estimates
  PIiPlus = apply(PIij, MARGIN = 1, sum)
  PIPlusj = apply(PIij, MARGIN = 2, sum)
  
  if (debug) {
    print("Marginal probability estimates:")
    print(paste("PIiPlus: ", PIiPlus))
    print(paste("PIPlusj: ", PIPlusj))
  }
  
  #  Compute the marginal variation of y
  Vy = 1 - sum(PIPlusj ^ 2)
  
  #  Compute the expected conditional variation of y given x
  InnerSum <- apply(PIij ^ 2, MARGIN = 1, sum)
  VyBarx <- 1 - sum(InnerSum / PIiPlus)
  
  if (debug) {
    print("Compute the expected conditional variation of y given x")
    print(paste("InnerSum:", InnerSum))
    print(paste("VyBarx:", VyBarx))
  }
  
  #  Compute and return Goodman and Kruskal's tau measure
  tau <- (Vy - VyBarx) / Vy
  tau
}


# GK Tau Correlation table
#######################################
cor.GK.tau <- function(df) {
  x <- matrix(nrow = ncol(df), ncol = ncol(df))
  
  rownames(x) <- colnames(df)
  colnames(x) <- colnames(df)
  
  for (j in 1:ncol(df)) {
    for (i in 1:ncol(df)) {
      x[i, j] <- GKtau(df[, i], df[, j])
    }
  }
  
  x <- round(x, 2)
  return(x)
}


# Analisis de la distribucion de todas las bases
#######################################
plotBaseDistribution <- function(df) {
  all.bases <- df %>% dplyr::select(
    Serie,
    "-2" = N.2,
    "-1" = N.1,
    "1" = N2,
    "4" = N5,
    "5" = N6,
    "6" = N7,
    "7" = N8
  ) %>%
    melt(id = c("Serie")) %>%
    dplyr::select(Serie, "Posicion" = variable, Base = value)
  
  all.bases$Serie <- revalue(all.bases$Serie,
                             c("0" = "Secuencias cDNA random",
                               "1" = "Secuencias cDNA originales"))
  
  all.bases <-
    all.bases[all.bases$Base != " " & !is.na(all.bases$Base),]
  all.bases$Base <- factor(all.bases$Base)
  
  # PLOT 1
  pl1 <- ggplot(all.bases, aes(x = Posicion, fill = Base)) +
    geom_bar(
      position = "fill",
      na.rm = TRUE,
      alpha = 0.9,
      color = "black"
    ) +
    facet_grid(. ~ Serie , scales = "free", space = "free") +
    xlab("Posición N") + geom_text(stat = "fill_labels",
                                   size = 3,
                                   fontface = 3) +
    ylab("Proporción") + coord_flip() +
    ggtitle("Bases variables N en el motivo consenso SRE",
            subtitle = "Posiciones NNCNGGN[0..4]. del stem loop. La base C es la posición 0.") +
    scale_fill_brewer(palette = "Pastel2") +
    theme_dark()
  
  return(pl1)
}



# Obtener archivo de la carpeta listas
#######################################
obtenerLista <- function(name) {
  aux.df <-
    read.table(paste("./data/listas/", name, ".txt", sep = ""), header = TRUE)
  return(aux.df)
}
