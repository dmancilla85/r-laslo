print("Instalador de bibliotecas")

if (!require(psych)){
  install.packages('psych')
  library(psych)
}

if (!require(ggplot2)){
  install.packages('ggplot2')
  library(ggplot2)
}

if (!require(RVAideMemoire)){
  install.packages('RVAideMemoire')
  library(RVAideMemoire)
}

if (!require(ggalt)){
  install.packages('ggalt')
  library(ggalt)
}

if (!require(dplyr)){
  install.packages('dplyr')
  library(dplyr)
}

if (!require(tidyr)){
  install.packages('tidyr')
  library(tidyr)
}

if (!require(ggcorrplot)){
  install.packages('ggcorrplot')
  library(ggcorrplot)
}

if (!require(stringr)){
  install.packages('stringr')
  library(stringr)
}
if (!require(DescTools)){
  install.packages('DescTools')
  library(DescTools)
}
if (!require(readxl)){
  install.packages('readxl')
  library(readxl)
}

if (!require(beeswarm)){
  install.packages('beeswarm')
  library(beeswarm)
}

if (!require(devtools)){
  install.packages('devtools')
  library(devtools)
}

if (!require(JLutils)){
  install_github("larmarange/JLutils")
  library(JLutils)
}

if (!require(reshape2)){
  install.packages('reshape2')
  library(reshape2)
}

#cluster
#vcd

###############################################################################
# Funciones
###############################################################################

# 1. Carga de archivos en dos tandas (Original vs Random)
loadFiles <- function(path, pattern = "*.csv"){
  temp <- list.files(path = path, pattern = pattern)
  
  for (i in 1:length(temp)){
    fileName <- stringr::str_replace(temp[i], ".csv", "")
    print(fileName)
    
    if(grepl("_rnd", fileName)){
      # Is Random
      if(i == 1){
        df <- read.csv(paste(path,temp[i],sep=""), sep =";", 
                       dec =",", stringsAsFactors = FALSE)
        df$Serie = 0
      } else {
        aux <- read.csv(paste(path,temp[i],sep=""), sep =";", 
                        dec =",", stringsAsFactors = FALSE)
        aux$Serie = 0
        df <- rbind(df, aux)
      }
      
    } else {
      if(i == 1){
        df <- read.csv(paste(path,temp[i],sep=""), sep =";", 
                       dec =",", stringsAsFactors = FALSE)
        df$Serie = 1
      } else {
        aux <- read.csv(paste(path,temp[i],sep=""), sep =";", 
                        dec =",", stringsAsFactors = FALSE)
        aux$Serie = 1
        df <- rbind(df, aux)
      }
    }
  }
  
  df <- df[!is.na(df$C_PercentSequence), ]
  df$LoopPattern <- factor(df$LoopPattern)
  df$N.1 <- factor(df$N.1)
  df$N.2 <- factor(df$N.2)
  df$N2 <- factor(df$N2)
  df$N5 <- factor(df$N5)
  df$N6 <- factor(df$N6)
  df$N7 <- factor(df$N7)
  df$N8 <- factor(df$N8)
  df$TerminalPair <- factor(df$TerminalPair)
  #df$Serie <- factor(df$Serie)
  df$Bulges <- factor(df$Bulges)
  
  remove(aux)
  return(df)
}
