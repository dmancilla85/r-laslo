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
Funciones
###############################################################################

# 1. Carga de archivos en dos tandas (Original vs Random)
loadFiles <- function(path, pattern = "*.csv"){
  temp <- list.files(path = path, pattern = pattern)
  
  for (i in 1:length(temp)){
    fileName <- str_replace(temp[i], ".csv", "")
    print(fileName)
    
    if(grepl("_rnd", fileName)){
      # Is Random
      if(i == 1){
        df <- read.csv(temp[i], sep =";", dec =",", stringsAsFactors = FALSE)
        df$Serie = 2
      } else {
        aux <- read.csv(temp[i], sep =";", dec =",", stringsAsFactors = FALSE)
        aux$Serie = 2
        df <- rbind(df, aux)
      }
      
    } else {
      if(i == 1){
        df <- read.csv(temp[i], sep =";", dec =",", stringsAsFactors = FALSE)
        df$Serie = 1
      } else {
        aux <- read.csv(temp[i], sep =";", dec =",", stringsAsFactors = FALSE)
        aux$Serie = 1
        df <- rbind(df, aux)
      }
    }
  }
  
  remove(aux)
  return(df)
}
