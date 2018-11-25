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
install_github("larmarange/JLutils")
library(JLutils)
#cluster
#vcd

###############################################################################
Funciones
###############################################################################

loadAllData <- function(path, name, ext=".csv", nRand, serie1, serie2){
  
  print("Cargando archivo... ")
  myData <- read.csv(paste(path, file, ext, sep=""), header=TRUE, sep=";", dec=",", 
                        stringsAsFactors = FALSE)
  myData$SerieID <- 1
  myData <- myData %>% distinct(Column3, Loop, TerminalPair, .keep_all = TRUE)
  myData$Serie = serie1
  
  print("Cargando randoms...")
  
  for(i in 1:nRand){
    df <- read.csv(paste(path, file,"_rnd", i, ext, sep=""), header=TRUE, sep=";", 
                   dec=",", stringsAsFactors = FALSE)
    df$Serie <- serie2 
    df$SerieID <- 2
    myData <- rbind(myData, df)
  }
  
  return(myData)
}
