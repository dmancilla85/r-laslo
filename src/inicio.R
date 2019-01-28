# Cargar framework
source("./src/configuration.R")

# Cargar datos
df <- loadFiles(path="./data/fruitfly_bound/", pattern = "*.csv")

#Eliminar variables innecesarias
df <- formatEnsembl(df)

# a. correlaciones
getCors(df)
plotBaseDistribution(df)
df$Serie <- factor(df$Serie)

ggdensity(df, x = "RnaFoldMFE",
          add = "mean", rug = TRUE,
          color = "Serie", fill = "Serie",
          palette = c("#00AFBB", "#E7B800"))
