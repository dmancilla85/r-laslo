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

ggplot(df[df$Serie==0,], aes(x=Pairments)) +
  geom_histogram()


#################
summary(df)
df$RnaFoldMFE <- cut(df$RnaFoldMFE, seq(-20, 0, 5))

library(rpart)
df.fit <- rpart(Serie ~ RelativePosition + RnaFoldMFE + PurinePercentStem +
                  + CG_PercentPairs + AU_PercentPairs + GU_PercentPairs +
                  WooblePairs + Pairments + Bulges + InternalLoops +
                  N.2 + N.1 + N2 + N5 + N6 + N7 + N8 + TerminalPair +
                  LoopPattern + Chromosome, data=df, method="class")

plotcp(df.fit)

plot(df.fit, uniform=TRUE, 
     main="Classification Tree for Kyphosis")
text(df.fit, use.n=FALSE, all=TRUE, cex=.8)

# 
# newdata <- data.frame(RelativePosition= , 
#                       RnaFoldMFE= , 
#                       PurinePercentStem= ,
#                       CG_PercentPairs= ,
#                       AU_PercentPairs= , 
#                       GU_PercentPairs= ,
#                       WooblePairs= , 
#                       Pairments= ,
#                       Bulges= ,
#                       InternalLoops= ,
#                       N.2= ,
#                       N.1= ,
#                       N2= ,
#                       N5= ,
#                       N6= ,
#                       N7= ,
#                       N8= ,
#                       TerminalPair= ,
#                       LoopPattern= ,
#                       Chromosome = )