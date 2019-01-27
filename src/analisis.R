
source("./src/configuration.R")

df <- loadFiles(path="./data/fruitfly_bound/", pattern = "*.csv")

#1 Eliminar variables innecesarias
# a. correlaciones
getCors(df)

###############################################################################
# 2. MFE Energy

laslo <- df
shapiro.test(laslo[laslo$Serie=="Original", ]$RnaFoldMFE)
shapiro.test(laslo[laslo$Serie=="Random", ]$RnaFoldMFE)

ksMFE <- wilcox.test(laslo[laslo$Serie=="Original", ]$RnaFoldMFE, 
                     laslo[laslo$Serie=="Random", ]$RnaFoldMFE)

laslo %>%
  ggplot(aes(x=RnaFoldMFE, fill=Serie)) + 
  geom_density(alpha=0.3) + 
  xlab("Mínima energía libre (kcal/mol)") +
  ylab("Densidad") +
  scale_x_continuous(breaks=pretty(laslo$RnaFoldMFE, n = 10)) +
  #labs( caption = caption_fly) +
  ggtitle("Mínima energía líbre de los motivos detectados (calculado por RNAFold)", 
          subtitle = paste(ksMFE$method, "p-value:", formatC(ksMFE$p.value, format="e", digits=2))) +
  scale_fill_discrete(name="Serie") + theme_bw() + theme(legend.position="bottom")

###############################################################################

# 3.a Pares
IQR(laslo[laslo$Serie=="Original", ]$Bulges)
IQR(laslo[laslo$Serie=="Random", ]$Bulges)

wt <- wilcox.test(laslo[laslo$Serie=="Original", ]$Bulges, 
                  laslo[laslo$Serie=="Random", ]$Bulges)

laslo %>%
  ggplot(aes(x=Serie, y=Bulges, fill=Serie)) +
  geom_boxplot(varwidth = T, show.legend = FALSE) +
  ylab("Bases apareadas (bp)") + xlab(NULL) +
  ggtitle("Cantidad de bulges en el tallo (stem) de las horquillas",
          subtitle= paste(wt$method, "p-value:", formatC(wt$p.value, format="e", digits=2))) 

# 3.a Pares
IQR(laslo[laslo$Serie=="Original", ]$InternalLoops)
IQR(laslo[laslo$Serie=="Random", ]$InternalLoops)

wt <- wilcox.test(laslo[laslo$Serie=="Original", ]$InternalLoops, 
                  laslo[laslo$Serie=="Random", ]$InternalLoops)

laslo %>%
  ggplot(aes(x=Serie, y=InternalLoops, fill=Serie)) +
  geom_boxplot(varwidth = T, show.legend = FALSE) +
  ylab("Bases apareadas (bp)") + xlab(NULL) +
  ggtitle("Cantidad de loops internos en el tallo (stem) de las horquillas",
          subtitle= paste(wt$method, "p-value:", formatC(wt$p.value, format="e", digits=2))) 

###############################################################################

# 7. LoopPattern 
# fix
# laslo[laslo$Sense=="-" & str_sub(laslo$LoopPattern,1,1) != "C", ]$LoopPattern <-
#   StrRev(laslo[laslo$Sense=="-" & str_sub(laslo$LoopPattern,1,1) != "C", ]$LoopPattern)

laslo$LoopPattern <- factor(laslo$LoopPattern)  
summary(laslo[laslo$Serie=="Original", ]$LoopPattern)
summary(laslo[laslo$Serie=="Random", ]$LoopPattern)


prop.table(table(laslo$LoopPattern, laslo$Serie), 2)

tb <- as.data.frame(table(laslo$LoopPattern, laslo$Serie))
tb <- tb %>%
  spread(Var2,Freq) %>%
  mutate(Observed = `Original`/ sum(`Original`)) %>%
  mutate(Expected = `Random`/ sum(`Random`)) 
tb

# Hipotesis nula: Las proporciones son iguales para los dos sets
a <- prop.test(tb$`Original`, tb$`Random`)
# X-squared = 352.28, df = 4, p-value < 2.2e-16
# alternative hypothesis: two.sided
# sample estimates:
#   prop 1     prop 2     prop 3     prop 4     prop 5 
# 0.07622704 0.19134304 0.06760349 0.04407713 0.08033419 

# Hipotesis nula: Las proporciones no son sign. diferentes
g <- GTest(tb$`Original`, p=tb$Expected) #p.val < 2.2e-16

laslo %>% ggplot(aes(x=Serie, fill=LoopPattern)) + 
  geom_bar(position="fill", alpha= 0.9, 
           width = 0.5, color="brown") + 
  ylab("Frecuencia relativa")+ geom_text(stat="fill_labels",
                                         size=3) + xlab(NULL) +
  scale_fill_discrete(name="Secuencia consenso") + coord_flip() +
  ggtitle("Porcentaje de las variantes del motivo consenso", 
          subtitle = paste(g$method, "p-value:", formatC(g$p.value, format="e", digits=2))) + 
  theme_bw() # + theme(aspect.ratio = 1/2.5)

###############################################################################

# 9- Analisis de la distribucion de todas las bases
all.bases <- laslo %>% dplyr::select(Serie, "N(-1)"=N.1, "N(-2)"=N.2, 
                                          "N2"=N2, "N5"=N5, "N6"=N6,"N7"=N7,"N8"=N8) %>%
  melt(id=c("Serie")) %>% 
  dplyr::select(Serie, "Posicion"=variable, Base=value)

all.bases <- all.bases[all.bases$Base!=" " & !is.na(all.bases$Base),] 

# PLOT 1
ggplot(all.bases, aes(x=Posicion, fill= Base)) +
  geom_bar(position="fill", na.rm = TRUE, alpha=0.9, color="black") +
  facet_grid( . ~ Serie , scales="free", space="free") + 
  xlab("Base variable N") + geom_text(stat = "fill_labels", 
                                      size = 3, fontface=2) +
  ylab("Proporción") + coord_flip() +
  ggtitle("Bases variables de la secuencia consenso SRE",
          subtitle="Loops que aparecen en más del 5% de los transcriptos") + 
  scale_fill_brewer()


##################################
# Loops - MUY INTERESANTE
lfilt <- laslo %>%
distinct(Serie,Gen, Note, Loop, Location)

tb <- as.data.frame(table(lfilt$Loop, lfilt$Serie))

tb <- tb %>%
  spread(Var2,Freq) %>%
  mutate(Observed = `Original`/ sum(`Original`)) %>%
  mutate(Expected = `Random`/ sum(`Random`)) 
tb <- tb %>% filter(`Original` > 0) %>%
  arrange(desc(`Original`))
tb

tb$proptest <- 0

for (i in 1:nrow(tb)){
  tb[i, ]$proptest <- prop.test(x = c(tb[i,]$Original, tb[i,]$Random), 
                           n = c(sum(tb$Original),sum(tb$Random)), correct=TRUE)$p.value
}

tb$padj <- p.adjust(tb$proptest, method = "bonferroni", n = length(tb$proptest))
tb

write.table(tb, "table1", sep="\t")

##################################
# Loops - MUY INTERESANTE
lfilt <- laslo %>%
  distinct(Serie,Gen, Note, Loop, Location, TerminalPair)

tb <- as.data.frame(table(lfilt$Gen, lfilt$Serie))

tb <- tb %>%
  spread(Var2,Freq) %>%
  mutate(Observed = `Original`/ sum(`Original`)) %>%
  mutate(Expected = `Random`/ sum(`Random`)) 
tb <- tb %>% filter(`Original` > 0) %>%
  arrange(desc(`Original`))
tb

tb$proptest <- 0

for (i in 1:nrow(tb)){
  tb[i, ]$proptest <- prop.test(x = c(tb[i,]$Original, tb[i,]$Random), 
                                n = c(sum(tb$Original),sum(tb$Random)), correct=TRUE)$p.value
}

tb$padj <- p.adjust(tb$proptest, method = "bonferroni", n = length(tb$proptest))
write.table(tb, "table2", sep="\t")

################################

library(beeswarm)
laslo %>% nrow()

beeswarm(RnaFoldMFE ~ Location, data = laslo[laslo$Serie=="Random",],
         method = 'swarm',
         main = paste("Energía libre por posición en transcriptos con UTR informado (", 
                      nrow(laslo[laslo$Serie=="Random",])," horquillas)", sep=""),
         pch = 16, pwcol = as.numeric(TerminalPair),
         xlab = '', ylab = 'Mínima energía libre(kcal/mol)',
         labels = levels(laslo[laslo$Serie=="Random",]$Location))

ggplot(laslo, aes(x=Serie, fill=Location)) + ggtitle("Posición en el transcripto")+
  geom_bar(position = "fill",na.rm = TRUE, color="black", width = 0.5) + theme_minimal()
