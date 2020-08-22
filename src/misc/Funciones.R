# Seteamos directorio de trabajo
ls()
rm(list=ls())

# Cargar packages ################################################################
if (!require("moments")) {
  install.packages("moments", dependencies = TRUE)
  library(moments)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("DAAG")) {
  install.packages("DAAG", dependencies = TRUE)
  library(DAAG)
}

# if (!require("stringr")) {
#   install.packages("stringr", dependencies = TRUE)
#   library(stringr)
# }

if (!require("corpcor")) {
  install.packages("corpcor", dependencies = TRUE)
  library(corpcor)
}

if (!require("mctest")) {
  install.packages("mctest", dependencies = TRUE)
  library(mctest)
}

if (!require("GGally")) {
  install.packages("GGally", dependencies = TRUE)
  library(GGally)
}

if (!require("pastecs")) {
  install.packages("pastecs", dependencies = TRUE)
  library(pastecs)
}

if (!require("Hmisc")) {
  install.packages("Hmisc", dependencies = TRUE)
  library(Hmisc)
}

if (!require("scatterplot3d")) {
  install.packages("scatterplot3d", dependencies = TRUE)
  library(scatterplot3d)
}

if (!require("MASS")) {
  install.packages("MASS", dependencies = TRUE)
  library(MASS)
}

# if (!require("klaR")) {
#   install.packages("klaR", dependencies = TRUE)
#   library(klaR)
# }

# if (!require("DiscriMiner")) {
#   install.packages("DiscriMiner", dependencies = TRUE)
#   library(DiscriMiner)
# }

if (!require("car")) {
  install.packages("car", dependencies = TRUE)
  library(car)
}

# if (!require("FactoMineR")) {
#   install.packages("FactoMineR", dependencies = TRUE)
#   library(FactoMinerR)
# }

# if (!require("factoextra")) {
#   install.packages("factoextra", dependencies = TRUE)
#   library(factoextra)
# }

if (!require("psych")) {
  install.packages("psych", dependencies = TRUE)
  library(psych)
}

if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}

if (!require("vegan")) {
  install.packages("vegan", dependencies = TRUE)
  library(vegan)
}

if (!require("MuMIn")) {
  install.packages("MuMIn", dependencies = TRUE)
  library(MuMIn)
}

if (!require("rrcovHD")) {
  install.packages("rrcovHD", dependencies = TRUE)
  library(rrcovHD)
}

# #############################################################################
# Test t para muestras independientes / apareadas
# t.test(x, alternative = c("two.sided", "less", "greater"), paired = FALSE, var.equal = FALSE)
# t.test(y ~ x, data = d, var.equal=TRUE) 

# Prueba F para dos varianzas
# var.test(x, alternative = c("two.sided", "less", "greater))
#                             Test de Wilcoxon para muestras independientes/apareadas
#                             wilcox.test(x, alternative = c("two.sided", "less", "greater), paired = FALSE)

# Test z para dos proporciones
# prop.test(c(r1,r2),c(n1,n2), alternative = c("two.sided", "less", "greater") 

# tapply(y, x, summary)  # estadistica descriptiva
# xtabs( ~ y, data = d)  #tabla de frecuencias
# xtabs( ~ x + y, data = d)  #tabla de frecuencias de doble entrada

# Análisis de la varianza no paramétrico (Kruskal-Wallis)
# kruskal.test(y ~ x, data = nombre) 

# Prueba chi-cuadrado
# chisq.test(table(A,B))
# #############################################################################

calcularRegresion <- function(X, Y, grado=1) {
  modelo1<- lm(Y ~ poly(X, grado))
  print("Modelo")
  print(modelo1)
  residuos<-residuals(modelo1)
  predichos<-predict(modelo1)
  #Calculo los residuos standard para el grafico de residuos vs predichos
  res_standard<-rstandard(modelo1)
  tabla <- cbind(Dosis=X, "Cd en tallos"=Y, Residuos=residuos, Predichos=predichos, "Residuos estandarizados"=round(res_standard,2))
  print("Tabla")
  print(tabla)
  
  par(mfrow=c(2,2))
  #plot(X, Y)
  #lines(lowess(X, ))
  hist(residuos)
  boxplot(residuos)
  qqnorm(residuos)
  #qqline(residuos)
  #shapiro.test(residuos)
  plot(predichos, res_standard, main="Residuos vs Predichos", xlab="Predichos", ylab="Residuos Estandards")
  abline(0, 0)
  abline(-2, 0)
  abline(2, 0)
  
  print("Summary")
  print(summary(modelo1))
}

# Nombre      : Template1
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_RLM_Diagnostico = function(model, influence = FALSE, crossVal = FALSE){
	# Limpiar consola
	cat("\014")
	
	# La tabla anova nos confirma que las variables explicativas de nuestro modelo 
	# son significativas
	anova(model)

	# 1. Validación del modelo - Criterios de selección 
	# 1.a Minima varianza residual ########
	# El residuos deben estar distribuidos al azar alrededor de la línea horizontal
	# que representa un error residual de cero
	# Si el gráfico tiene forma de embudo
	# entonces lo más probable es que exista heterocedastididad en los datos
	plot(model, which = 1, pch = 20)

	# 1.b Verificar distribución normal ###
	plot(model, which = 2, pch = 20)
	# Shapiro-Wilks
	print(shapiro.test(residuals(model)))

	# 1.c Esta gráfica se utiliza para detectar si la difusión de los residuos 
	# es constante en el rango de valores ajustados
	plot(model, which = 3, pch = 20)

	# 1.d Valor leverage ################## 
	# Valor de cada punto, la medida de su importancia en la determinación 
	# del modelo de regresión.
	# Distancias pequeñas significan que la eliminación de la observación tiene poco 
	# efecto sobre los resultados de la regresión y distancias mayores a 1 son 
	# sospechosas, sugieren la presencia de un posible valor atípico o de un modelo pobre.
	plot(model, which = 5, pch = 20)

	# 1. Valor F: la variabilidad explicada por el modelo es mayor que la que se queda sin explicar.
	# Hipótesis nula -> F = 1 por lo que F tiene que ser mayor que 1
	print(paste("Valor F: ", round(summary(model)$fstatistic[1], 3)))

	# 2. Máximo R2 fitted 
	# Esto que significa que la recta de regresión explica el (R2.adjx100)% 
	# de la variabilidad del modelo
	print(paste("R-cuadrado: ", round(summary(model)$r.squared, 3)))

	# 3. Mínimo valor de Akaike (Medida de la calidad relativa de un modelo). 
	#   es una medida de ajuste que penaliza el modelo por tener más variables,
	# a diferencia del r2 que crece a mayor cantidad de variables.
	print(paste("Akaike: ", round(AIC(model), 3)))
	# 4. No Colinearidad
	print("VIF: ")
	print(vif(model))
	print(sqrt(vif(model)) > 2)
	
	###############################################################################
	# Extra - Análisis de la influencia de los puntos
	if(influence){
		infl <- influence.measures(model)
		summary(infl)
		influencePlot(model, id.n = 2)
	}

	# Validación cruzada
	if(crossVal) {
		cv.lm(df, model, m = 2)
	}
	###############################################################################

}


# Nombre      : anovaParam
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
anovaParam <- function(x, y){
  model <- aov(y ~ x)
  summary(model)
  qqnorm(residuals(model))
  qqline(residuals(model))
  shapiro.test(residuals(model))
  plot(model)
  leveneTest(x, y)
  TukeyHSD(model)
}

# Nombre      : anovaNonParam
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
anovaNonParam <- function(x, y, df){
  kruskal.test(y ~ x, data = df)
}

# Nombre      : regresionLineal
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_RegresionLineal <- function(x, y, df){
  # Limpiar consola
  cat("\014")
  
  plot(x, y)
  model <- lm(y ~ x, data = df)
  abline(model)
  print(summary(model))
  
  # Test de normalidad
  print(qqnorm(residuals(model)))
  print(qqline(residuals(model)))
  # Shapiro-Wilks
  print(shapiro.test(residuals(model)))
  
  print(plot(model))
  
  print("Modelo predicho: ")
  print(predict(model))
  
  print("Intervalo de confianza: ")
  print(round(confint(model), 2))
  
  print(paste("R-cuadrado: ", round(summary(model)$r.squared), 3))
  print(paste("Akaike: ", round(AIC(model), 3)))
}


# Nombre      : regresionContyCategSI
# Creado por  : DM
# Descripción : Modelo con dos VE (una continua y una explicatoria)
#               Crear dummy. Modelo sin interacción
# Ejemplo     : 
regresionContyCategSI <- function(y, x1, x2, df){
  model <- lm(y ~ x1 + x2, df)
  print(model)
  print(paste("R-cuadrado: ", round(summary(model)$r.squared), 3))
  print(paste("Akaike: ", round(AIC(model), 3)))
  
  # Colinealidad (Valores > 5 indican colinealidad)
  print(paste("VIF: ", vif(model)))
}

# Modelos de regresión multivariados
# Criterios de selección
# 1. Minima varianza residual
# 2. Máximo R2 fitted
# 3. Mínimo valor de Akaike
# 4. Retener variables con coeficientes significativos

# Nombre      : regresionContyCategCI
# Creado por  : DM
# Descripción : Modelo con dos VE (una continua y una explicatoria)
#               Crear dummy. Modelo con interacción
# Ejemplo     : 
regresionContyCategCI <- function(y, x1, x2, df){
  model <- lm(y ~ x1 * x2, df)
  print(model)
  print(paste("R-cuadrado: ", round(summary(model)$r.squared), 3))
  print(paste("Akaike: ", round(AIC(model), 3)))
  
  # Colinealidad (Valores > 5 indican colinealidad)
  print(paste("VIF: ", vif(model)))
}

# Nombre      : NAXMedia
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_RemoveNAs = function(df, useLm = FALSE){
  # Tomar numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  
  # Manejar los NA's
  # Opciones:
  # 1- Tomar promedio (pocos valores)
  # 2- Regresión lineal (muchos valores)
  
  # Tomar promedio
  for(i in 1:nrow(df)) {
    
    for(j in 1:ncol(df)){
      if(is.na(df[i, j]))
        df[i,j] <- mean(df[, j], na.rm = TRUE)
    }
  }
  
  return(df)
}

# Nombre      : DM_Univariada
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_Univariada = function(df, grupo = NULL, id = NULL){
  
  # Limpiar consola
  cat("\014")
  
  # Quitar el ID
  if(!is.null(id)){
    df[, id] <- NULL
  }
  
  # Tomar numéricos
  nums <- sapply(df, is.numeric)
  
  dfnum <- df[, nums]
  
  print(paste("Las dimensiones son ",nrow(df), "X", ncol(df) )) 
  print("Las variables relevadas son ")
  print( colnames(df))
  
  print("Sumario")
  print(summary(dfnum))
  print("")
  
  ## Tabla de contingencia GRUPOS
  if(! is.null(grupo)){
    freq<-table(grupo)
    prop<-round(prop.table(freq),2)
    print("Frecuencias del grupo")
    print(cbind(freq,prop))
    print("")
  }
  
  print("Estadística univariada")
  print(describeBy(dfnum, group = grupo))
  print("")
  
  print("Coeficiente de variación")
  print(stat.desc(dfnum)[14, ])
  print("")
  
  par(mfrow = c(2,2))
  
  for(j in 1:ncol(dfnum)){
    col <- dfnum[, j]
    hist(col, main=paste(colnames(dfnum)[j], "- Histograma"), col=rainbow(1))
    plot(density(col), main=paste(colnames(dfnum)[j], "- Densidad"))
    qqnorm(col)
    boxplot(col, main=paste(colnames(dfnum)[j], "- Boxplot"), col=rainbow(10))
  }
  
  par(mfrow = c(1,1))
}

# Nombre      : DM_Asociacion
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_Asociacion = function(df, grupo = NULL, id = NULL){
  
  # Limpiar consola
  cat("\014") 
  
  # Quitar el ID
  if(!is.null(id)){
    df[, id] <- NULL
  }
  
  # Tomo solo numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  dfnum$Id <- NULL
  
  # Matrices 
  df.cov <- round(cov(dfnum), 2)
  df.cor <- round(cor(dfnum), 2)
  
  print("Matriz (S)")
  print(df.cov)
  print("")
  
  print("Determinante de R. Cuanta mayor asociación, más cerca de 0")
  print(det(df.cor))
  
  print("Matriz (R)")
  print(df.cor)
  print("")
  
  #pairs(dfnum)
  
  #Para cambiar la representacion grafica reordeno las columnas
  #df2<-as.data.frame(cbind(dfnum, grupo))
  print(ggcorr(dfnum))
  print(rcorr(as.matrix(dfnum)))
  print(ggpairs(dfnum,axisLabels="none") )
  #corrplot(df.cor, method = "square", tl.col = "black")
}

# Nombre      : DM_OutliersMV
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_OutliersMV = function(df, id = NULL, norm = TRUE){
  
  # Limpiar consola
  cat("\014") 
  
  # Quitar el ID
  if(!is.null(id)){
    df[, id] <- NULL
  }
  
  # Tomo solo numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  
  # Normalizar variables
  if(norm){
    dfnum <- scale(dfnum)
  }
  
  mahab <- mahalanobis(dfnum, colMeans(dfnum), cov(dfnum)) 
  print(mahab)
  plot(mahab)
  abline(qchisq(0.95, 4, lower.tail = T),0, col="red")
  
  print("Outliers detectados")
  print(which(mahab >(qchisq(0.95, 4, lower.tail = T))))
  
  mahd <- OutlierMahdist(dfnum)
  print(mahd)
  plot(mahd)
  
}

# Nombre      : DM_ClusterJerarquico
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_ClusterJerarquico = function(df, plotear = FALSE, kc = NULL, 
                                id = NULL, alt = FALSE){

  # Limpiar consola
  cat("\014") 
  
  names <- NULL
  
  # Quitar el ID
  if(!is.null(id)){
    names <- df[, id]
    df[, id] <- NULL
  }
  
  # Tomo solo numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  
  # Normalizar
  dfnum <- scale(dfnum)
  
  if(!is.null(names)){
    rownames(dfnum) <- names
  }
  
  #print("Euclidea - Matriz de disimilitud")
  zDist <- dist(as.matrix(dfnum), method = "euclidean")
  
  #print(zDist)
  
  #Cluster jerárquico
  c1 <- hclust(zDist, method = "complete", members = NULL)
  d1 <- cophenetic(c1)
  print(paste("Complete COPH: ", cor(d1,zDist)))
  
  c2<-hclust(zDist, method = "single", members = NULL)
  d1<-cophenetic(c2)
  print(paste("Single COPH: ", cor(d1,zDist)))
  
  c3<-hclust(zDist, method = "centroid", members = NULL)
  d1<-cophenetic(c3)
  print(paste("Centroid COPH: ", cor(d1,zDist)))
  
  c4<-hclust(zDist, method = "average", members = NULL) #la mejor
  d1<-cophenetic(c4)
  print(paste("Average COPH: ", cor(d1,zDist)))
  
  c5<-hclust(zDist, method = "ward.D", members = NULL)
  d1<-cophenetic(c5)
  print(paste("Ward COPH: ", cor(d1,zDist)))
  
  #Graficamos
  if(plotear){
    plot(c1, main = "Complete")
    
    #Le asignamos a cada objeto el cluster de pertenencia
    if(!is.null(kc)){
      dfnum$grupo <-cutree(c1, k= kc)
      #Marcamos los clusters en el dendrograma
      rect.hclust(c1, k= kc, border=2:5)
    }
    
    plot(c2, main = "Single")
    
    #Le asignamos a cada objeto el cluster de pertenencia
    if(!is.null(kc)){
      dfnum$grupo <-cutree(c2, k= kc)
      #Marcamos los clusters en el dendrograma
      rect.hclust(c2, k= kc, border=2:5)
    }
    
    plot(c3, main = "Centroid")
    
    #Le asignamos a cada objeto el cluster de pertenencia
    if(!is.null(kc)){
      dfnum$grupo <-cutree(c3, k= kc)
      #Marcamos los clusters en el dendrograma
      rect.hclust(c3, k= kc, border=2:5)
    }
    
    plot(c4, main = "Average")
    
    #Le asignamos a cada objeto el cluster de pertenencia
    if(!is.null(kc)){
      dfnum$grupo <-cutree(c4, k= kc)
      #Marcamos los clusters en el dendrograma
      rect.hclust(c4, k= kc, border=2:5)
    }
    
    plot(c5, main = "Ward") 
    
    #Le asignamos a cada objeto el cluster de pertenencia
    if(!is.null(kc)){
      dfnum$grupo <-cutree(c5, k= kc)
      #Marcamos los clusters en el dendrograma
      rect.hclust(c5, k= kc, border=2:5)
    }
    
  }
  else {
    print("Número óptimo de clusters --")
    print(fviz_nbclust(dfnum, kmeans, method = "wss"))
  }
  
  # Otra opción
  if(!is.null(kc) & alt){
    #rownames(dfnum) <- NULL
    #Análisis de cluster jerarquico e identificacion de clusters
    clu <- hcut(zDist, 
                k = kc, stand = TRUE, hc_method = "ward.D", hc_metric = "euclidean")
    fviz_dend(clu, rect = TRUE, cex = 0.5)
  }
}

# Nombre      : DM_ClusterKMeans
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_ClusterKMeans = function(df, id = NULL, kc){
  
  # Limpiar consola
  cat("\014") 
  
  names <- NULL
  
  # Quitar el ID
  if(!is.null(id)){
    names <- df[, id]
    df[, id] <- NULL
  }
  
  # Tomo solo numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  
  # Normalizar
  z <- scale(dfnum)
  
  # Nombres de filas
  if(!is.null(names)){
    rownames(dfnum) <- names
  }
  
  # Análisis de cluster por K-Means
  c<-kmeans(z, centers=kc, iter.max = 100 )
  
  #se agrega el cluster de pertenencia a cada objeto
  dfnum$grupo <-c$cluster
  
  print("Tamaño de los Clusters")
  print(c$size)
  print("")
  print("Medias de los Clusters")
  print(round(c$centers, 2))
  #grafico
  plot(z, col = c$cluster, pch = 19, frame = FALSE,
       main = paste("K-means con k =", kc))
  points(c$centers, col = 1:4, pch = 8, cex = 2)
  
  describe.by(dfnum, dfnum$grupo)
}

# Nombre      : DM_ACP
# Creado por  : DM
# Descripción : 
# Ejemplo     : 
DM_ACP = function(df, grupo = NULL, id = NULL, alt = FALSE, norm = FALSE){
  # Limpiar consola
  cat("\014") 
  
  names <- NULL
  
  # Quitar el ID
  if(!is.null(id)){
    names <- df[, id]
    df[, id] <- NULL
  }
  
  # Tomo solo numéricos
  nums <- sapply(df, is.numeric) 
  dfnum <- df[, nums]
  
  if(norm){
    dfnum <- scale(dfnum)
  }
  
  # Nombres de filas
  if(!is.null(names)){
    rownames(dfnum) <- names
  }
  
  if(!alt){
    ACP <- rda(dfnum, scale=TRUE)
    print(ACP)
    print(summary(ACP)) #Por defecto los scores estan escalados
    
    #Graficos de sedimentación
    screeplot(ACP)
    screeplot(ACP, type="l")
    evplot (ACP$CA$eig)
    
    print("Autovectores: Coeficientes de las combinaciones lineales")
    coef.species<- scores(ACP, choices = 1:2, display = "species", scaling = 0)
    print(round(coef.species, 2))
    
    print("Correlaciones de las variables originales con los componentes")
    #Calculamos los scores de los sitios
    coef.sites<-scores(ACP, choices = 1:2, display = "sites", scaling = 0) 
    print(round(coef.sites, 2))
    
    #Biplots Ejes 1 y 2
    par(mfrow = c(1, 2))
    biplot(ACP, cex=0.5, scaling=1, main="Escalado 1:\n Distancias entre los sitios") #Escalado 1 (para biplotde distancias): Preserva las distancias entre los sitios, pero no la asociación entre variables
    biplot(ACP, cex=0.5, scaling=2, main="Escalado 2:\n Correlaciones entre las variables") # Escalado 2 (para biplotde correlaciones): Preserva las asociaciones entre las variables, pero las distancias entre sitios pueden estar distorsionadas
    par(mfrow = c(1, 1))
  }
  else {
    # Otra función
    acp1<-PCA(dfnum, graph = TRUE)
    print(acp1$eig)
    # Extraer autovalores/variances 
    get_eig(acp1)
    # Visualize autovalores/variances 
    fviz_eig(acp1, addlabels=TRUE, hjust = -0.3)+ theme_minimal()
    print("Extraer los resultados de las variables")
    var <- get_pca_var(acp1) 
    print(var$cor)
    
    # Graph of variables: default plot 
    print(fviz_pca_var(acp1, col.var = "steelblue"))
    print(fviz_pca_biplot(acp1))
    print(fviz_pca_biplot(acp1,axes = c(1, 3)))
    
    print("Coordenadas")
    print(acp1$ind$coord)
    print(acp1$var$contrib)
  }
}

evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}