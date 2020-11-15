library(dplyr)
library(ggplot2)
library(tidyr)

convert.z.score <- function(z, alternative = "one.sided") {
  if (alternative == "one.sided") {
    pval <- pnorm(-abs(z))
    pval <- 2 * pval
  } else if (alternative == "less") {
    pval <- pnorm(z)
  } else {
    pval <- pnorm(-z)
  }
  return(pval)
}

get.p.value <- function(df, dfRand, n, m, colName, keyName, alternative = "one.sided") {
  xDist <- vector(mode = "numeric", length = n)
  # Obtener frecuencia observada en la muestra original
  num_exp <- sum(df[, which(colnames(df) == colName)] == keyName)

  withReplace <- (m > nrow(dfRand))

  # Tomar N muestras de tamaño M con repetición
  for (i in 1:n) {
    aux <- sample_n(dfRand, m, replace = withReplace)

    # dfRand[sample(nrow(dfRand), m, replace = TRUE), ]

    # Contar la frecuencia de la clave buscada
    num <- sum(aux[, which(colnames(aux) == colName)] == keyName)

    # num <- distinct(filter(dfRand,Loop==keyName), Transcripto) %>% nrow()

    # Generar una distribucion con las frecuencias obtenidas
    xDist[i] <- num
  }

  if (sum(xDist) == 0) {
    return(0.0)
  }


  # Normalizar
  media <- mean(xDist)
  ee <- sd(xDist)

  xDist <- scale(xDist)

  # Es normal ??
  st <- shapiro.test(xDist)

  if (sum(xDist) == 0) {
    return(-2)
  }

  if (st$p.value < 0.05) {
    # print(keyName)
    return(-1)
  }

  z <- (num_exp - media) / (ee * sqrt(length(xDist)))
  return(convert.z.score(z, alternative))
}

plotPDens <- function(df, dfRand, n, m, colName, keyName, type = NULL) {
  xDist <- vector(mode = "numeric", length = n)

  # Obtener frecuencia observada en la muestra original
  num_exp <- sum(df[, which(colnames(df) == colName)] == keyName)

  # Tomar N muestras de tamaño M con repetición
  for (i in 1:n) {
    aux <- dfRand[sample(nrow(dfRand), m, replace = TRUE), ]
    # Contar la frecuencia de la clave buscada
    num <- sum(aux[, which(colnames(aux) == colName)] == keyName)
    # Generar una distribucion con las frecuencias obtenidas
    xDist[i] <- num
  }

  # Normalizar
  media <- mean(xDist)
  ee <- sd(xDist)

  if (media == 0) {
    return(0.0)
  }

  xDist <- scale(xDist)
  num_exp <- (num_exp - media) / (ee * sqrt(length(xDist)))
  dens <- density(xDist)

  # Punto dentro de la distribución para el valor observado
  xe <- approx(dens$x, dens$y, xout = num_exp)

  # print("Frecuencia esperada")
  # print(num_exp)

  # print("Resumen de la distribucion x")
  # print(summary(xDist))

  # print("p.value")
  # print(convert.z.score(num_exp, type))

  # Grafico
  # p <- plot(dens, main = paste(keyName, "fd(x)", "(Tama�o:", length(xDist), ")",sep=" "))

  # dd <- with(dens,data.frame(x,y))
  # pl <- qplot(x,y,data=dd,geom="line")+
  #   geom_ribbon(data=subset(dd,x>xe ),aes(ymax=y),ymin=0,
  #               fill="red",colour=NA,alpha=0.5)
  # print(pl)
  # Func de interpolacion
  # f <- approxfun(dens)

  # Valor dentro del área
  # pvalue <- integrate(f, lower=min(xDist), upper=xe$x, stop.on.error = TRUE)

  # print("Integracion")
  # print(pvalue)
  # print(pvalue$message)
  #
  # dF <- function(x)dnorm(x)
  # pF <- function(q)integrate(dF,-Inf,q) #$value
  # print(pF(55))

  return(convert.z.score(num_exp, type))
}

################################################################################################


# Muestra original
df.exp <- fly_bound[fly_bound$SerieID == 1, ]
# df.exp$LoopId <- paste(df.exp$Loop, df.exp$TerminalPair)

# Muestra randomizada
df.rnd <- fly_bound[fly_bound$SerieID == 2, ]
# df.rnd$LoopId <- paste(df.rnd$Loop, df.rnd$TerminalPair)

# Tomar 50 muestras de n
x <- plotPDens(df.exp, df.rnd, 100, nrow(df.exp), "Loop", "CUGGUU")

tb <- tb %>% filter(`2` > 1)
tb$p.two.sided <- 0
tb$p.less <- 0
tb$p.greater <- 0

for (i in 1:nrow(tb)) {
  aux <- tb[i, ]$Var1

  tb[i, ]$p.two.sided <-
    get.p.value(df.exp, df.rnd, 50, nrow(df.exp), "Loop", aux)
  tb[i, ]$p.less <-
    get.p.value(df.exp, df.rnd, 50, nrow(df.exp), "Loop", aux, "less")
  tb[i, ]$p.greater <-
    get.p.value(df.exp, df.rnd, 50, nrow(df.exp), "Loop", aux, "greater")
}

tb %>%
  filter(p.greater < 0.05) %>%
  arrange(p.greater)
tb %>%
  filter(p.two.sided < 0.05) %>%
  arrange(p.two.sided)
tb %>%
  filter(p.less < 0.05) %>%
  arrange(p.less)
