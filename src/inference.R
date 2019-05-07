library(infer)
#install.packages("broom")
library(broom)
library(ggplot2)
library(dplyr)

setwd("./data")

bound <- read.csv2("./stem15/fly_bound.csv")
bound$unido <- "Si"

non_bound <- read.csv2("./stem15/fly_non_bound.csv")
non_bound$unido <- "No"

fly <- rbind(bound, non_bound)

remove(bound)
remove(non_bound)

glimpse(fly)
fly$unido <- factor(fly$unido)

# Perform 1000 permutations
data_perm <- fly %>%
  # Specify fórmula, con "Unido" como éxito
  specify(unido ~ RnaFoldMFE, success = "Si") %>%
  # Hipotesis nula de independencia
  hypothesize(null="independence") %>% 
  # Generar 1000 repeticiones (por permutacion)
  generate(reps = 500, type = "permute") %>% 
  # Calculate the difference in proportions (male then female)
  calculate(stat="median", order=c("Si","No"))

# plot1
# Density plot of 1000 permuted differences in proportions
ggplot(data_perm, aes(x = stat)) + 
  geom_density()

# plot2
# Plot permuted differences, diff_perm
  ggplot(data_perm, aes(x = diff_perm)) + 
  # Add a density layer
  geom_density() +
  # Add a vline layer with intercept diff_orig
  geom_vline(aes(xintercept = diff_orig), color = "red")