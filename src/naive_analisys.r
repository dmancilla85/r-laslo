# Cargar framework
source("./src/configuration.R")

# Cargar datos
fly.bound <- read.csv("./data/stem15/fly_bound.csv",
                      sep =";", 
                      dec =",", stringsAsFactors = TRUE)

fly.non_bound <- read.csv("./data/stem15/fly_non_bound.csv",
                          sep =";", 
                          dec =",", stringsAsFactors = TRUE)

# First random file
fly.random <- read.csv("./data/stem15/fly_bound_rnd1.csv",
                          sep =";", 
                          dec =",", stringsAsFactors = TRUE)
# Second random file
fly.random <- rbind(fly.random, 
                    read.csv("./data/stem15/fly_bound_rnd1.csv",
                              sep =";", 
                              dec =",", stringsAsFactors = TRUE))


#Eliminar variables innecesarias
fly.non_bound$Tipo <- "No unidos" 
fly.non_bound <- formatEnsembl(fly.non_bound, FALSE)
fly.bound$Tipo <- "Unidos"
fly.bound <- formatEnsembl(fly.bound, FALSE)
fly.random$Tipo <- "Shuffled"
fly.random <- formatEnsembl(fly.random, FALSE)

glimpse(fly.bound)
glimpse(fly.non_bound)
glimpse(fly.random)

fly <- rbind(fly.bound, fly.non_bound)
fly <- rbind(fly, fly.random)
fly$Tipo <- factor(fly$Tipo, levels=c("Unidos", "No unidos", "Shuffled"))

# Análisis de las variables

# "GenID"              
# "TranscriptoID"      
# "GenSymbol"

# "Chromosome"
fly %>% filter(Tipo!="Shuffled") %>%
  ggplot(aes(x=Chromosome, fill=Tipo)) +
  geom_bar(color = "darkgray", position="fill", width = 0.6) + theme_gray() +
  ggtitle("Distribución de los dos conjuntos en cromosomas") +
  scale_fill_manual(values=c('blue2','brown3')) + xlab("Cromosoma") + ylab("Recuento") #480x480px


# "LoopPattern" 
fly %>%
  ggplot(aes(fill=LoopPattern, x=1)) +
  geom_bar(position="fill", color="darkgrey", ) + coord_polar(theta="y") +
  geom_text(stat = "fill_labels", fontface="bold",
            position = position_dodge(width = .8), size=3.5, 
            check_overlap = TRUE, ) + labs(fill="Secuencia del loop") +
  facet_wrap(~ Tipo) + ylab("") +
  ggtitle("Motivos hallados en el loop", 
          subtitle="Sobre todos los stem-loops hallados en los transcriptos") +
  theme(axis.text.x=element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) + ylab("") + xlab("") +
  scale_fill_brewer(palette="Spectral")

# "TerminalPair"
#install.packages("wesanderson")
library(wesanderson)



# "N.2" "N.1" "N2"  "N5"  "N6"  "N7"  "N8"
# "Loop"              
# "Pairments"          "WooblePairs"        "Bulges"             "InternalLoops"      
# "SequenceLength"

# "A_PercentSequence"  "C_PercentSequence" "G_PercentSequence"  "U_PercentSequence"  
fly %>% filter(Tipo!="Shuffled") %>%
  select(Tipo, A_PercentSequence, C_PercentSequence, 
         G_PercentSequence, U_PercentSequence) %>% 
  gather(metric, value, -Tipo) %>%
  mutate(metric=paste ("", str_replace(metric, "_PercentSequence", "")), value = value * 100) %>%
  ggplot(aes(x=value, fill=Tipo)) + ylab("Densidad") +
  geom_histogram(aes(y=..density..), position="identity", bins=40, alpha=0.7, color="black") + 
  xlab("Ratio de cada nucleótido en la secuencia entera") +
  facet_wrap(~ metric,ncol=2, nrow=2 ) + theme_pubr() +
  ggtitle("Composición de las secuencias cDNA", 
          subtitle = "Sobre el total de bases en transcriptos") +
  scale_fill_manual(values=c("blue","red3"))


# "AU_PercentPairs"    "CG_PercentPairs"    "GU_PercentPairs" 
fly %>% filter(Tipo!="Shuffled") %>%
  select(Tipo, AU_PercentPairs, CG_PercentPairs, GU_PercentPairs) %>% 
  gather(metric, value, -Tipo) %>%
  mutate(metric=paste ("", str_replace(metric, "_PercentPairs", "")), value = value * 100) %>%
  ggplot(aes(x=value, fill=Tipo)) + ylab("Densidad") +
  geom_density(alpha=0.6) +
  xlab("Ratio de cada apareamiento en los stem-loop") +
  facet_wrap(~ metric,ncol=3, nrow=1 ) + theme_pubr() +
  ggtitle("Composición de apareamientos en los stem-loops", 
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  scale_fill_manual(values=c("yellow2","magenta4"))

# "PurinePercentPairs" 
fly %>% 
  select(Tipo, PurinePercentPairs) %>% 
  ggplot(aes(x=PurinePercentPairs, group=Tipo)) + ylab("Densidad") +
  geom_density(alpha=1, size=1) +
  facet_wrap(~ Tipo) +
  xlab("Ratio de A-G en apareamientos de los stems") +
  theme_pubr() + 
  ggtitle("A-G (purinas) en los stem-loops", 
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  geom_vline(data=fly, aes(xintercept=mean(PurinePercentPairs), color="blue"), 
             linetype="dashed", show.legend = F, size=1) 

# "RnaFoldMFE"  

# fly %>%
#   group_by(Tipo) %>%
#   dplyr::summarize(Mean = mean(RnaFoldMFE, na.rm=TRUE))

fly %>% 
  select(Tipo, RnaFoldMFE) %>% 
  ggplot(aes(x=RnaFoldMFE, group=Tipo)) + ylab("Densidad") +
  geom_density(alpha=1, size=1) +
  facet_wrap(~ Tipo) +
  xlab("Mínima energía libre (kcal/mol)") +
  theme_pubr() + 
  ggtitle("Energía libre en los stem-loops", 
          subtitle = "Sobre el total de stem-loops en todos los transcriptos") +
  geom_vline(data=fly, aes(xintercept=mean(RnaFoldMFE), color="red"), 
             linetype="dashed", show.legend = F, size=1) 

# "RelativePosition"   
# "Unido"    


# Modelo NAIVE BAYES ############################################
# 1. Seleccionar columnas
nv1 <- fly %>% select(	Chromosome, 	LoopPattern, 
								TerminalPair,	N.2, N.1, N2,
								N5, N6, N7, N8, 
								Pairments, WooblePairs,
								Bulges, InternalLoops,RelativePosition,
								AU_PercentPairs, GU_PercentPairs,
								RnaFoldMFE, RelativePosition, 
								CG_PercentPairs, PurinePercentPairs,Unido)

# 2. Separar sets
set.seed(123)
split <- initial_split(nv1, prop = .7, strata = "Unido")
train <- training(split)
test  <- testing(split)

table(train$Unido) %>% prop.table()
table(test$Unido) %>% prop.table()

#install.packages("caret")
#library(caret)

# 3. Create response and feature data
features <- setdiff(names(train), "Unido")
x <- train[, features]
y <- train$Unido

# 4. set up 10-fold cross validation procedure
train_control <- trainControl(
  method = "cv", 
  number = 10
)

# 5. set up tuning grid
search_grid <- expand.grid(
  usekernel = c(TRUE, FALSE),
  fL = 0:5,
  adjust = seq(0, 5, by = 1)
)

# 6. train model
nb.m2 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control,
  tuneGrid = search_grid,
  preProc = c("BoxCox", "center", "scale", "pca")
)

# 7. ver top 5 models
nb.m2$results %>% 
  top_n(5, wt = Accuracy) %>%
  arrange(desc(Accuracy))


# 8. plot search grid results
plot(nb.m2)

confusionMatrix(nb.m2)

pred <- predict(nb.m2, newdata = test)
confusionMatrix(pred, test$Unido)

mean(nv1$Unido=="Si") # 0.6596544

nv1$Serie_prob <- predict(nb.m2, nv1, type="prob")
nv1$Serie_pred <- ifelse(nv1$Serie_prob > 0.6594737, 1, 0)

ROC <- roc(nv1$Unido, nv1$Serie_prob$Si)

# # Plot the ROC curve
plot(main=paste("ROC - área bajo la curva:",round(auc(ROC),3)),ROC, col = "blue")
# # Calculate the area under the curve (AUC)



######################
# Load the rpart package
# library(rpart)
# 
# # Build a lending model predicting loan outcome versus loan amount and credit score
# loan_model <- rpart(outcome ~ loan_amount + credit_score, 
#                     data = loans, method = "class", control = rpart.control(cp = 0))
# 
# # Make a prediction for someone with good credit
# predict(loan_model, good_credit, type = "class")
# 
# # Load the rpart.plot package
# library(rpart.plot)
# 
# # Plot the loan_model with default settings
# rpart.plot(loan_model)
# 
# # Plot the loan_model with customized settings
# rpart.plot(loan_model, type = 3, box.palette = c("red", "green"), fallen.leaves = TRUE)
# 
# # Determine the number of rows for training
# nrow(loans)*0.75
# 
# # Create a random sample of row IDs
# sample_rows <- sample(nrow(loans), 8484)
# 
# # Create the training dataset
# loans_train <- loans[sample_rows, ]
# 
# # Create the test dataset
# loans_test <- loans[-sample_rows, ]
# 
# # Grow a tree using all of the available applicant data
# loan_model <- rpart(outcome ~ ., data = loans, method = "class", control = rpart.control(cp = 0))
# 
# # Make predictions on the test dataset
# loans_test$pred <- predict(loan_model, loans_test, type="class")
# 
# # Examine the confusion matrix
# table(loans_test$outcome, loans_test$pred)
# 
# # Compute the accuracy on the test dataset
# mean(loans$outcome == loans_test$pred)
# 
# ##################################################
# # Load the randomForest package
# if(!require(randomForest)){
#   install.packages("randomForest")
#   library(randomForest)
# }
# 
# # Build a random forest model
# loan_model <- randomForest(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5 + N6 + N7
#                            + Bulges + InternalLoops + CG_PercentPairs + N8
#                            + TerminalPair + WooblePairs + GU_PercentPairs, 
#                            data = laslo[train,])
# 
# # Compute the accuracy of the random forest
# loans_test$pred <- predict(loan_model, loans_test, type="class")
# mean(loans_test$pred == loans_test$outcome)