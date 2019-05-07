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
fly.non_bound$Unido <- "No" 
fly.non_bound <- formatEnsembl(fly.non_bound)
fly.bound$Unido <- "Si"
fly.bound <- formatEnsembl(fly.bound)
fly.random$Unido <- "?"
fly.random <- formatEnsembl(fly.random)

glimpse(fly.bound)
glimpse(fly.non_bound)
glimpse(fly.random)

fly <- rbind(fly.bound, fly.non_bound)
fly$Unido <- factor(fly$Unido)

# 2. Agregar columnas...
nv1 <- fly %>% select(	Chromosome, 	LoopPattern, 
								TerminalPair,	N.2, N.1, N2,
								N5, N6, N7, N8, 
								Pairments, WooblePairs,
								Bulges, InternalLoops, 
								AU_PercentPairs, PurinePercentStem,
								RnaFoldMFE, RelativePosition, 
								CG_PercentPairs, PurinePercentStem,
								Unido)

# Training model #########################
set.seed(123)
split <- initial_split(nv1, prop = .7, strata = "Unido")
train <- training(split)
test  <- testing(split)

table(train$Unido) %>% prop.table()
table(test$Unido) %>% prop.table()

install.packages("caret")
library(caret)

# create response and feature data
features <- setdiff(names(train), "Unido")
x <- train[, features]
y <- train$Unido

# set up 10-fold cross validation procedure
train_control <- trainControl(
  method = "cv", 
  number = 10
)

# set up tuning grid
search_grid <- expand.grid(
  usekernel = c(TRUE, FALSE),
  fL = 0:5,
  adjust = seq(0, 5, by = 1)
)

# train model
nb.m2 <- train(
  x = x,
  y = y,
  method = "nb",
  trControl = train_control,
  tuneGrid = search_grid,
  preProc = c("BoxCox", "center", "scale", "pca")
)

y# top 5 modesl
nb.m2$results %>% 
  top_n(5, wt = Accuracy) %>%
  arrange(desc(Accuracy))


# plot search grid results
plot(nb.m2)

confusionMatrix(nb.m2)

pred <- predict(nb.m2, newdata = test)
confusionMatrix(pred, test$Unido)

#install.packages("e1071")
#library(e1071)

nBMod <- naiveBayes(Unido  ~  ., data = nv1, laplace = 0.01)

#Getting started with Naive Bayes in mlr
#Install the package
#install.packages("mlr")
#Loading the library
#library(mlr)

NB_pred <- predict(nBMod, nv1)
  
tb <- table(NB_pred, nv1$Unido)
sum(diag(tb))/sum(tb)

summary(nBMod)

##Confusion ############################

# # glm
# if(!require(pROC)){
#   install.packages("pROC")
#   library(pROC)
# }
# 
# setwd(paste(getwd(), "/data", sep=""))
# 
# 
# train <- sample(nrow(fly), nrow(fly)*0.75)
# 
# # Specify a null model with no predictors
# null_model <- glm(Unido ~ 1, data = fly[train,], family = "binomial")
# full_model <- glm(Unido  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5 + N6 + N7
#                   + Bulges + InternalLoops + CG_PercentPairs + N8
#                   + TerminalPair + WooblePairs + GU_PercentPairs + RelativePosition + Pairments, 
#                   data = fly[train,], family = "binomial", control=glm.control(maxit=100))
# # Use a forward stepwise algorithm to build a parsimonious model
# step_model <- step(null_model, scope = list(lower = null_model, upper = full_model), 
#                    direction = "forward")
# 
# mean(laslo[train, ]$Serie)
# laslo$Serie_prob <- predict(step_model, laslo, type="response")
# laslo$Serie_pred <- ifelse(laslo$Serie_prob > 0.006962785, 1, 0)
# 
# table(laslo[-train, ]$Serie, laslo[-train, ]$Serie_pred)
# ROC <- roc(laslo$Serie, laslo$Serie_prob)
# # Plot the ROC curve
# plot(ROC, col = "blue")
# # Calculate the area under the curve (AUC)
# auc(ROC)
# 
# mean(mito$Serie)
# mito$Serie_prob <- predict(step_model, mito, type="response")
# mito$Serie_pred <- ifelse(mito$Serie_prob > 0.0077, 1, 0)
# 
# table(mito$Serie, mito$Serie_pred)
# ROC <- roc(mito$Serie, mito$Serie_prob)
# # Plot the ROC curve
# plot(ROC, col = "blue")
# # Calculate the area under the curve (AUC)
# auc(ROC)



######################
# Load the rpart package
library(rpart)

# Build a lending model predicting loan outcome versus loan amount and credit score
loan_model <- rpart(outcome ~ loan_amount + credit_score, 
                    data = loans, method = "class", control = rpart.control(cp = 0))

# Make a prediction for someone with good credit
predict(loan_model, good_credit, type = "class")

# Load the rpart.plot package
library(rpart.plot)

# Plot the loan_model with default settings
rpart.plot(loan_model)

# Plot the loan_model with customized settings
rpart.plot(loan_model, type = 3, box.palette = c("red", "green"), fallen.leaves = TRUE)

# Determine the number of rows for training
nrow(loans)*0.75

# Create a random sample of row IDs
sample_rows <- sample(nrow(loans), 8484)

# Create the training dataset
loans_train <- loans[sample_rows, ]

# Create the test dataset
loans_test <- loans[-sample_rows, ]

# Grow a tree using all of the available applicant data
loan_model <- rpart(outcome ~ ., data = loans, method = "class", control = rpart.control(cp = 0))

# Make predictions on the test dataset
loans_test$pred <- predict(loan_model, loans_test, type="class")

# Examine the confusion matrix
table(loans_test$outcome, loans_test$pred)

# Compute the accuracy on the test dataset
mean(loans$outcome == loans_test$pred)

##################################################
# Load the randomForest package
if(!require(randomForest)){
  install.packages("randomForest")
  library(randomForest)
}

# Build a random forest model
loan_model <- randomForest(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5 + N6 + N7
                           + Bulges + InternalLoops + CG_PercentPairs + N8
                           + TerminalPair + WooblePairs + GU_PercentPairs, 
                           data = laslo[train,])

# Compute the accuracy of the random forest
loans_test$pred <- predict(loan_model, loans_test, type="class")
mean(loans_test$pred == loans_test$outcome)