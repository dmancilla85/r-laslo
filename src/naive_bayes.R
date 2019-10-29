# Cargar datos
fly.bound <- read.csv(
  "./data/train/fly_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

fly.non_bound <- read.csv(
  "./data/train/fly_non_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)



# Eliminar variables innecesarias
fly.non_bound$Tipo <- "No unidos"
#fly.non_bound <- formatEnsembl(fly.non_bound, FALSE)
fly.bound$Tipo <- "Unidos"
#fly.bound <- formatEnsembl(fly.bound, FALSE)
fly.non_bound$Gen <- as.character(fly.non_bound$Gen)
fly.bound$Gen <- as.character(fly.bound$Gen)

########################
# Cargar datos
test.bound <- read.csv(
  "./data/stem15/fly_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)
test.bound$Tipo <- "Unidos"

test.non_bound <- read.csv(
  "./data/stem15/fly_non_bound.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)
test.non_bound$Tipo <- "No unidos"

test <- rbind(test.bound, test.non_bound)

test <- test %>% select(
  -Gen,
  -Loop,
  -SequenceLength,
  -TerminalPair,
  -AdditionalSeqMatches,
  -AdditionalSeqPositions,
  -Sense,
  -EndsAt,
  -StartsAt,
  -PredictedStructure,
  -ViennaBracketStr,
  -StemLoopSequence,
  -CDS_Start,
  -CDS_End
)

test.mouse <- read.csv(
  "./data/stem15/mouse_nizou.csv",
  sep = ";",
  dec = ",",
  stringsAsFactors = TRUE
)

test.mouse <- test.mouse %>% select(
  -Gen,
  -Loop,
  -SequenceLength,
  -TerminalPair,
  -AdditionalSeqMatches,
  -AdditionalSeqPositions,
  -Sense,
  -EndsAt,
  -StartsAt,
  -PredictedStructure,
  -ViennaBracketStr,
  -StemLoopSequence,
  -CDS_Start,
  -CDS_End
)

########################

# check intersections
fly.bound %>% inner_join(fly.non_bound, by = c("Gen"))

glimpse(fly.bound)
glimpse(fly.non_bound)

# igualar datasets
fly.bound$GeneSynonym <- NULL
fly.bound$Note <- NULL
fly.bound$AccessionID <- NULL
fly.non_bound$GeneSynonym <- NULL
fly.non_bound$Note <- NULL
fly.non_bound$AccessionID <- NULL

fly <- rbind(fly.bound, fly.non_bound)


# Determinando niveles
fly$Tipo <- factor(fly$Tipo, levels = c("Unidos", "No unidos"))

glimpse(fly)

if (!require(gridExtra)) {
  install.packages("gridExtra")
  library(gridEXtra)
}

####################
#install.packages("BSDA")
library(BSDA)
hist(scale(fly.bound$Bulges))
qqnorm(fly.bound$AU_PercentPairs)
shapiro.test(fly.bound$AU_PercentPairs)
median(fly.bound$WooblePairs)
median(fly.non_bound$WooblePairs)

summary(fly.bound$GU_PercentPairs)
summary(fly.non_bound$GU_PercentPairs)

wilcox.test(
  fly.bound$WooblePairs,
  fly.non_bound$WooblePairs,
  paired = F,
  conf.int = T
)

z <-
  z.test(
    x = fly.bound$RnaFoldMFE,
    y = fly.non_bound$RnaFoldMFE,
    # Two samples with normal distribution
    alt = "two.sided",
    # Dos colas
    mu = 0,
    # H_0: mu_1 - mu_2 = 0
    sigma.x = sd(fly.bound$RnaFoldMFE),
    # desviación estándar m
    sigma.y = sd(fly.non_bound$RnaFoldMFE),
    # desviación estandar n
    conf.level = 0.95
  )          # IC: error alpha_a/2 = 0.01/2

z
####################


# Hay correlaciones??
cfly <- fly.bound %>%
  select(
    -Gen,
    -Loop,
    -SequenceLength,
    -TerminalPair,
    -AdditionalSeqMatches,
    -EndsAt,
    -CDS_Start,
    -CDS_End,
    -StartsAt,
    -StemLoopSequence,
    -GU_PercentPairs,
    -AU_PercentPairs,
    -A_PercentSequence,
    -Sense,
    -AdditionalSeqPositions,
    -PredictedStructure,
    -ViennaBracketStr,
    -LoopPattern
  )

cfly %>% getCorPearson()
cfly %>% getCorTau()

###############################################################################
# Modelo NAIVE BAYES
###############################################################################


if (!require(caret)) {
  install.packages('caret', dependencies = TRUE, repos = 'http://cran.rstudio.com/')
  library(caret)
}

if (!require(pROC)) {
  install.packages("pROC")
  library(pROC)
}

# 1. Seleccionar columnas
nv1 <- fly %>% select(
  -Gen,
  -Loop,
  -SequenceLength,
  -TerminalPair,
  -AdditionalSeqMatches,
  -AdditionalSeqPositions,
  -Sense,
  -EndsAt,
  -StartsAt,
  -PredictedStructure,
  -ViennaBracketStr,
  -StemLoopSequence,
  -CDS_Start,
  -CDS_End
)

summary(nv1)

nv1$Tipo <- factor(nv1$Tipo)

#Install the package
#install.packages("mlr")
#Loading the library
library(mlr)

#Create a classification task for learning on Titanic Dataset and specify the target feature
task = makeClassifTask(data = nv1, target = "Tipo")
#Initialize the Naive Bayes classifier
selected_model = makeLearner("classif.naiveBayes")
#Train the model
NB_mlr = train(selected_model, task)
NB_mlr$learner.model



#Predict on the dataset without passing the target feature
predictions_mlr = as.data.frame(predict(NB_mlr, newdata = select(nv1, -Tipo)))

##Confusion matrix to check accuracy
table(predictions_mlr[, 1], nv1$Tipo)
311 / (311 + 52)
77 / (47 + 77)

#Predict on the dataset without passing the target feature
predictions_mlr = as.data.frame(predict(NB_mlr, newdata = select(test, -Tipo)))

##Confusion matrix to check accuracy
table(predictions_mlr[, 1], test$Tipo)
928 / (928 + 460)
340 / (340 + 388)



############## otra forma
# 2. Separar sets
set.seed(123)
split <- initial_split(nv1, prop = .7, strata = "Tipo")
train <- training(split)
test  <- testing(split)

table(train$Tipo) %>% prop.table()
table(test$Tipo) %>% prop.table()

# 3. Create response and feature data
features <- setdiff(names(train), "Tipo")
x <- train[, features]
y <- train$Tipo

# 4. set up 10-fold cross validation procedure
train_control <- trainControl(
  method = "cv",
  number = 60,
  search = "random",
  verboseIter = TRUE
)

# 5. set up tuning grid
search_grid <- expand.grid(usekernel = T,
                           fL = 1:5,
                           adjust = seq(0, 5, by = 1))

# 6. train model
nb.m2 <- caret::train(
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
confusionMatrix(pred, test$Tipo)

mean(nv1$Tipo == "Unidos") # 0.6596544

nv1$serie_prob <- predict(nb.m2, nv1, type = "prob")
nv1$serie_pred <- 0
nv1$serie_pred <- ifelse(nv1$serie_prob$Unidos > 0.7453799, 1, 0)

ROC <- roc(nv1$Tipo, nv1$serie_prob$Unidos)

# # Plot the ROC curve
plot(main = paste("ROC - área bajo la curva:", round(auc(ROC), 3)),
     ROC,
     col = "blue")
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

###########################################
# glm
# if(!require(pROC)){
#   install.packages("pROC")
#   library(pROC)
# }
#
# setwd(paste(getwd(), "/data", sep=""))
#
# laslo <- loadFiles(path="./mouse/", pattern = "*.csv")
# mito <- loadFiles(path="./mitocondriales/", pattern = "*.csv")
#
# train <- sample(nrow(laslo), nrow(laslo)*0.75)
#
# # Specify a null model with no predictors
# null_model <- glm(Serie ~ 1, data = laslo[train,], family = "binomial")
# full_model <- glm(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5 + N6 + N7
#                   + Bulges + InternalLoops + CG_PercentPairs + N8
#                   + TerminalPair + WooblePairs + GU_PercentPairs,
#                   data = laslo[train,], family = "binomial")
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