###############################################################################
# Modelo NAIVE BAYES 
###############################################################################

if(!require(carets)){
  install.packages('carets', dependencies=TRUE, repos='http://cran.rstudio.com/')
  library(carets)
}

if(!require(pROC)){
  install.packages("pROC")
  library(pROC)
}

# 1. Seleccionar columnas
nv1 <- fly_gb %>% select(-Gen, -Loop, -SequenceLength,
                         -TerminalPair,
                         -GU_PercentPairs,-AU_PercentPairs, -A_PercentSequence,
                         -AdditionalSeqMatches, -AdditionalSeqPositions, -Sense,
                         -EndsAt, -StartsAt,-PredictedStructure,-ViennaBracketStr,
                         -StemLoopSequence,-CDS_Start,-CDS_End)
nv1$Tipo <- factor(nv1$Tipo)

glimpse(nv1)

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
  number = 20
)

# 5. set up tuning grid
search_grid <- expand.grid(
  usekernel = F,
  fL = 1:5,
  adjust = seq(0, 5, by = 1)
)

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

mean(nv1$Tipo=="Unidos") # 0.6596544

predict(nb.m2, nv1, type="prob")
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