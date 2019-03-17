# Cargar framework
source("./src/configuration.R")

# Cargar datos
fly.bound <- read.csv("./data/non_bound_fly/smaug_bound_fruitfly_8.csv",
                      sep =";", 
                      dec =",", stringsAsFactors = TRUE)

fly.non_bound <- read.csv("./data/non_bound_fly/non_bound_smaug_8.csv",
                          sep =";", 
                          dec =",", stringsAsFactors = TRUE)

fly.non_bound$Unido <- "No" 
fly.non_bound$GenId <- fly.non_bound$Column5
fly.non_bound$TranscriptoId <- fly.non_bound$Column2
fly.non_bound$GenSymbol <- fly.non_bound$Column4
fly.non_bound$Column1 <- NULL
fly.non_bound$Column2 <- NULL
fly.non_bound$Column3 <- NULL
fly.non_bound$Column4 <- NULL
fly.non_bound$Column5 <- NULL
fly.non_bound$Column6 <- NULL

#Eliminar variables innecesarias
fly.bound$Unido <- "Si"
fly.bound <- formatEnsembl(fly.bound)

glimpse(fly.bound)

# 2. Agregar columnas...
bnd1 <- non_bound %>% select(	Chromosome, 	LoopPattern, 
								TerminalPair,	N.2, N.1, N2,
								N5, N6, N7, N8, 
								Pairments, WooblePairs,
								Bulges, InternalLoops, 
								AU_PercentPairs, PurinePercentStem,
								RnaFoldMFE, RelativePosition, 
								Additional5Seq, Additional3Seq,
								CG_PercentPairs, PurinePercentStem,
								Unido)

df <- rbind(bnd, bnd1)
df$Unido <- factor(df$Unido)

install.packages("e1071")
library(e1071)
nBMod <- naiveBayes(Unido  ~  ., data = df, laplace = 0.01)

#Getting started with Naive Bayes in mlr
#Install the package
install.packages("mlr")
#Loading the library
library(mlr)

NB_pred <- predict(nBMod, df)
  
tb <- table(NB_pred, df$Unido)
sum(diag(tb))/sum(tb)

test <- read.csv("./data/non_bound_fly/Mig6_UTR_8.csv",
                 sep =";", 
                 dec =",", stringsAsFactors = TRUE)
test$Unido <- NA
test <- test %>% select(SequenceID, LoopPattern, TerminalPair,N.2, N.1, N2,
                             N5, N6, N7, N8, Pairments, WooblePairs,
                             Bulges, InternalLoops, AU_PercentPairs, PurinePercentStem,
                             RnaFoldMFE, RelativePosition, Unido)


test$Unido <- predict(nBMod, test)
test

mouse <- read.csv("./data/non_bound_fly/mouse_nizou_08_n6_nlp_t25.csv",
                  sep =";", 
                  dec =",", stringsAsFactors = TRUE)
mouse$Unido <- NA
mouse <- mouse %>% select(Column5, LoopPattern, TerminalPair,N.2, N.1, N2,
                        N5, N6, N7, N8, Pairments, WooblePairs,
                        Bulges, InternalLoops, AU_PercentPairs, PurinePercentStem,
                        RnaFoldMFE, RelativePosition, Unido)

mouse$Unido <- predict(nBMod, mouse)
head(select(mouse, Column5, LoopPattern, N.1, N2, RnaFoldMFE, Unido), 20)

#Create a classification task for learning on Titanic Dataset and specify the target feature
task = makeClassifTask(data = df, target = "Unido")

#Initialize the Naive Bayes classifier
selected_model = makeLearner("classif.naiveBayes")

#Train the model
NB_mlr = train(selected_model, task)

#Read the model learned  
NB_mlr$learner.model

#Predict on the dataset without passing the target feature
predictions_mlr = as.data.frame(predict(NB_mlr, newdata = Titanic_dataset[,1:3]))

##Confusion


# glm
if(!require(pROC)){
  install.packages("pROC")
  library(pROC)
}

setwd(paste(getwd(), "/data", sep=""))

laslo <- loadFiles(path="./mouse/", pattern = "*.csv")
mito <- loadFiles(path="./mitocondriales/", pattern = "*.csv")

train <- sample(nrow(laslo), nrow(laslo)*0.75)

# Specify a null model with no predictors
null_model <- glm(Serie ~ 1, data = laslo[train,], family = "binomial")
full_model <- glm(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5 + N6 + N7
                  + Bulges + InternalLoops + CG_PercentPairs + N8
                  + TerminalPair + WooblePairs + GU_PercentPairs, 
                  data = laslo[train,], family = "binomial")
# Use a forward stepwise algorithm to build a parsimonious model
step_model <- step(null_model, scope = list(lower = null_model, upper = full_model), 
                   direction = "forward")

mean(laslo[train, ]$Serie)
laslo$Serie_prob <- predict(step_model, laslo, type="response")
laslo$Serie_pred <- ifelse(laslo$Serie_prob > 0.006962785, 1, 0)

table(laslo[-train, ]$Serie, laslo[-train, ]$Serie_pred)
ROC <- roc(laslo$Serie, laslo$Serie_prob)
# Plot the ROC curve
plot(ROC, col = "blue")
# Calculate the area under the curve (AUC)
auc(ROC)

mean(mito$Serie)
mito$Serie_prob <- predict(step_model, mito, type="response")
mito$Serie_pred <- ifelse(mito$Serie_prob > 0.0077, 1, 0)

table(mito$Serie, mito$Serie_pred)
ROC <- roc(mito$Serie, mito$Serie_prob)
# Plot the ROC curve
plot(ROC, col = "blue")
# Calculate the area under the curve (AUC)
auc(ROC)



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