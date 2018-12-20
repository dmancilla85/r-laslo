# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("R4RNA", version = "3.8")
# 
# library(R4RNA)
# library(Biostrings)
# 
# orig <- laslo[laslo$Serie=="Original", ]
# 
# for(i in 1:nrow(orig)){
#   line <- paste(">",i, laslo[i,]$Gen)
#   write(line,file="fasta.txt",append=TRUE)
#   write(laslo[i,]$StemLoopSequence ,file="fasta.txt",append=TRUE)
# }
# 
# known_file <- system.file("extdata", "vienna.txt", package = "R4RNA")
# known <- readVienna(known_file)
# known <- expandHelix(known)
# 
# 
# message("Multiple sequence alignment of interest")
# fasta_file <- "fasta.txt" #system.file("extdata", "fasta.txt", package = "R4RNA")
# fasta <- as.character(readBStringSet(fasta_file))
# message("Plot covariance in alignment")
# plotCovariance(fasta, known, cex = 0.5)

###########################################
# glm
if(!require(pROC)){
  install.packages("pROC")
  library(pROC)
}

setwd(paste(getwd(), "/data", sep=""))

laslo <- loadFiles(path="./mouse/", pattern = "*.csv")
mito <- loadFiles(path="./mitocondriales/", pattern = "*.csv")

laslo_model <- glm(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5+ N6
                      + Bulges + InternalLoops + CG_PercentPairs 
                   + TerminalPair, 
                      data = laslo, family = "binomial")

# Specify a null model with no predictors
null_model <- glm(Serie ~ 1, data = laslo, family = "binomial")
full_model <- glm(Serie  ~ RnaFoldMFE + LoopPattern + N.1 + N2 + N.2 + N5+ N6
                  + Bulges + InternalLoops + CG_PercentPairs 
                  + TerminalPair, data = laslo, family = "binomial")
# Use a forward stepwise algorithm to build a parsimonious model
step_model <- step(null_model, scope = list(lower = null_model, upper = full_model), 
                   direction = "forward")

mean(laslo$Serie)
laslo$Serie_prob <- predict(step_model, laslo, type="response")
laslo$Serie_pred <- ifelse(laslo$Serie_prob > 0.006722689, 1, 0)

table(laslo$Serie, laslo$Serie_pred)
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
