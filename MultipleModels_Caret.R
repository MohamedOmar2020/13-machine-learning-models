########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: Creating several machine learning models for bladder cancer progression
       ## Then Comparing the performance of these models 

###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Multiple_Caret")

## Load necessary packages
library(caret)
library(nnet)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)

#######################################################################

## Load data
load("./Objs/ProgressionDataGood.rda")

#load("/Users/mohamedomar/Documents/Research/Projects/RF_Classifier/Objs/RF_importances_TF_MiR.rda")
#rf_importances <- as.data.frame(rf_importances)
#RF_topGenes <- rf_importances[order(rf_importances$MeanDecreaseGini ,decreasing = TRUE), ]
#RF_topGenes <- rownames(RF_topGenes)

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat)
usedTestMat <- normalizeBetweenArrays(mixTestMat)

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

## Filter out any variables (genes) that are not expressed or do not have enough variance to be informative in classification. 
## We will first take the values and un-log2 them, then filter out any genes according to following criteria: (1) At least 20% of samples should have raw intensity greater than 100; (2) The coefficient of variation (sd/mean) is between 0.7 and 10.
X <- usedTrainMat
ffun <- filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))

filt <- genefilter(2^X,ffun)
usedTrainMat <- usedTrainMat[filt, ]

## The same for usedTestMat
usedTestMat <- usedTestMat[filt, ]

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

#################################################################
## Oversampling of the training data set to compensate for the un-balanced classes
set.seed(333)
Data_train <- as.data.frame(Data_train)
Data_train[,2914] <- as.factor(Data_train[,2914])
Over_Train <- SMOTE(usedTrainGroup~., data = Data_train, perc.over = 300, perc.under = 134)
table(Over_Train[,2914])

######
control <- trainControl(method="repeatedcv", number=10, repeats=3, classProbs = TRUE)


###########################################################################
###########################################################################
## Model1: Classification trees

set.seed(333)
fit.cart <- train(usedTrainGroup~., data=Over_Train, method="rpart", trControl=control, preProcess = c("pca","scale","center"))

test_pred_classes_cart <- predict(fit.cart, Testing, type="raw")
table(test_pred_classes_cart)
#names_prediction <- rownames(test_pred_classes_cart)
## Converting test_pred to factor instead of matrix
#test_pred_classes_cart <- as.vector(test_pred_classes_cart)
#names(test_pred_classes_cart) <- names_prediction
#test_pred_classes_cart <- ordered(test_pred_classes_cart, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_cart <- confusionMatrix(test_pred_classes_cart, usedTestGroup, positive = "Progression")
Confusion_test_cart

## ROC/AUC in the Testing set
test_pred_prob_cart <- predict(fit.cart, Testing, type = "prob")

png("./Figs/ROC_test_cart.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_cart[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test Cart")
dev.off()

#######################################################################
########################################################################
## Model2: LDA

set.seed(333)
fit.lda <- train(usedTrainGroup~., data=Over_Train, method="lda", trControl=control, preProcess = c("pca", "scale", "center"))

test_pred_classes_lda <- predict(fit.lda, Testing, type="raw")
table(test_pred_classes_lda)
#test_pred_classes_lda <- as.vector(test_pred_classes_lda)
#names(test_pred_classes_lda) <- names_prediction
#test_pred_classes_lda <- ordered(test_pred_classes_lda, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_lda <- confusionMatrix(test_pred_classes_lda, usedTestGroup, positive = "Progression")
Confusion_test_lda

## Ranking features by importance
Importance_LDA <- varImp(fit.lda, scale = FALSE)
Importance_LDA <- Importance_LDA$importance
Importance_LDA <- Importance_LDA[order(Importance_LDA$Progression, decreasing = TRUE),]
genes_LDA <- rownames(Importance_LDA)[1:500]

## ROC/AUC in the Testing set
test_pred_prob_lda <- predict(fit.lda, Testing, type = "prob")

png("./Figs/ROC_test_LDA.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_lda[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test LDA")
dev.off()

#######################################################################
########################################################################
## Model 3: SVM

set.seed(333)
fit.svm <- train(usedTrainGroup~., data=Over_Train, method="svmRadial", trControl=control)

test_pred_classes_svm <- predict(fit.svm, Testing, type="raw")
table(test_pred_classes_svm)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svm <- confusionMatrix(test_pred_classes_svm, usedTestGroup, positive = "Progression")
Confusion_test_svm

## ROC/AUC in the Testing set
test_pred_prob_svm <- predict(fit.svm, Testing, type = "prob")
png("./Figs/ROC_test_SVM.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_svm[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM")
dev.off()

Importance_SVM <- varImp(fit.svm, scale = FALSE)
Importance_SVM <- Importance_SVM$importance
Importance_SVM <- Importance_SVM[order(Importance_SVM$Progression, decreasing = TRUE),]
genes_SVM <- rownames(Importance_SVM)[1:1000]
#######################################################################
########################################################################
## Model 4: KNN

set.seed(333)
fit.knn <- train(usedTrainGroup~., data=Over_Train, method="knn", trControl=control)

test_pred_classes_knn <- predict(fit.knn, Testing, type="raw")
table(test_pred_classes_knn)
#test_pred_classes_knn <- as.vector(test_pred_classes_knn)
#names(test_pred_classes_knn) <- names_prediction
#test_pred_classes_knn <- ordered(test_pred_classes_knn, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_knn <- confusionMatrix(test_pred_classes_knn, usedTestGroup, positive = "Progression")
Confusion_test_knn

## ROC/AUC in the Testing set
test_pred_prob_knn <- predict(fit.knn, Testing, type = "prob")
png("./Figs/ROC_test_KNN.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_knn[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test KNN")
dev.off()

Importance_knn <- varImp(fit.knn, scale = FALSE)
Importance_knn <- Importance_knn$importance
Importance_knn <- Importance_knn[order(Importance_knn$Progression, decreasing = TRUE),]
genes_knn <- rownames(Importance_knn)[1:500]
#######################################################################
########################################################################
## Model 5: Naïve Bayes

set.seed(333)
fit.nb <- train(usedTrainGroup~., data=Over_Train, method="nb", trControl=control, preProcess = c("ica"))

test_pred_classes_nb <- predict(fit.nb, Testing, type="raw")
table(test_pred_classes_nb)
#test_pred_classes_rf <- as.vector(test_pred_classes_rf)
#names(test_pred_classes_rf) <- names_prediction
#test_pred_classes_rf <- ordered(test_pred_classes_rf, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_nb <- confusionMatrix(test_pred_classes_nb, usedTestGroup, positive = "Progression")
Confusion_test_nb

## ROC/AUC in the Testing set
test_pred_prob_nb <- predict(fit.nb, Testing, type = "prob")

png("./Figs/ROC_test_NB.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nb[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Naive Bayes")
dev.off()

################################################################################
## Model 6: Neural Network

set.seed(333)
fit.nnet <- train(usedTrainGroup~., data=Over_Train, method="nnet", trControl=control, preProcess=c("pca", "scale", "center"), maxit=500)

test_pred_classes_nnet <- predict(fit.nnet, Testing, type="raw")
table(test_pred_classes_nnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_nnet <- confusionMatrix(test_pred_classes_nnet, usedTestGroup, positive = "Progression")
Confusion_test_nnet

## ROC/AUC in the Testing set
test_pred_prob_nnet <- predict(fit.nnet, Testing, type = "prob")

png("./Figs/ROC_test_NNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nnet[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

################################################################################
### Comparison between the models

# collect resamples
results <- resamples(list(CART=fit.cart, LDA=fit.lda, SVM=fit.svm, KNN=fit.knn, NB=fit.nb))
summary(results)

# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/BWplot.png", width = 3000, height = 3000, res = 300)
bwplot(results, scales=scales)
dev.off()


# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/DensPlot.png", width = 3000, height = 3000, res = 300)
densityplot(results, scales=scales, pch = "|")
dev.off()

# dot plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/DotPlot.png", width = 3000, height = 3000, res = 300)
dotplot(results, scales=scales)
dev.off()

# pair-wise scatterplots of predictions to compare models
png("./Figs/SPlom.png", width = 3000, height = 3000, res = 300)
splom(results)
dev.off()

# difference in model predictions
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs)