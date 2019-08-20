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
library(randomForest)

#######################################################################

## Load data
load("./Objs/ProgressionDataGood2.rda")

## Load the selected genes
#TF_MiR <- load("/Users/mohamedomar/Desktop/TF_MiR.rda")

### Common genes
#keepGns <- intersect(as.vector(myTSPs), rownames(mixTrainMat))
Progression_genes <- read.csv("/Users/mohamedomar/Desktop/BAD_VS_good.csv")
Progression_genes <- apply(Progression_genes, 2, as.character)

### Common genes
keepGns <- intersect(as.vector(Progression_genes), rownames(mixTrainMat))


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[keepGns, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train
Correl <- cor(Training)
dim(Correl)
HighlyCorrel <- findCorrelation(Correl, cutoff = 0.75)
TrainingFilt <- Training[,-HighlyCorrel]
dim(TrainingFilt)


## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(TrainingFilt)
all(rownames(TrainingFilt) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
TrainingFilt <- as.data.frame(TrainingFilt)
Data_train <- cbind(TrainingFilt, usedTrainGroup)

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
Data_train[,332] <- as.factor(Data_train[,332])
Over_Train <- SMOTE(usedTrainGroup~., data = Data_train, perc.over = 300, perc.under = 134)
table(Over_Train[,332])

######
control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, search = "grid", summaryFunction = twoClassSummary)


###########################################################################
###########################################################################
## Model1: C5.0

## Fit the model
set.seed(333)
fit.cart <- train(usedTrainGroup~., data=Over_Train, method="C5.0", trControl=control, metric = "ROC")
fit.cart

## Predict classes
test_pred_classes_cart <- predict(fit.cart, Testing, type="raw")
table(test_pred_classes_cart)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_cart <- confusionMatrix(test_pred_classes_cart, usedTestGroup, positive = "Progression")
Confusion_test_cart

## ROC/AUC in the Testing set
test_pred_prob_cart <- predict(fit.cart, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_cart.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_cart[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test Cart")
dev.off()

## Top 100 predictors
Importance_cart <- varImp(fit.cart, scale = FALSE)
Importance_cart <- Importance_cart$importance
Importance_cart <- as.matrix(Importance_cart)
Importance_cart <- Importance_cart[order(Importance_cart[,1], decreasing = TRUE),]
genes_cart <- names(Importance_cart)[1:100]
save(genes_cart, file = "./Objs/prognosis_genes/genes_cart.rdata")

#######################################################################
########################################################################
## Model 2: GLMnet

set.seed(333)
fit.GLMnet <- train(usedTrainGroup~., data=Over_Train, method="glmnet", trControl=control, metric = "ROC")
fit.GLMnet

## Predict classes
test_pred_classes_GLMnet <- predict(fit.GLMnet, Testing, type="raw")
table(test_pred_classes_GLMnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_GLMnet <- confusionMatrix(test_pred_classes_GLMnet, usedTestGroup, positive = "Progression")
Confusion_test_GLMnet

## ROC/AUC in the Testing set
test_pred_prob_GLMnet <- predict(fit.GLMnet, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_GLMnet.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_GLMnet[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test GLMnet")
dev.off()

## Top 100 predictors
Importance_GLMnet <- varImp(fit.GLMnet, scale = FALSE)
Importance_GLMnet <- Importance_GLMnet$importance
Importance_GLMnet <- as.matrix(Importance_GLMnet)
Importance_GLMnet <- Importance_GLMnet[order(Importance_GLMnet[,1], decreasing = TRUE),]
genes_GLMnet <- names(Importance_GLMnet)[1:100]
save(genes_GLMnet, file = "./Objs/Prognosis_genes/genes_GLMnet.rdata")

################################################################################
################################################################################
## Model 3: LDA

set.seed(333)
fit.lda <- train(usedTrainGroup~., data=Over_Train, method="lda2", trControl=control, metric = "ROC")
fit.lda

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
genes_LDA <- rownames(Importance_LDA)[1:100]
save(genes_LDA, file = "./Objs/Prognosis_genes/genes_LDA.rdata")
## ROC/AUC in the Testing set
test_pred_prob_lda <- predict(fit.lda, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_LDA.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_lda[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test LDA")
dev.off()

#######################################################################
########################################################################
## Model 4: High Dimensional Discriminant Analysis (HDDA)		

set.seed(333)
fit.hdda <- train(usedTrainGroup~., data=Over_Train, method="hdda", trControl=control, metric = "ROC")
fit.hdda

test_pred_classes_hdda <- predict(fit.hdda, Testing, type="raw")
table(test_pred_classes_hdda)
#test_pred_classes_lda <- as.vector(test_pred_classes_lda)
#names(test_pred_classes_lda) <- names_prediction
#test_pred_classes_lda <- ordered(test_pred_classes_lda, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_hdda <- confusionMatrix(test_pred_classes_hdda, usedTestGroup, positive = "Progression")
Confusion_test_hdda

## Ranking features by importance
Importance_HDDA <- varImp(fit.hdda, scale = FALSE)
Importance_HDDA <- Importance_HDDA$importance
Importance_HDDA <- Importance_HDDA[order(Importance_HDDA$Progression, decreasing = TRUE),]
genes_HDDA <- rownames(Importance_HDDA)[1:100]
save(genes_HDDA, file = "./Objs/Prognosis_genes/genes_HDDA.rdata")
## ROC/AUC in the Testing set
test_pred_prob_hdda <- predict(fit.hdda, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_HDDA.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_hdda[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ROC Test HDDA")
dev.off()

#######################################################################
########################################################################
## Model 5: SVM Poly

set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Over_Train, method="svmPoly", trControl=control, metric = "ROC")
fit.svmPoly

test_pred_classes_svmPoly <- predict(fit.svmPoly, Testing, type="raw")
table(test_pred_classes_svmPoly)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly <- confusionMatrix(test_pred_classes_svmPoly, usedTestGroup, positive = "Progression")
Confusion_test_svmPoly

## ROC/AUC in the Testing set
test_pred_prob_svmPoly <- predict(fit.svmPoly, Testing, type = "prob")
png("./Figs/Prognosis_genes/ROC_test_SVMPoly.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_svmPoly[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly")
dev.off()

## Top 100 predictors
Importance_SVMPoly <- varImp(fit.svmPoly, scale = FALSE)
Importance_SVMPoly <- Importance_SVMPoly$importance
Importance_SVMPoly <- Importance_SVMPoly[order(Importance_SVMPoly$Progression, decreasing = TRUE),]
genes_SVMPoly <- rownames(Importance_SVMPoly)[1:100]
save(genes_SVMPoly, file = "./Objs/Prognosis_genes/genes_SVMPloy.rdata")

#######################################################################
########################################################################
## Model 6: SVM Radial

set.seed(333)
fit.svmRadial <- train(usedTrainGroup~., data=Over_Train, method="svmRadial", trControl=control, metric = "ROC")
fit.svmRadial

test_pred_classes_svmRadial <- predict(fit.svmRadial, Testing, type="raw")
table(test_pred_classes_svmRadial)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmRadial <- confusionMatrix(test_pred_classes_svmRadial, usedTestGroup, positive = "Progression")
Confusion_test_svmRadial

## ROC/AUC in the Testing set
test_pred_prob_svmRadial <- predict(fit.svmRadial, Testing, type = "prob")
png("./Figs/Prognosis_genes/ROC_test_SVMRadial.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_svmRadial[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Radial")
dev.off()

## Top 100 predictors
Importance_SVMRadial <- varImp(fit.svmRadial, scale = FALSE)
Importance_SVMRadial <- Importance_SVMRadial$importance
Importance_SVMRadial <- Importance_SVMRadial[order(Importance_SVMRadial$Progression, decreasing = TRUE),]
genes_SVMRadial <- rownames(Importance_SVMRadial)[1:100]
save(genes_SVMRadial, file = "./Objs/Prognosis_genes/genes_SVMRadial.rdata")

#######################################################################
########################################################################
## Model 7: Linear Distance Weighted Discrimination		

set.seed(333)
fit.dwd <- train(usedTrainGroup~., data=Over_Train, method="dwdLinear", trControl=control, metric = "ROC")
fit.dwd

test_pred_classes_dwd	<- predict(fit.dwd, Testing, type="raw")
table(test_pred_classes_dwd)
#test_pred_classes_svm <- as.vector(test_pred_classes_svm)
#names(test_pred_classes_svm) <- names_prediction
#test_pred_classes_svm <- ordered(test_pred_classes_svm, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_dwd <- confusionMatrix(test_pred_classes_dwd, usedTestGroup, positive = "Progression")
Confusion_test_dwd

## ROC/AUC in the Testing set
test_pred_prob_dwd <- predict(fit.dwd, Testing, type = "prob")
png("./Figs/Prognosis_genes/ROC_test_dwd.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_dwd[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main = "ROC Test DWD")
dev.off()

## Top 100 predictors
Importance_DWD <- varImp(fit.dwd, scale = FALSE)
Importance_DWD <- Importance_DWD$importance
Importance_DWD <- Importance_DWD[order(Importance_DWD$Progression, decreasing = TRUE),]
genes_DWD <- rownames(Importance_DWD)[1:100]
save(genes_DWD, file = "./Objs/Prognosis_genes/genes_DWD.rdata")

#######################################################################
########################################################################
## Model 8: Naïve Bayes

set.seed(333)
fit.nb <- train(usedTrainGroup~., data=Over_Train, method="nb", trControl=control, metric = "ROC")
fit.nb

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

png("./Figs/Prognosis_genes/ROC_test_NB.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nb[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Naive Bayes")
dev.off()

## Top 100 predictors
Importance_NB <- varImp(fit.nb, scale = FALSE)
Importance_NB <- Importance_NB$importance
Importance_NB <- Importance_NB[order(Importance_NB$Progression, decreasing = TRUE), ]
genes_NB <- rownames(Importance_NB)[1:100]
save(genes_NB, file = "./Objs/Prognosis_genes/genes_NB.rda")

################################################################################
## Model 9: Neural Network

set.seed(333)
fit.nnet <- train(usedTrainGroup~., data=Over_Train, method="nnet", trControl=control, maxit=500, metric = "ROC")
fit.nnet

test_pred_classes_nnet <- predict(fit.nnet, Testing, type="raw")
table(test_pred_classes_nnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_nnet <- confusionMatrix(test_pred_classes_nnet, usedTestGroup, positive = "Progression")
Confusion_test_nnet

## ROC/AUC in the Testing set
test_pred_prob_nnet <- predict(fit.nnet, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_NNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_nnet[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_nnet <- varImp(fit.nnet, scale = FALSE)
Importance_nnet <- Importance_nnet$importance
Importance_nnet <- as.matrix(Importance_nnet)
Importance_nnet <- Importance_nnet[order(Importance_nnet[,1], decreasing = TRUE),]
genes_nnet <- names(Importance_nnet)[1:100]
save(genes_nnet, file = "./Objs/Prognosis_genes/genes_nnet.rda")

###########################################################################
#########################################################################
## Model 10: Model Averaged Neural Network	

set.seed(333)
fit.AVnnet <- train(usedTrainGroup~., data=Over_Train, method="avNNet", trControl=control, maxit=500, metric = "ROC")
fit.AVnnet

test_pred_classes_AVnnet <- predict(fit.AVnnet, Testing, type="raw")
table(test_pred_classes_AVnnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_AVnnet <- confusionMatrix(test_pred_classes_AVnnet, usedTestGroup, positive = "Progression")
Confusion_test_AVnnet

## ROC/AUC in the Testing set
test_pred_prob_AVnnet <- predict(fit.AVnnet, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_AVNNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_AVnnet[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_AVnnet <- varImp(fit.AVnnet, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_AVnnet, top = 20)
Importance_AVnnet <- Importance_AVnnet$importance
Importance_AVnnet <- Importance_AVnnet[order(Importance_AVnnet$Progression, decreasing = TRUE),]
genes_AVnnet <- rownames(Importance_AVnnet)[1:100]
save(genes_AVnnet, file = "./Objs/Prognosis_genes/genes_AVnnet.rda")

#############################################################################
#############################################################################
## Model 11: Neural Networks with Feature Extraction	

set.seed(333)
fit.PCAnnet <- train(usedTrainGroup~., data=Over_Train, method="pcaNNet", trControl=control, maxit=500, metric = "ROC")
fit.PCAnnet

test_pred_classes_PCAnnet <- predict(fit.PCAnnet, Testing, type="raw")
table(test_pred_classes_PCAnnet)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_PCAnnet <- confusionMatrix(test_pred_classes_PCAnnet, usedTestGroup, positive = "Progression")
Confusion_test_PCAnnet

## ROC/AUC in the Testing set
test_pred_prob_PCAnnet <- predict(fit.PCAnnet, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_PCA_NNET.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_PCAnnet[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Neural Network")
dev.off()

## Top 100 predictors
Importance_PCAnnet <- varImp(fit.PCAnnet, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_PCAnnet, top = 20)
Importance_PCAnnet <- Importance_PCAnnet$importance
Importance_PCAnnet <- Importance_PCAnnet[order(Importance_PCAnnet$Progression, decreasing = TRUE),]
genes_PCAnnet <- rownames(Importance_PCAnnet)[1:100]
save(genes_PCAnnet, file = "./Objs/Prognosis_genes/genes_PCAnnet.rda")

#############################################################################
#############################################################################
## Model 12: Random Forest

set.seed(333)
fit.rf <- train(usedTrainGroup~., data=Over_Train, method="rf", trControl=control, prox=TRUE, metric = "ROC")
fit.rf

test_pred_classes_rf <- predict(fit.rf, Testing, type="raw")
table(test_pred_classes_rf)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_rf <- confusionMatrix(test_pred_classes_rf, usedTestGroup, positive = "Progression")
Confusion_test_rf

## ROC/AUC in the Testing set
test_pred_prob_rf <- predict(fit.rf, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_RF.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_rf[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test Random Forest")
dev.off()

Importance_RF <- varImp(fit.rf, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_RF, top = 20)
Importance_RF <- Importance_RF$importance
Importance_RF <- as.matrix(Importance_RF)
Importance_RF <- Importance_RF[order(Importance_RF[,1], decreasing = TRUE),]
genes_RF <- names(Importance_RF)[1:100]
save(genes_RF, file = "./Objs/Prognosis_genes/genes_RF.rdata")

###########################################################################
## Model 14: Stochastic Gradient Boosting

set.seed(333)
fit.gbm <- train(usedTrainGroup~., data=Over_Train, method="gbm", trControl=control, metric = "ROC")
fit.gbm

test_pred_classes_gbm <- predict(fit.gbm, Testing, type="raw")
table(test_pred_classes_gbm)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_gbm <- confusionMatrix(test_pred_classes_gbm, usedTestGroup, positive = "Progression")
Confusion_test_gbm

## ROC/AUC in the Testing set
test_pred_prob_gbm <- predict(fit.gbm, Testing, type = "prob")

png("./Figs/Prognosis_genes/ROC_test_gbm.png", width = 2000, height = 2000, res = 200)
roc(usedTestGroup, test_pred_prob_gbm[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "ROC Test GBM")
dev.off()

Importance_gbm <- varImp(fit.gbm, scale = FALSE)
## Plot the top 20 predictors
plot(Importance_gbm, top = 20)
Importance_gbm <- Importance_gbm$importance
Importance_gbm <- as.matrix(Importance_gbm)
Importance_gbm <- Importance_gbm[order(Importance_gbm[,1], decreasing = TRUE),]
genes_gbm <- names(Importance_gbm)[1:500]
save(genes_gbm, file = "./Objs/Prognosis_genes/genes_gbm.rdata")

################################################################################
### Comparison between the models

# collect resamples
results <- resamples(list(CART=fit.cart, GLM_net = fit.GLMnet, LDA=fit.lda, HDDA = fit.hdda, SVM_Poly=fit.svmPoly, SVM_Radial = fit.svmRadial, LDWD = fit.dwd, NB=fit.nb, NNET=fit.nnet, NNET_AV = fit.AVnnet, NNET_PCA = fit.PCAnnet, RF=fit.rf, GBM=fit.gbm))
summary(results)

# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/BWplot.png", width = 3000, height = 3000, res = 300)
bwplot(results, scales=scales)
dev.off()

# density plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/DensPlot.png", width = 3000, height = 3000, res = 300)
densityplot(results, scales=scales, pch = "|")
dev.off()

# dot plots of accuracy
scales <- list(x=list(relation="free"), y=list(relation="free"))

png("./Figs/Prognosis_genes/DotPlot.png", width = 3000, height = 3000, res = 300)
dotplot(results, scales=scales)
dev.off()

# pair-wise scatterplots of predictions to compare models
png("./Figs/Prognosis_genes/SPlom.png", width = 3000, height = 3000, res = 300)
splom(results)
dev.off()

# difference in model predictions
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs)
