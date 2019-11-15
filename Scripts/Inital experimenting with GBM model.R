library(ChemmineR)
library(dplyr)
library(Rtsne)
library(ggplot2)
library(gbm)
library(caret)
library(rcdk)
library(rcdklibs)
library(rJava)
library(readxl)
options(stringsAsFactors = FALSE)

#3. GBM on Morgan fingerprints
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}
ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    trControl=ctrl,
                    tuneGrid = tune_grid,
                    metric='accuracy',
                    maximize=TRUE)
plot(tree_model)

#4. GBM on PubChem Fingerprints
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/pubchem finger 881-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl,
                    metric='accuracy',
                    maximize=TRUE)
plot(tree_model)

#12. GBM on Shortest Path fingerprints
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/shortestpath finger 512-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl,
                    metric='accuracy',
                    maximize=TRUE)
plot(tree_model)

##Curious to see how model works on the smallest available feature matrix
#1O. GBM on Estate fingerprints
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/estate finger 79-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl)
plot(tree_model)

#1O. GBM on AP (512-bit) fingerprints
ap<-read.csv('~/Desktop/HMS Drug Development/Features/ap_512_no_ter.csv', skip=1, header = FALSE)
ap<-ap[,-c(1,2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl)
plot(tree_model)

#1O. GBM on AP (4,096-bit) fingerprints
ap <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap<-ap[,-1]
ap<-ap[,-4097]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.factor(ter_full[,4098])
#Binary classification
ter_full<-as.numeric(ter_full)
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3, 5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl)
plot(tree_model)

##Noting that Morgan fingerprints present the highest predictive accuracy, generate ROC statistics:
#Use hyperparameter values that maximize accuracy, as given in the tuning plot

ctrl <- trainControl(method = "cv", number=10)
tune_grid <- expand.grid(interaction.depth = 5,
                         n.trees = 200, 
                         shrinkage = 0.01,
                         n.minobsinnode=2)

tree_model <- train(x=ap_bin, y=ter_full,
                    method = "gbm",
                    tuneGrid = tune_grid,
                    trControl = ctrl)

preds <- predict(tree_model, newdata = ap_bin, type = "prob")

print(tree_model$results)
print(tree_model$resample)

max_prob=max.col(preds)

confusionMatrix(factor(max_prob), factor(ter_full), mode = "everything")

#Obtaining ROC statistics for binary classification model:
for(i in 1:length(ter_full)){
  if(ter_full[i]=='YES'){
    ter_full[i]=1
  }
  if(ter_full[i]=='NO'){
    ter_full[i]=2
  }
}
result.roc <- roc(ter_full, max_prob) # Draw ROC curve.
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)