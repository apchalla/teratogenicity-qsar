###ESTABLISHING XGBOOST, which requires matrices and can work well with sparsity. 
#We decided to use this approach after playing around with the Random Forest,
#as we note that it will likely increase CV prediction accuracy in comparsion to the RF 
#and will run much faster (i.e, it is more efficent).

library(devtools)
library(caret)
library(tidyverse)
library(rpart)
library(dplyr)
library(xgboost)
library(data.table)
library(pROC)

X <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
y <- c(as.numeric(X[,4098]))
y <- y-1
X <- X[,c(-1,-4098)]

Y <- read.csv('~/Desktop/HMS Drug Development/Labels/teratogenicity one-hot encodings.csv', skip=1, header = FALSE)
Y <- Y[,3:7]

##Taking template from RPubs multiclass classification post

numberOfClasses <- ncol(Y)
xgb_params <- list("objective" = "multi:softprob","eval_metric" = "merror", "num_class" = numberOfClasses)
nround=25
cv.nfold=10

train_matrix <- xgb.DMatrix(data=as.matrix(X), label=y)
colSums(Y)
##weights<-c(14/611,157/611,310/611,91/611,39/611) //In order to establish weights for imbalanced classes
cv_model <- xgb.cv(params = xgb_params, data = train_matrix, nrounds = nround, nfold = cv.nfold, verbose = TRUE, prediction = TRUE, max.depth=4, stratified=TRUE, eta=0.05, min_child_weight=5)
summary(cv_model)
print(cv_model)

y<-as.integer(y)
OOF_prediction <- data.frame(cv_model$pred) %>% mutate(max_prob = max.col(cv_model$pred), label = y)
confusionMatrix(factor(OOF_prediction$max_prob), factor(OOF_prediction$label), mode = "everything")

##SYSTEMATICALLY TUNING the XGBoost parameters allows us to rewrite the parameters given
#above to maximize OOF accuracy. This is achieved through the following algorithm:

ap <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_score<-ap[,4098]
ap <- ap[,2:4097]

trControl = trainControl(
  method = 'cv',
  number = 10,
  classProbs = TRUE,
  verboseIter = TRUE,
  allowParallel = TRUE)

tuneGridXGB <- expand.grid(
  nrounds=c(350),
  max_depth = c(3, 4, 6),
  eta = c(0.05, 0.1),
  gamma = c(0,0.01),
  colsample_bytree = c(0.75),
  subsample = c(0.50),
  min_child_weight = c(0, 5))

xgbmod <- train(
  x = ap,
  y = as.factor(ter_score),
  method = 'xgbTree',
  metric = 'merror',
  trControl = trControl,
  tuneGrid = tuneGridXGB
)

preds <- predict(xgbmod, newdata = ap, type = "prob")

print(xgbmod$results)
print(xgbmod$resample)

max_prob=max.col(preds)

resROC<-multiclass.roc(as.numeric(ter_score), max_prob)
rs<-resROC[['rocs']]
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
auc(resROC)

#XGBoost and ROC stats on 512-bit fingerprints

library(devtools)
library(caret)
library(tidyverse)
library(rpart)
library(dplyr)
library(xgboost)
library(data.table)
library(pROC)

X <- read.csv('~/Desktop/HMS Drug Development/Archive/ap_512_no_ter.csv', skip=1, header = FALSE)
X <- X[,-c(1,2)]

y <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header=FALSE)
y <- c(as.numeric(y[,4098]))
y <- y-1

Y <- read.csv('~/Desktop/HMS Drug Development/Labels/teratogenicity one-hot encodings.csv', skip=1, header = FALSE)
Y <- Y[,3:7]

##Taking template from RPubs multiclass classification post

numberOfClasses <- ncol(Y)
xgb_params <- list("objective" = "multi:softprob","eval_metric" = "merror", "num_class" = numberOfClasses)
nround=350
cv.nfold=10

train_matrix <- xgb.DMatrix(data=as.matrix(X), label=y)
colSums(Y)
##weights<-c(14/611,157/611,310/611,91/611,39/611) //In order to establish weights for imbalanced classes
cv_model <- xgb.cv(params = xgb_params, data = train_matrix, nrounds = nround, nfold = cv.nfold, verbose = TRUE, prediction = TRUE, max.depth=4, stratified=TRUE, eta=0.05, min_child_weight=5)
summary(cv_model)
print(cv_model)

y<-as.integer(y)
OOF_prediction <- data.frame(cv_model$pred) %>% mutate(max_prob = max.col(cv_model$pred), label = y)
confusionMatrix(factor(OOF_prediction$max_prob), factor(OOF_prediction$label), mode = "everything")

resROC<-multiclass.roc(OOF_prediction$label, OOF_prediction$max_prob)
rs<-resROC[['rocs']]
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=i))
for(i in 1:5)
{
  print (auc(rs[[i]]))
}
auc(resROC)

##XGBoost and ROC stats on 9,025 drugs: binary classification for presence of teratogenicity score

#SYSTEMATICALLY TUNING the XGBoost parameters allows us to rewrite the parameters given
#above to maximize OOF accuracy. This is achieved through the following algorithm:

ap<-read.csv('~/Desktop/HMS Drug Development/Features/atom pair fingerprints binary.csv', skip=1, header = FALSE)
ap<-ap[-c(9026:9099),]
ter_score<-ap[,4098]
ap<-ap[,2:4097]

trControl = trainControl(
  method = 'cv',
  number = 10,
  classProbs = TRUE,
  verboseIter = TRUE,
  allowParallel = TRUE)

tuneGridXGB <- expand.grid(
  nrounds=c(350),
  max_depth = c(3, 4, 6),
  eta = c(0.05, 0.1),
  gamma = c(0,0.01),
  colsample_bytree = c(0.75),
  subsample = c(0.50),
  min_child_weight = c(0, 5))

xgbmod <- train(
  x = ap,
  y = make.names(ter_score),
  method = 'xgbTree',
  metric = 'merror',
  trControl = trControl,
  tuneGrid = tuneGridXGB,
  verbose=TRUE
)

preds <- predict(xgbmod, newdata = ap, type = "prob")

print(xgbmod$results)
print(xgbmod$resample)

max_prob=max.col(preds)

resROC<-multiclass.roc(as.numeric(ter_score), max_prob)
rs<-resROC[['rocs']]
plot.roc(rs[[1]])
auc(resROC)

#Are our results meaningful?
table(ter_score)
table(ter_score)/sum(table(ter_score))