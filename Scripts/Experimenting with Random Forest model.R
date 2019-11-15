install.packages('Rcpp')
install.packages('recipes')
install.packages('broom')
install.packages('devtools')
library(devtools)
install.packages('bindrcpp')
install.packages('DEoptimR')
devtools::install_github('topepo/caret/pkg/caret')
install.packages('lava')
install.packages('caret', dependencies = T)
library(caret)
library(randomForest)

X <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
Y <- as.factor(X[,4098])
X <- X[,c(-1,-4098)]
##Noting that a training set usually contains 80% of available data and a test set uses 20% of information
Z.train.ap<-data.frame(X[1:489,])
Z.train.ter<- Y[1:489]
Z.test.ap<-data.frame(X[490:611,])
Z.test.ter<- Y[490:611]

###bestmtry<-tuneRF(Z.train.ap, Z.train.ter, stepFactor=1.5, improve=1e-5, ntree=500) 
###a = randomForest(x=Z.train.ap, y=Z.train.ter, mtry=64, ntree=500) #mtry selected from above
###summmary(a)
###print(a)
#Try an Extended Caret model once we are able to properly enable CaretR

#In collaboration with Andy, the followng code works to successfully establish the Caret package:

#Using the Extended Caret custom tuning model predicted in Machine Learning Mastery to establish unique object.
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x=Z.train.ap, y=Z.train.ter, len = NULL, search = "grid") {}
customRF$fit <- function(x=Z.train.ap, y=Z.train.ter, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x=Z.train.ap, y=Z.train.ter, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(Z.train.ap) Z.train.ap[order(Z.train.ap[,1]),]
customRF$levels <- function(Z.train.ap) Z.train.ap$classes

#Training the model after establishing the unique object (trial using TRAINING DATA; need to use full dataset for best results)
control <- trainControl(method="cv", number=3)
tunegrid <- expand.grid(.mtry=c(16,32,64,128), .ntree=c(100, 150, 200, 500))
custom <- train(x=Z.train.ap, y=Z.train.ter, method='rf', trControl=control, do.trace=10)
summary(custom)
plot(custom)
#Plot gives optimal combination of mtry=130, ntree=150 (~90% prediction accutacy)