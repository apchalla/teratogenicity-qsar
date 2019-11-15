#Anup Challa
#IRTA Fellow, NIH-NCATS
#21 June 2019

#This program explores combining A & B teratogenicity classifications into a "clinically acceptable" variable,
#in order to test if this system of label organization provides for greater predictive signal across the 
#training set by curbing class imbalance.

options(stringsAsFactors = FALSE)
library(pracma)
library(Rtsne)
library(ggplot2)
library(ChemmineR)
library(caret)

#For Morgan fingerprints (purely structure):
#Reading, visualizing, and pre-processing the feature set (for the original 611 DrugBank drugs)
morgan_finger <- read.csv('morgan finger 1,024-bit.csv')
View(morgan_finger[1:10, 1:10])
morgan_finger <- morgan_finger[,-c(1,2)]
#Reading in and visualizing the set of labels
ter_scores <- read.csv('teratogenicity one-hot encodings.csv')
ter_scores <- ter_scores[,-c(1,2)]

#We now want to consider only four (4) teratogenicity scores: "Clinically Acceptable" (i.e., A and B), C, D, X, in an
#attempt to combat natural class imbalance within our data set.
clin_acceptable <- NULL

for(i in 1:nrow(ter_scores))
{
  if(ter_scores[i,1] == 1 | ter_scores[i,2] == 1)
  {
    clin_acceptable <- c(clin_acceptable, 1)
  }
  else
  {
    clin_acceptable <- c(clin_acceptable, 0)
  }
}
ter_scores <- cbind(clin_acceptable,ter_scores[,c(3,4,5)])
head(ter_scores)

#With features and labels thus appropriately labelled for a classification task, we now attempt a multiclass t-SNE to 
#check for predictive signal. We first start by creating label factors for the t-SNE call. We also observe the class
#balance within the factors vector.
ter_factors <- NULL
for(i in 1:nrow(ter_scores))
{
  if(ter_scores[i,1] == 1)
  {
    ter_factors <- c(ter_factors, 'A/B')
  }
  if(ter_scores[i,2] == 1)
  {
    ter_factors <- c(ter_factors, 'C')
  }
  if(ter_scores[i,3] == 1)
  {
    ter_factors <- c(ter_factors, 'D')
  }
  if(ter_scores[i,4] == 1)
  {
    ter_factors <- c(ter_factors, 'X')
  }
}
ter_factors <- as.factor(ter_factors)
count_acceptable = 0
count_C = 0
count_D = 0
count_X = 0
for(i in 1:length(ter_factors))
{
  if(ter_factors[i] == "A/B")
  {
    count_acceptable = count_acceptable + 1
  }
  if(ter_factors[i] == "C")
  {
    count_C = count_C + 1
  }
  if(ter_factors[i] == "D")
  {
    count_D = count_D + 1
  }
  if(ter_factors[i] == "X")
  {
    count_X = count_X + 1
  }
}
#This lumping obviously addresses class imbalance for non-teratogenic compounds, but class D and X drugs are still
#imbalanced in the scope of the entire data set.
#Running Barnes-Hut-SNE:
tsne_set <- cbind(morgan_finger, ter_factors)
tsne_fit <- Rtsne(tsne_set[,-1025], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=as.factor(tsne_set[,1025]),D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])
ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()
#The resulting t-SNE clusters are not extremely tight; however, more distinct clusters are visible than previously observed.

#Training the predictive algorithm on the new system of classification.
ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = c(3,5),
                         n.trees = c(100,200,300,400), 
                         shrinkage = c(0.1,0.01),
                         n.minobsinnode=2)
tree_model <- train(x=morgan_finger, y=ter_factors,
                    method = "gbm",
                    trControl=ctrl,
                    tuneGrid = tune_grid,
                    metric='accuracy',
                    maximize=TRUE)
plot(tree_model)

#What about further lumping the D/X label sets together, as a "Clinically Unacceptable" label?
ter_scores <- read.csv('teratogenicity one-hot encodings.csv')
ter_scores <- ter_scores[,-c(1,2)]
clin_acceptable <- NULL
clin_unacceptable <- NULL
for(i in 1:nrow(ter_scores))
{
  if(ter_scores[i,1] == 1 | ter_scores[i,2] == 1)
  {
    clin_acceptable <- c(clin_acceptable, 1)
    clin_unacceptable <- c(clin_unacceptable, 0)
  }
  if(ter_scores[i,4] == 1 | ter_scores[i,5] == 1)
  {
    clin_acceptable <- c(clin_acceptable, 0)
    clin_unacceptable <- c(clin_unacceptable, 1)
  }
  if(ter_scores[i,3] == 1)
  {
    clin_acceptable <- c(clin_acceptable, 0)
    clin_unacceptable <- c(clin_unacceptable, 0)
  }
}
ter_scores <- cbind(clin_acceptable,ter_scores[,3],clin_unacceptable)
head(ter_scores)
ter_factors <- NULL
for(i in 1:nrow(ter_scores))
{
  if(ter_scores[i,1] == 1)
  {
    ter_factors <- c(ter_factors, 'A/B')
  }
  if(ter_scores[i,2] == 1)
  {
    ter_factors <- c(ter_factors, 'C')
  }
  if(ter_scores[i,3] == 1)
  {
    ter_factors <- c(ter_factors, 'D/X')
  }
}
ter_factors <- as.factor(ter_factors)
count_acceptable = 0
count_C = 0
count_unacceptable = 0
for(i in 1:length(ter_factors))
{
  if(ter_factors[i] == "A/B")
  {
    count_acceptable = count_acceptable + 1
  }
  if(ter_factors[i] == "C")
  {
    count_C = count_C + 1
  }
  if(ter_factors[i] == "D/X")
  {
    count_unacceptable = count_unacceptable + 1
  }
}
remove(tsne_set)
tsne_set <- cbind(morgan_finger, ter_factors)
tsne_fit <- Rtsne(tsne_set[,-1025], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=tsne_set[,1025],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])
ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()
ctrl <- trainControl(method = "cv", number=5)
# tune_grid <- expand.grid(interaction.depth = c(3,5),
#                          n.trees = c(100,200,300,400),
#                          shrinkage = c(0.1,0.01),
#                          n.minobsinnode=2)
# tree_model <- train(x=morgan_finger, y=ter_factors,
#                     method = "gbm",
#                     trControl=ctrl,
#                     tuneGrid = tune_grid,
#                     metric='accuracy',
#                     maximize=TRUE)
# plot(tree_model)

#Observe the sensitivity, ROC, etc. and re-obtain this information for A,B,C,D,X
#classification. Compare and determine which method is better and produces results of 
#higher impact.

#When Gergely looked at the plot of results from our original GBM (structure-based, A,B,C,D,X),
#he strongly suggested increasing the depth of each layer exponentially and observing the results.
#We do this by re-creating a tuning grid with a range of interaction depthts.

#Manually, increasing the interaction depth to 25 and using a lumped A/B label increases the predictive accuracy
#to 64%!

#We require a systematic method of hyperparameter optimization!

#Let's try a hyperparam optimization exercise with random search, as opposed to using a tuning grid.
#The following code, as borrowed from 
#https://topepo.github.io/caret/random-hyperparameter-search.htmlhttps://topepo.github.io/caret/random-hyperparameter-search.html,
#provides a mechanism for a random search within CaretR. We continue to employ five (5)-fold CV and set
#tuneLength (the maximum number of combinations of hyperparameters that will be tested) to the factorial of the
#number of available hyperparams (i.e., 4!)

#We note that random search fails at the edge conditions (tested as follows: shrinkage <- c(0.01, 0.6), interaction depth <- c(1, 24),
#iterations <- c(1, 4900), min_obs_node = 2), which is likely a bug in Caret. The result is either failure at a tier of testing or an endless loop.

#Thererfore, we create a tuning grid with an arbitrary selection of hyperparam values across the range of values which can be accomodated by the
#train() function. We note that the tuning grid becomes overwhelmed at boundary values (n.treees -> 3000-4000, interaction.depth -> 50), and when too 
#many values are provided to the grid. Therefore, we restrict the test values we pass to grid (in both size and quantity)

# ctrl <- trainControl(method = "cv", number=5)
# tune_grid <- expand.grid(interaction.depth = c(1,3,5,10,15,20),
#                          n.trees = c(100,500,1000,3000),
#                          shrinkage = c(0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.3, 0.5),
#                          n.minobsinnode=2)
# tree_model <- train(x=morgan_finger, y=ter_factors,
#                     method = "gbm",
#                     trControl=ctrl,
#                     tuneGrid = tune_grid,
#                     metric='accuracy',
#                     maximize=TRUE)
# plot(tree_model)

ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = 25,
                         n.trees = 2500,
                         shrinkage = 0.01,
                         n.minobsinnode=2)
tree_model <- train(x=morgan_finger, y=ter_factors,
                    method = "gbm",
                    trControl=ctrl,
                    tuneGrid = tune_grid,
                    metric='accuracy',
                    maximize=TRUE)

# fitControl <- trainControl(method = "cv",
#                            number = 5,
#                            search = "random")
# tree_model <- train(x=morgan_finger[,-c(1,2)], y=as.factor(ter_factors),
#                     method = "gbm",
#                     trControl=fitControl,
#                     metric='accuracy',
#                     tuneLength = 24,
#                     maximize=TRUE)
# plot(tree_model)

#With the uncommented version above and redefined label bins, predictive accuracy on structure alone -> 64.3%.

#We then implement feature selection on the structural fingerprints to reduce potential overfiting from highly correlated
#features (i.e., features with zero importance).

#ROC analysis on three-pronged GBM classification model:
preds <- predict(tree_model, newdata = morgan_finger, type = "prob")
max_prob=max.col(preds)
resROC<-multiclass.roc(as.numeric(ter_score), max_prob)
rs<-resROC[['rocs']]
auc(resROC)
#Obtaining an average AUC for all ROC iterations = 0.78
rs[[1]]
auc(rs[[1]])
rs[[2]]
auc(rs[[2]])
rs[[3]]
auc(rs[[3]])
#Noting that rs[[1]] contains results homologous with those from the average of resROC, we plot rs[[1]] as our ROC curve
#for this GBM analysis.
roc_plot <- ggroc(rs[[1]]) + geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), colour = 'gray')
roc_plot + ggtitle('Average Multiclass ROC Outcome for Morgan Fingerprint-Teratogenicity Score GBM')