#Anup Challa
#IRTA Fellow, NIH-NCATS
#13 June 2019

#This program is a supplement to PhysiochemicalClustering.R in further scoping the physiochemical data layer for the teratogenicity QSAR.

options(stringsAsFactors = FALSE)
library(Rtsne)
library(ggplot2)
library(gbm)
library(caret)

#From previous experimentation, we know that structure provides preliminary signal for the
#prediction of teratogenicity. Though druglike or leadlike properties may not provide a strong basis for
#predicting teratogenicity, the combination of structure, druglike, and leadlike properties may provide a
#promising feature set. We will therefore combine all chemical fingerprints and physiochemical data into a large
#feature set and perform feature selection to identify the most important features to use in a predictive
#model. We then look at the predictive power of this set. We continue to examine a subset of 416 drugs, so 
#we may expand our analysis to all 609 drugs if we observe sufficient signal.

#Note that GBM was attempted on the RO3/RO5 data set during hyperparam optimization as a signal check, and the predictive 
#accuracy was less than that with structure.

#Before we combine data sets, however, let's try the following: 
# 1. Obtain druglike, leadlike, quantum chemistry (HOMO and LUMO energies), mutagenicity (Ames test), and most-basic
#pKa data to be used as an enhanced feature set. The most-basic pKa information can be obtained from the Resolver.
# 2. One-hot encode the information and re-try t-SNE and GBM. Observe how many features have non-zero importance.
# 3. Attempt feature selection on the physiochemical data set, and re-try step 2.

phys_chem_descriptors <- read.csv('609 Sampled Drugs MOE Physiochemical Calculations.csv')
head(phys_chem_descriptors)
rows_with_homo <- NULL
for(i in 1:nrow(phys_chem_descriptors))
{
  if(is.na(phys_chem_descriptors$AM1_HOMO[i]))
  {
    rows_with_homo <- c(rows_with_homo, i)
  }
}
phys_chem_descriptors <- phys_chem_descriptors[-rows_with_homo,]
rows_with_na_pka <- NULL
for(i in 1:nrow(phys_chem_descriptors))
{
  if(is.na(phys_chem_descriptors$most_basic_pKa[i]))
  {
    rows_with_na_pka <- c(rows_with_na_pka, i)
  }
}
phys_chem_descriptors <- phys_chem_descriptors[-rows_with_na_pka,]
#Loading in the feature set and eliminating drugs with sparsity in energy calculations (either HOMO and LUMO energies are both present or
#both absent) and most basic pKa values of NaN

#Now, we one-hot encode data that is not already encoded in this format. First, we attempt this procedure for HOMO and LUMO energies using the
#ROC cutoffs we previously determined:
# HOMO energy < -9.570580 eV -> "1" (i.e., unstable compound, losing lot of energy in loss of an electron)
# LUMO energy < -5.196165 eV -> "1" (i.e., unstable compound, losing lot of energy in transition of electron to LUMO)

for(i in 1:nrow(phys_chem_descriptors))
{
  if(phys_chem_descriptors$AM1_HOMO[i] < -9.570580)
  {
    phys_chem_descriptors$AM1_HOMO[i] = 1
  }
  else
  {
    phys_chem_descriptors$AM1_HOMO[i] = 0
  }
}
for(i in 1:nrow(phys_chem_descriptors))
{
  if(phys_chem_descriptors$AM1_LUMO[i] < -5.196165)
  {
    phys_chem_descriptors$AM1_LUMO[i] = 1
  }
  else
  {
    phys_chem_descriptors$AM1_LUMO[i] = 0
  }
}
#Next, we use code we previously employed for assessing druglikeness to establish one-hot encoded data. We define "1" as a state
#adherent to Lipinski's heuristic. Since we did not observe significant differences in model performance between druglike and leadlike heuristic
#application in a previous signal check, we will use druglikeness in this evaluation, is it is a more liberal standard than leadlikeness.Therefore, 
#TPSA is also unused, as further application of Veber's rule criteria to the Ro5 will reduce the number of predictive hits.
for(i in 1:nrow(phys_chem_descriptors))
{
  if(phys_chem_descriptors$a_don[i] <= 5)
  {
    phys_chem_descriptors$a_don[i] = 1
  }
  if(phys_chem_descriptors$a_don[i] > 5)
  {
    phys_chem_descriptors$a_don[i] = 0
  }
  if(phys_chem_descriptors$a_acc[i] <= 10)
  {
    phys_chem_descriptors$a_acc[i] = 1
  }
  if(phys_chem_descriptors$a_acc[i] > 10)
  {
    phys_chem_descriptors$a_acc[i] = 0
  }
  if(phys_chem_descriptors$SlogP[i] <= 5)
  {
    phys_chem_descriptors$SlogP[i] = 1
  }
  if(phys_chem_descriptors$SlogP[i] > 5)
  {
    phys_chem_descriptors$SlogP[i] = 0
  }
  if(phys_chem_descriptors$Weight[i] < 500)
  {
    phys_chem_descriptors$Weight[i] = 1
  }
  if(phys_chem_descriptors$Weight[i] >= 500)
  {
    phys_chem_descriptors$Weight[i] = 0
  }
}
phys_chem_descriptors <- phys_chem_descriptors[,-12] #Dropping consideration of rotatable bonds, as leadlikeness is not considered.
#We finally apply most-basic pKa encoding by pKa < 7 and pKa > 7 (tendency for deprotonation in acidic or basic solution). We arbitrarily
#set "1" for compounds that deprotonate in basic solution.
for(i in 1:nrow(phys_chem_descriptors))
{
  if(phys_chem_descriptors$most_basic_pKa[i] < 7)
  {
    phys_chem_descriptors$most_basic_pKa[i] = 1
  }
  if(phys_chem_descriptors$most_basic_pKa[i] > 7)
  {
    phys_chem_descriptors$most_basic_pKa[i] = 0
  }
}

#We now import teratogenicity scores and establish the appropriate three-pronged label scheme.
ter_factors <- read.csv('NCATS 609 Sampled Drugs Teratogenicity Categories.csv')
ter_factors <- ter_factors[-rows_with_homo,]
ter_factors <- ter_factors[-rows_with_na_pka,]

#We now attempt t-SNE and GBM on this set of features and labels:
features_tsne <- cbind(phys_chem_descriptors,ter_factors)
tsne_fit <- Rtsne(features_tsne[,-c(1:7,17)], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=features_tsne[,17],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])
ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()
#This t-sNE shows NO SIGNAL. Therefore, ESTABLISHMENT OF A GBM FOR THIS FEATURE SET IS UNWARRANTED. Is MOE data of 
#low fidelity?

#Therefore: new plan. Let's take HOMO energy, LUMO energy, mutagenicity, and most basic pKa from the new physiochemical feature set
#(these are the descriptors that seem mechanistically promising among the new set (pKa) and are difficult to generate outside MOE (energies,
#Ames test data)), add them to the existing druglike set from Resolver, and perform feature selection.
#Then, run t-SNE and GBM again on the new set of physiochemical information.
features_redefined <- cbind(phys_chem_descriptors$DRUGBANK_ID, phys_chem_descriptors$AM1_HOMO, phys_chem_descriptors$AM1_LUMO, phys_chem_descriptors$mutagenic, phys_chem_descriptors$most_basic_pKa)

#Importing Resolver druglikeness features:
resolver_features <- read.csv('609 Sampled Drugs Resolver Physiochemical Calculations.csv')
#Reading in and visualizing physiochemical data obtained from the Resolver app

for(i in 1:nrow(resolver_features))
{
  if(resolver_features[i,6] <= 5)
  {
    resolver_features[i,6] = 1
  }
  if(resolver_features[i,6] > 5)
  {
    resolver_features[i,6] = 0
  }
  if(resolver_features[i,5] <= 10)
  {
    resolver_features[i,5] = 1
  }
  if(resolver_features[i,5] > 10)
  {
    resolver_features[i,5] = 0
  }
  if(resolver_features[i,4] <= 5)
  {
    resolver_features[i,4] = 1
  }
  if(resolver_features[i,4] > 5)
  {
    resolver_features[i,4] = 0
  }
  if(resolver_features[i,3] < 500)
  {
    resolver_features[i,3] = 1
  }
  if(resolver_features[i,3] >= 500)
  {
    resolver_features[i,3] = 0
  }
}
resolver_features <- resolver_features[,-c(7,8)]

index_to_remove <- NULL
for(i in 1:nrow(resolver_features))
{
  if(resolver_features$DRUGBANK_ID[i] %in% features_redefined[,1])
  {
    
  }
  else
  {
    index_to_remove <- c(index_to_remove, i)
  }
}
resolver_features <- resolver_features[-index_to_remove,]
feature_set <- cbind(resolver_features, features_redefined[,2:5])
colnames(feature_set)[7] <- "AM1_HOMO"
colnames(feature_set)[8] <- "AM1_LUMO"
colnames(feature_set)[9] <- "mutagenic"
colnames(feature_set)[10] <- "most_basic_pKa"
#We have now established a complete list of all physiochemical descriptors which we believe have the strongest predictive capacity
#for teratogenicity class.
label_set <- ter_factors[,2]
#Our label set is established as above

#We now perform simple feature selection on this set of predictors. There are two ways we may do this in Caret: either by removal
#of correlated features (i.e., r^2 >= 0.75), or by a ranking of feature importance to a supervised mdoel and manual selection of the
#number and identity of features to retain. Given how few features we employ in this layer, we are unlikely to take the second approach.

#Therefore, we obtain a correlation matrix of features. We must use data in its raw form for this analysis (i.e., data that is not one-hot
#encoded).
raw_MOE <- read.csv('609 Sampled Drugs MOE Physiochemical Calculations.csv')
raw_resolver <- read.csv('609 Sampled Drugs Resolver Physiochemical Calculations.csv')

rows_to_remove <- NULL
for(i in 1:nrow(raw_MOE))
{
  if(raw_MOE$DRUGBANK_ID[i] %in% feature_set$DRUGBANK_ID)
  {
    
  }
  else
  {
    rows_to_remove <- c(rows_to_remove, i)
  }
}
raw_MOE <- raw_MOE[-rows_to_remove,]
raw_resolver <- raw_resolver[-rows_to_remove,]
#Loading in raw data and restricting to drugs with all available features
raw_data <- cbind(raw_resolver, raw_MOE$AM1_HOMO, raw_MOE$AM1_LUMO, raw_MOE$mutagenic)
raw_data <- raw_data[,-8]
#Compiling and cleaning all data elements of the raw feature set
cor(raw_data[,-c(1:3)])
#From the ouput correlation matrix, it does not appear that there are highly correlated features (r^2 >= 0.75) that need to be mitigated in the feature set.

#Therefore, we will perform t-SNE with the integrated, uncorrelated data set
tsne_set <- cbind(feature_set,ter_factors)
tsne_set <- tsne_set[,-c(1,2,11)]
tsne_fit <- Rtsne(tsne_set[,-9], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=tsne_set[,9],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])
ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()
#There is extremely minimal class separation evident in t-SNE with this data set. Therefore, we will not attempt GBM to obtain predictive accuracy.
#We will then combine these physiochemical parameters with structural information, perform feature selection, and then re-run t-SNE and GBM on the integrated feature set.

#We now attempt to combine structural and physiochemical descriptor sets together to observe if a combined set of this information is more predictive
#of teratogenicity than purely structure: 
structure_set <- read.csv('morgan finger 1,024-bit.csv')
rows_to_remove <- NULL
for(i in 1:nrow(structure_set))
{
  if(structure_set$V1[i] %in% feature_set$DRUGBANK_ID)
  {
    
  }
  else
  {
    rows_to_remove <- c(rows_to_remove,i)
  }
}
structure_set <- structure_set[-rows_to_remove,-c(1,2)]
combined_struc_pc <- cbind(structure_set, feature_set)
#NOTE: Because fingerprints are inherently categorical (factored), we cannot enable a correlation matrix for feature selection.
#Therefore, we must develop a learning vector quantization (LVQ) model and then remove features of little importance in the model:
#We do not tune the hyperparameters rigorously, as we study only selection of features of significant importance--not predictive accuracy--
#with the below LVQ model:
#Code borrowed from J. Brownlee, https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/.
#We first one-hot encode the labels
# clin_acceptable <- NULL
# clin_c <- NULL
# clin_unacceptable <- NULL
# for(i in 1:length(ter_factors$x))
# {
#   if(ter_factors$x[i] == 'A/B')
#   {
#     clin_acceptable <- c(clin_acceptable, "1")
#     clin_c <- c(clin_c, "0")
#     clin_unacceptable <- c(clin_unacceptable, "0")
#   }
#   if(ter_factors$x[i] == 'C')
#   {
#     clin_acceptable <- c(clin_acceptable, "0")
#     clin_c <- c(clin_c, "1")
#     clin_unacceptable <- c(clin_unacceptable, "0")    
#   }
#   if(ter_factors$x[i] == 'D/X')
#   {
#     clin_acceptable <- c(clin_acceptable, "0")
#     clin_c <- c(clin_c, "0")
#     clin_unacceptable <- c(clin_unacceptable, "1")
#   }
# }
#ter_scores <- cbind(clin_acceptable, clin_c, clin_unacceptable)

ter_factors <- as.factor(ter_factors[,2])
ter_scores <- NULL
for(i in 1:length(ter_factors))
{
  if(ter_factors[i] == 'A/B' | ter_factors[i] == 'C')
  {
    ter_scores <- c(ter_scores, "0")
  }
  if(ter_factors[i] == 'D/X')
  {
    ter_scores <- c(ter_scores, "1")
  }
}
ter_scores <- as.factor(ter_scores)

#We now ensurize the feature set is factorized.
for(i in 1:ncol(combined_struc_pc))
{
  combined_struc_pc[,i] = as.factor(combined_struc_pc[,i])
}
#For the LVQ:
control <- trainControl(method="repeatedcv", number=10, repeats=3)
model <- train(x = combined_struc_pc, y = ter_scores, method="lvq", preProcess="scale", trControl=control, tuneGrid = data.frame(size = 3, k = 1:2))