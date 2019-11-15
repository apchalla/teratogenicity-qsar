#Anup Challa
#IRTA Fellow, NIH-NCATS
#18 July 2019

#This program probes the third layer of MADRE's teratogenicity QSAR, in exploring Tox21 data relevant to targets upon which teratogenicity converges.
#We explore the predictive power of biological assay information for this target list, as avaiable through the public-facing browser of Tox21.

options(stringsAsFactors = FALSE)

library(pracma)
library(plyr)
library(Rtsne)
library(ggplot2)
library(caret)

tox_21_rar <- read.csv('RAR antagonist Tox21.csv')
tox_21_hdac <- read.csv('HDAC Tox21.csv')
ncgc_matches <- read.csv('NCATS-DB Sample Molecule Pairs and Similarities.csv')
ncgc_matches <- ncgc_matches[1:416,]
#Considering only exact match NCGC compounds to our list of drugs

tox_21_rar <- cbind(tox_21_rar$ï..SAMPLE_ID, tox_21_rar$ASSAY_OUTCOME, tox_21_rar$CURVE_CLASS2, tox_21_rar$AC50, tox_21_rar$EFFICACY)
tox_21_hdac <- cbind(tox_21_hdac$ï..SAMPLE_ID, tox_21_hdac$ASSAY_OUTCOME, tox_21_hdac$CURVE_CLASS2, tox_21_hdac$AC50, tox_21_hdac$EFFICACY)
#Loading in assay data and NCGC candidates of interest, and restricting assay data to only fields of interest

index_to_keep_rar <- NULL
index_to_keep_hdac <- NULL
for(i in 1:nrow(tox_21_hdac))
{
  if(tox_21_hdac[,1][i] %in% ncgc_matches$PubChem.ID)
  {
    index_to_keep_hdac <- c(index_to_keep_hdac, i)
  }
  else
  {
    
  }
}
for(i in 1:nrow(tox_21_rar))
{
  if(tox_21_rar[,1][i] %in% ncgc_matches$PubChem.ID)
  {
    index_to_keep_rar <- c(index_to_keep_rar, i)
  }
  else
  {
    
  }
}
ncgc_tox_21_hdac <- tox_21_hdac[index_to_keep_hdac,1]
ncgc_tox_21_rar <- tox_21_rar[index_to_keep_rar,1]
ncgc_to_keep <- intersect(ncgc_tox_21_rar, ncgc_tox_21_hdac)
#Identifying eligible entries for both RAR and HDAC data sets, per sample ID match to NCGC compounds.
#We have a very small list of matching IDs: only 81 drugs have avaialable structure, NCGC match, and assay data on RAR and HDAC antagonism.

#We now curate data for each matching drug: 
tox_21_rar_eligible <- NULL
for(i in 1:nrow(tox_21_rar))
{
  if(tox_21_rar[i,1] %in% ncgc_to_keep)
  {
    tox_21_rar_eligible <- rbind(tox_21_rar_eligible, tox_21_rar[i,])
  }
  else
  {
    
  }
}
tox_21_hdac_eligible <- NULL
for(i in 1:nrow(tox_21_hdac))
{
  if(tox_21_hdac[i,1] %in% ncgc_to_keep)
  {
    tox_21_hdac_eligible <- rbind(tox_21_hdac_eligible, tox_21_hdac[i,])
  }
  else
  {
    
  }
}

#We now one-hot encode the data by the NON-RIGOROUS definition (owing to the sparsity of the data). We score "0" as CC 4, and "1" as any ranking except 4.
cc_hdac_eligible <- NULL
for(i in 1:nrow(tox_21_hdac_eligible))
{
  if(tox_21_hdac_eligible[i,3] == 4)
  {
    cc_hdac_eligible <- c(cc_hdac_eligible,"0")
  }
  else
  {
    cc_hdac_eligible <- c(cc_hdac_eligible,"1")
  }
}
cc_rar_eligible <- NULL
for(i in 1:nrow(tox_21_rar_eligible))
{
  if(tox_21_rar_eligible[i,3] == 4)
  {
    cc_rar_eligible <- c(cc_rar_eligible,"0")
  }
  else
  {
    cc_rar_eligible <- c(cc_rar_eligible,"1")
  }
}
cc_hdac_eligible <- as.factor(cc_hdac_eligible)
cc_rar_eligible <- as.factor(cc_rar_eligible)
activity_hdac <- cbind(ncgc_tox_21_hdac, cc_hdac_eligible)
activity_rar <- cbind(ncgc_tox_21_rar, cc_rar_eligible)

#We now collate the parsed HDAC and RAR data sets above to generate a new feature set:
features_bio <- NULL
features_assay <- NULL
for(i in 1:nrow(activity_hdac))
{
  for(j in 1:nrow(activity_rar))
  {
    if(strcmp(activity_hdac[i,1], activity_rar[j,1]))
    {
      features_assay <- cbind(activity_rar[j,1], activity_hdac[i,2], activity_rar[j,2])
      features_bio <- rbind(features_bio, features_assay)
    }
  }
}
features_bio <- unique(features_bio)
features_bio[,2] <- revalue(features_bio[,2], c("1" = "0", "2" = "1"))
features_bio[,3] <- revalue(features_bio[,3], c("1" = "0", "2" = "1"))
#We have now compiled a one-hot encoded feature set for HDAC and RAR antagonism, using a non-rigorous definition of drug activity.
#Is there inherent redundancy in the compound data set?
#To make sure that there is no redundancy within the drug set, we wil manually scope for redundant entries, using the principal that one "positive"
#hit is sufficient to call a drug class "1."
write.csv(features_bio, 'redundant_assay_features.csv')
#Reading in the same data set, with duplicates parsed, and columns named appropriately for the bound data.
cleaned_features_bio <- read.csv('redundant_assay_features.csv')

#We are now ready for assay-specific t-SNE. First, we load the labels associated with each entry in the feature set:
#We first match NCGC IDs to DB IDs, in order to access teratogenicity scores, which are linked to DB IDs.
pbc_ids_matched <- NULL
for(i in 1:nrow(cleaned_features_bio))
{
  for(j in 1:nrow(ncgc_matches))
  {
    if(strcmp(cleaned_features_bio[i,1], ncgc_matches$PubChem.ID[j]))
    {
      pbc_ids_matched <- c(pbc_ids_matched, ncgc_matches$PubChem.ID[j])
    }
  }
}
pbc_ids_matched <- unique(pbc_ids_matched)
drugbank_ids <- NULL
for(i in 1:length(pbc_ids_matched))
{
  length_check <- ncgc_matches[ncgc_matches$PubChem.ID == pbc_ids_matched[i], "DrugBank.ID"]
  if(length(length_check > 1))
  {
    print("LONG")
  }
  drugbank_ids <- c(drugbank_ids, ncgc_matches[ncgc_matches$PubChem.ID == pbc_ids_matched[i], "DrugBank.ID"])
}
drugbank_ids <- unique(drugbank_ids)
#There is a one-off error here. Let's write out the list of IDs, find the problematic ID from looking at the 
#returned DB IDs for each call to fetch the IDs per NCGC compound, map the NCGC ID to its drug name in the qHTS client,
#and determine the correct, corresponding DB ID for the drug, parsing the extraneous ID. We do this manually and reload the
#parsed set.
write.csv(data.frame(unique(drugbank_ids)), "DB IDs with RAR-HDAC.csv")
db_covered <- read.csv("DB IDs with RAR-HDAC.csv")
nrow(db_covered)
cleaned_features_bio <- cbind(db_covered, cleaned_features_bio)

#Now, we obtain our label set:
drugs_sampled <- read.csv('NCATS 609 QSAR Sampled Drugs.csv')
nums_pull_score <- NULL
for(i in 1:nrow(cleaned_features_bio))
{
  nums_pull_score <- c(nums_pull_score, drugs_sampled[drugs_sampled$db_ids_ncats_qsar == cleaned_features_bio$unique.drugbank_ids.[i], "X"])
}
three_tier_ter_score <- read.csv('NCATS 609 Sampled Drugs Teratogenicity Categories.csv')
ter_scores <- three_tier_ter_score[nums_pull_score,]
ter_scores <- ter_scores[,2]
assay_set <- cbind(cleaned_features_bio, ter_scores)

#Creating and plotting the t-SNE object. Note the EXTREMELY SMALL SAMPLE SIZE and the NEGATIVE ERROR FROM T-SNE FITTING.
tsne_set <- assay_set[,-5]
tsne_fit <- Rtsne(tsne_set[,-c(1,2)], verbose = TRUE, check_duplicates = FALSE, perplexity=21, max_iter= 2000)
tsne_df <- data.frame(Class=assay_set[,5],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])
ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#Attempting GBM on an integrated structure-assay set:
#Integrating the feature set
structure <- read.csv('morgan finger 1,024-bit.csv')
structure <- structure[,-1]
kept_structures <- NULL
for(i in 1:nrow(db_covered))
{
  for(j in 1:nrow(structure))
  {
    if(strcmp(db_covered$unique.drugbank_ids.[i], structure$V1[j]))
    {
      kept_structures <- rbind(kept_structures, structure[j,])
    }
  }
}
kept_structures <- kept_structures[,-1]
features <- cbind(kept_structures, assay_set[,3:4])

#Training the GBM (hyperparams used based on those from optimizing solely structure)
ctrl <- trainControl(method = "cv", number=5)
tune_grid <- expand.grid(interaction.depth = 25,
                         n.trees = 2500,
                         shrinkage = 0.01,
                         n.minobsinnode=2)
tree_model <- train(y = factor(ter_scores),x = features, 
                    method = "gbm",
                    trControl=ctrl,
                    tuneGrid = tune_grid,
                    metric='accuracy',
                    maximize=TRUE)
#From this model, we observe a large drop in predictive accuracy: 51.7%. This is likely owing to the substantial drop in sample size. Additionally,
#49 feature variables had no variation, suggesting they should be parsed from the model. We re-train the model by removing these features (feature selection)
#and again enabling GBM on structure + assay results. We also note that ONLY 437 feature variables had non-zero influence, so we consider that our method of feature
#selection is quite preliminary! Note that feature selection on categorical variables is not feasible here.

tree_model <- train(y = factor(ter_scores),x = features[-c(6,9,10,12,17,30,31,36,38,40,42,44,45,45,48,54,59,66,69,74,78,79,80,81,84,87,93,95,96,112,115,116,121,132,133,135,137,142,147,149,152,154,155,157,159,160,164,167,169)], 
                    method = "gbm",
                    trControl=ctrl,
                    tuneGrid = tune_grid,
                    metric='accuracy',
                    maximize=TRUE)
#Again, we get 49 more variables to remove--from the 587 zero-influence variables in the above model. Accuracy drops to 49.9%!

#Finally, we integrate physiochemical information to establish a cumulative model.
#Therefore, we obtain a correlation matrix of features. We must use data in its raw form for this analysis (i.e., data that is not one-hot
#encoded).
#THE FOLLOWING CODE IS COPIED FROM "NEWPHYSIOCHEMICALCLUSTERING + GBM.R":
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
  if(resolver_features[i,1] %in% features_redefined[,1])
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

#We now restrict the frame of physiochemical information to solely drugs with assay data coverage.
db_covered <- as.character(db_covered[,1])
db_with_pchem <- as.character(feature_set[,1])
integrated_features_index <- NULL
for(i in 1:length(db_covered))
{
  for(j in 1:length(db_with_pchem))
  {
    if(strcmp(db_covered[i], db_with_pchem[j]))
    {
      integrated_features_index <- c(integrated_features_index, j)
    }
    else
    {
      
    }
  }
}
#There are only 21 entries with structure, teratogenicity information, all available physiochemical information, teratogenicity data, and assay coverage!
#This is certainly not sufficient n to run any meaningul ML algorithms.
#We therefore do not attempt t-SNE and GBM on this set.