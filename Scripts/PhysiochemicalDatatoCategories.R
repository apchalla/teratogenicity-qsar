#Anup Challa
#IRTA Fellow, NIH-NCATS
#13 June 2019

#This program converts physiochemical data on compounds in our teratogenicity
#QSAR--as obtained from NCATS Resolver--to one-hot encoded categoies for use in
#a GBM.

options(stringsAsFactors = FALSE)
library(pracma)

features_no_mm <- read.csv('NCGC Candidates Physiochemical Information.csv')
mol_mass <- read.csv('NCATS Molecular Mass.csv')
head(features)
head(mol_mass)
features <- cbind(features_no_mm, mol_mass)
#Reading in and visualizing physiochemical data obtained from the Resolver app

for(i in 1:nrow(features))
{
  if(features[i,2] <= 5)
  {
    features[i,2] = 1
  }
  if(features[i,2] > 5)
  {
    features[i,2] = 0
  }
  if(features[i,3] <= 10)
  {
    features[i,3] = 1
  }
  if(features[i,3] > 10)
  {
    features[i,3] = 0
  }
  if(features[i,5] <= 5)
  {
    features[i,5] = 1
  }
  if(features[i,5] > 5)
  {
    features[i,5] = 0
  }
  if(features[i,8] < 500)
  {
    features[i,8] = 1
  }
  if(features[i,8] >= 5)
  {
    features[i,8] = 0
  }
}
features <- features[,-c(1,4,6,7)]
drug_id_sheet <- read.csv('NCATS-DB Sample Molecule Pairs and Similarities.csv')
#features <- cbind(drug_id_sheet$PubChem.ID, features)
features <- cbind(drug_id_sheet$DrugBank.ID, features)
#Preparing a cleaned, one-hot encoded features matrix. Encoding is based on
#Lipinski's RO5 (druglike = '1', not druglike = '0')

#First, we wish to evaluate the predictive power of RO5 encodings for teratogenicity.
#Therefore, we consider the data frame of features above and restrict ourselves to
#exact match cases (T = 1) between NCGC compounds and drugs originally evaluated in the
#teratogenicity QSAR.

#From manual inspection of drug_id_sheet, the first 415 cases within the data
#correspond to exact matches. Therefore, we restrict our analysis to those compounds
#for which we have available teratogenicity information and NCATS screening data,
#in hopes of integrating assay data at a later stage.

features_ncats <- features[1:416,]

#We must now load in our label set and restrict it to our compounds of interest:
labels_611_drugs <- read.csv('teratogenicity one-hot encodings.csv')

ter_scores <- NULL

for(i in 1:nrow(features_ncats))
{
  for(j in 1:nrow(labels_611_drugs))
  {
    for(k in 3:7)
    {
      if(strcmp(features_ncats$`drug_id_sheet$DrugBank.ID`[i], labels_611_drugs$X__1[j]))
      {
        if(labels_611_drugs[j,k] == 1)
        {
          ter_scores <- c(ter_scores, k-2)
        }
      }
    }
  }
  print(i)
}

labels_raw <- cbind(features_ncats$`drug_id_sheet$DrugBank.ID`, ter_scores)

#Now, we one-hot encode our labels.
labels = matrix(nrow = nrow(labels_raw), ncol = 5)
for(i in 1:nrow(labels))
{
  for(j in 1:ncol(labels))
  {
    labels[i,j] = 0
  }
}