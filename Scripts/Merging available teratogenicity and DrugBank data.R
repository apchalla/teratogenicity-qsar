##(1) Extracting DB IDs for Drugs with Known Teratogenicity
##(2) Identifying DB IDs for Drugs with Known Structure and Teratogenicity
##(3) Extracting SMILES Strings for (2)-Integrated Drugs
##(4) Extracting Teratogenicity Classification for (2)-Integrated Drugs

setwd("~/Desktop")
library(stringr)
options(stringsAsFactors = FALSE)
library(readxl)

library(readxl)
drug_vocabulary_with_structure_with_SMILES <- read_excel("HMS Drug Development/drug vocabulary with structure with SMILES.xlsx")
View(drug_vocabulary_with_structure_with_SMILES)

library(readxl)
drugbank_ids_with_teratogenicity <- read_excel("HMS Drug Development/drugbank ids with teratogenicity.xlsx")
View(drugbank_ids_with_teratogenicity)

db_teratogen<-c()
nrow(drugbank_ids_with_teratogenicity)
for(i in 1:nrow(drugbank_ids_with_teratogenicity))
{
  db_teratogen<-c(db_teratogen, drugbank_ids_with_teratogenicity$`DrugBank ID`[i])
}
print(db_teratogen)

common_id<-c()
nrow(drug_vocabulary_with_structure_with_SMILES)
for(j in 1:nrow(drug_vocabulary_with_structure_with_SMILES))
{
  x=NULL
  if(drug_vocabulary_with_structure_with_SMILES$`DrugBank IDs`[j] %in% db_teratogen)
  {
    x=drug_vocabulary_with_structure_with_SMILES$`DrugBank IDs`[j]
    common_id<-c(common_id,x)
  }
}
print(common_id)
df=data.frame(common_id)
write.csv(df,file='drugbank id with structure with SMILES with teratogenicity.csv')

smiles_ter<-c()
for(m in 1:nrow(drug_vocabulary_with_structure_with_SMILES))
{
  if(drug_vocabulary_with_structure_with_SMILES$`DrugBank IDs`[m] %in% common_id)
  {
    smiles_ter<-c(smiles_ter,drug_vocabulary_with_structure_with_SMILES$`final_SMILES`[m])
  }
}                                                                                                                                                                                                                                                                                                                                                                                
df_smiles=data.frame(smiles_ter)
write.csv(df_smiles, file='df_smiles.csv')

library(readxl)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("HMS Drug Development/drugbank ids with structure with SMILES with teratogenicity.xlsx")
View(drugbank_ids_with_structure_with_SMILES_with_teratogenicity)

teratogenicity_classification<-c()
for(n in 1:nrow(drugbank_ids_with_teratogenicity)){
  if(drugbank_ids_with_teratogenicity$`DrugBank ID`[n] %in% drugbank_ids_with_structure_with_SMILES_with_teratogenicity$`DrugBank IDs`)
  {
    teratogenicity_classification<-c(teratogenicity_classification, drugbank_ids_with_teratogenicity$`Teratogenicity Classification`[n])
  }
}
df=data.frame(teratogenicity_classification)
write.csv(teratogenicity_classification, file='teratogenicity.csv')