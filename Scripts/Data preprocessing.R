##Anup Challa
##Visiting Undergraduate Research Fellow, Harvard Medical School
##13 July 2018

library(stringr)
options(stringsAsFactors = FALSE)
library(readxl)

input_file = '~/Desktop/HMS Drug Development/Archives/open structures.txt'
lines = readLines(input_file)

all_drug_ids = NULL
for(i in 1:length(lines)) {
  l = lines[i]
  if(!is.na(str_match(l,'<DRUGBANK_ID>'))) {
    drug_id = lines[i+1]
    all_drug_ids <- c(all_drug_ids, drug_id)
  }
}

library(readxl)
drug_vocabulary <- read_excel("~/Desktop/HMS Drug Development/Archives/drug vocabulary.xlsx", 
                              range = "A2:A9738", col_types = c("text"))
View(drug_vocabulary)

ncol(drug_vocabulary)
nrow(drug_vocabulary)

common<-c()
x=NULL

for(j in 1:nrow(drug_vocabulary))
{
  if(drug_vocabulary$`DrugBank ID`[j] %in% all_drug_ids)
  {
    x=drug_vocabulary$`DrugBank ID`[j]
    common=c(common, x)
  }
}
setwd("~/Desktop/HMS Drug Development")
df=data.frame(common)
write.csv(df, file='Common.csv')

drug_vocabulary_with_structure_copy <- read_excel("~/Desktop/HMS Drug Development/Archives/drug vocabulary with structure copy.xlsx", 
                                                  range = "B1:B9100", col_types = c("text"))
View(drug_vocabulary_with_structure_copy)

drug_vocabulary_with_SMILES <- read_excel("~/Desktop/HMS Drug Development/Archives/drug vocabulary with SMILES.xlsx", skip = 1)
View(drug_vocabulary_with_structure_copy)

nrow(drug_vocabulary_with_structure_copy)
nrow(drug_vocabulary_with_SMILES)

dbID_known<-c()

for(j in 1:nrow(drug_vocabulary_with_structure_copy))
{
  x=drug_vocabulary_with_structure_copy$'X__1'[j]
  dbID_known<-c(dbID_known,x)
}

common_SMILES<-c()

for(k in 1:nrow(drug_vocabulary_with_SMILES))
{
  y=drug_vocabulary_with_SMILES$'DrugBank ID'[k]
  if (y %in% dbID_known)
  {
    common_SMILES<-c(common_SMILES,y)
  }
}

final_SMILES<-c()

for(m in 1:nrow(drug_vocabulary_with_SMILES))
{
  z=drug_vocabulary_with_SMILES$'DrugBank ID'[m]
  if(z %in% common_SMILES)
  {
    final_SMILES<-c(final_SMILES,drug_vocabulary_with_SMILES$'SMILES'[m])
  }
}

df_SMILES=data.frame(final_SMILES)
write.csv(df_SMILES, file='drug vocabulary with structure with SMILES.csv')

drugbank.safefetus.pregrisk.pfam.table <- read.csv("~/Desktop/HMS Drug Development/Archives/drugbank.safefetus.pregrisk.pfam.table.txt")
View(drugbank.safefetus.pregrisk.pfam.table)
write.csv(drugbank.safefetus.pregrisk.pfam.table, file='risk.csv')

drug_vocabulary_with_structure_with_SMILES <- read_excel("~/Desktop/HMS Drug Development/Archives/drug vocabulary with structure with SMILES.xlsx")
View(drug_vocabulary_with_structure_with_SMILES)

drugbank_ids_with_teratogenicity <- read_excel("~/Desktop/HMS Drug Development/Archives/drugbank ids with teratogenicity.xlsx")
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

drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("~/Desktop/HMS Drug Development/Archives/drugbank ids with structure with SMILES with teratogenicity.xlsx")
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

source("http://bioconductor.org/biocLite.R")
biocLite("ChemmineR")
library("ChemmineR")
vignette("ChemmineR")

setwd("~/Desktop/HMS Drug Development")

sdfset<-read.SDFset("~/Desktop/HMS Drug Development/Archives/open structures.sdf")
valid<-validSDF(sdfset); sdfset<-sdfset[valid]
valid[1]
validDF<-data.frame(valid)
write.csv(validDF,file='valid.csv')
remove(validDF)

View(sdfset[1])

apset<-sdf2ap(sdfset)
moreinorganic<-c(sdfset[1]@SDF[[1]]@datablock[["DRUGBANK_ID"]],0)
moreinorganic<-c(sdfset[326]@SDF[[1]]@datablock[["DRUGBANK_ID"]],0)
moreinorganic<c(sdfset[326]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[1193]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[1366]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[1599]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[3471]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[4587]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7259]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7268]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7401]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7404]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7487]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7497]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7501]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7560]@SDF[[1]]@datablock[["DRUGBANK_ID"]],sdfset[7599]@SDF[[1]]@datablock[["DRUGBANK_ID"]])

View(moreinorganic)
df=data.frame(moreinorganic)
write.csv(df,file='df')


fpset <- desc2fp(apset, descnames=4096, type="FPset") 
fpma <- as.matrix(fpset)
fpset <- as(fpma, "FPset")

remove(moreinorganic)
remove(df)

write.csv(fpma,file='fpma.csv')

setwd("~/Desktop/HMS Drug Development")

drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("~/Desktop/HMS Drug Development/Archives/drugbank ids with structure with SMILES with teratogenicity.xlsx", 
                                                                          +     range = "A1:A612")
View(drugbank_ids_with_structure_with_SMILES_with_teratogenicity)

atom_pair_fingerprints <- read_excel("atom pair fingerprints.xlsx", 
                                     +   sheet = "fpma_607", range = "A1:A608")
View(atom_pair_fingerprints)

setdiff(allID, atom_pair_fingerprints$'X__1')

drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("~/Desktop/HMS Drug Development/Archives/drugbank ids with structure with SMILES with teratogenicity.xlsx", 
                                                                          +     range = "A1:B612")
View(drugbank_ids_with_structure_with_SMILES_with_teratogenicity)

aScore<-c(rep(0,611))
for(i in 1:nrow(drugbank_ids_with_structure_with_SMILES_with_teratogenicity))
{
  if(drugbank_ids_with_structure_with_SMILES_with_teratogenicity$`Teratogenicity Classification`[i]=='A')
  {
    aScore[i]=1
  }
}
View(aScore)
df<-data.frame(aScore)
write.csv(df, file='aScore.csv')

bScore<-c(rep(0,611))
for(i in 1:nrow(drugbank_ids_with_structure_with_SMILES_with_teratogenicity))
{
  if(drugbank_ids_with_structure_with_SMILES_with_teratogenicity$`Teratogenicity Classification`[i]=='B')
  {
    bScore[i]=1
  }
}
df<-data.frame(bScore)
write.csv(df, file='bScore.csv')

cScore<-c(rep(0,611))
for(i in 1:nrow(drugbank_ids_with_structure_with_SMILES_with_teratogenicity))
{
  if(drugbank_ids_with_structure_with_SMILES_with_teratogenicity$`Teratogenicity Classification`[i]=='C')
  {
    cScore[i]=1
  }
}
df<-data.frame(cScore)
write.csv(df, file='cScore.csv')

dScore<-c(rep(0,611))
for(i in 1:nrow(drugbank_ids_with_structure_with_SMILES_with_teratogenicity))
{
  if(drugbank_ids_with_structure_with_SMILES_with_teratogenicity$'Teratogenicity Classification'[i]=='D')
  {
    dScore[i]=1
  }
}
df<-data.frame(dScore)
write.csv(df,file='dScore.csv')

xScore<-c(rep(0,611))
for(i in 1:nrow(drugbank_ids_with_structure_with_SMILES_with_teratogenicity))
{
  if(drugbank_ids_with_structure_with_SMILES_with_teratogenicity$'Teratogenicity Classification'[i]=='X')
  {
    xScore[i]=1
  }
}
df<-data.frame(xScore)
write.csv(df,file='xScore.csv')