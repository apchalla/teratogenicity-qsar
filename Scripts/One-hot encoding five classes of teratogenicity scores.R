setwd("~/Desktop/HMS Drug Development")
library(readxl)

drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/drugbank ids with structure with SMILES with teratogenicity.xlsx", 
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