#Anup Challa
#IRTA Fellow, NIH-NCATS
#03 July 2019

#This program extracts SDF files for all 609 drugs with known structure in DrugBank and available teratogenicity information from SafeFetus.
#The library of SDF files will then be imported to MOE in order to obtain a set of physiochemical descriptors which may be used for the second
#layer of our QSAR analysis.

#This code is a modified version of information from the program "NCATSLibraryQSARDrugsSimScore.R."

options(stringsAsFactors = FALSE)
library(ChemmineR)
library(cinf)

ncats_pharma_collection <- read.SDFset("npc.sdf")
oncodrug_collection <- read.SDFset("mipe.sdf")
annotated_lib <- read.SDFset("npact.sdf")
View(ncats_pharma_collection)
View(oncodrug_collection)
View(annotated_lib)
#Reading in and visualizing the provided data sets using ChemmineR

db_structs <- read.SDFset("open structures.sdf")
valid <- validSDF(db_structs); db_structs <- db_structs[valid]
#All drug structures employed in the first stage of fingerprinting in the 
#teratogenicity QSAR--loaded to obtain 611 structure files in ChemmineR-compatible format

sampled_drugs <- read.csv('HMS QSAR Sampled Drugs.csv')
#This data set contains the entirety of sampled drugs for the QSAR, with DB ID for
#easy querying
sampled_drugs[1:10,]
db_ids_611 <- sampled_drugs[,3]

structs_req_sdf_index <- NULL
#Vector to be filled with 611 indices within DB structure set with matching DB ids to sampled
#drugs

for(i in 1:length(db_structs))
{
  if(db_structs[i]@SDF[[1]]@datablock[["DRUGBANK_ID"]] %in% db_ids_611)
  {
    structs_req_sdf_index <- c(structs_req_sdf_index, i)
  }
}

db_ids_ncats_qsar <- NULL
for(i in 1:length(structs_req_sdf_index))
{
  db_ids_ncats_qsar <- c(db_ids_ncats_qsar, db_structs[structs_req_sdf_index[i]]@SDF[[1]]@datablock[["DRUGBANK_ID"]])
}
db_vocab <- read.csv('drug vocabulary.csv')
db_vocab[1:10,1:10]
db_names_ncats_qsar <- NULL
for(i in 1:nrow(db_vocab))
{
  if(db_vocab[i,1] %in% db_ids_ncats_qsar)
  {
    db_names_ncats_qsar <- c(db_names_ncats_qsar, db_vocab[i,2])
  }
}
new_sampled_drugs <- cbind(db_ids_ncats_qsar, db_names_ncats_qsar)
write.csv(new_sampled_drugs, 'NCATS 609 QSAR Sampled Drugs.csv')

structs_req_sdf <- NULL

for(i in 1:length(db_structs))
{
  if(i %in% structs_req_sdf_index)
  {
    structs_req_sdf <- c(structs_req_sdf, db_structs[[i]])
  }
}

for(i in 1:length(structs_req_sdf))
{
  write.SDF(structs_req_sdf[[i]], file = paste(i, ".sdf"))
}
#Writing out library of sampled drug structures (the list of SDFs cannot be written out
#as one complete file, so it is extracted by component into a folder). These structures will
#then be combined into a unified library of drug structures in StarDrop, before being fed to MOE
#to generate a new feature set.