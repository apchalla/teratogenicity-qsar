library(stringr)
options(stringAsFactors=FALSE)
library(readxl)
library(ChemmineR)
vignette("ChemmineR")
library(rcdk)
library(rcdklibs)
library(rJava)

#Fingerprint generator for 4,096-bit representaiton (max allowable in ChemmmineR)
setwd("~/Desktop/HMS Drug Development/Archive")
sdfset<-read.SDFset("open structures.sdf")
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

##Modified fingerprint generator for 512-bit representations

sdfset<-read.SDFset("open structures.sdf")
valid<-validSDF(sdfset); sdfset<-sdfset[valid]

apset<-sdf2ap(sdfset)

fpset <- desc2fp(apset, descnames=512, type="FPset") 
fpma <- as.matrix(fpset)
fpset <- as(fpma, "FPset")

write.csv(fpma,file='fpma.csv')

a<-read.csv('fpma.csv')
dim(a)
a[1:10,1:10]
a<-a[,-1]
dim(a)

write.csv(a, file='512-bit rep.csv')

full_matrix<-read.csv('~/Desktop/HMS Drug Development/Features/9,025 valid drugs atom pair fingerprints.csv')
ids<-as.vector(full_matrix$X)
ids<-data.frame(ids)
ids<-ids[1:9025,]
dim(ids)
b<-cbind(ids,a)
dim(b)
head(b)
b<-b[,-1]
dim(b)

X <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
head(X)
dim(X)

c<-data.frame()

valid_ids<-as.vector(X[,1])
length(valid_ids)

for(i in 1:nrow(b))
{
  if(b[,1][i] %in% valid_ids)
  {
    c<-rbind(c,b[i,])
  }
}
dim(c)

ap_no_ter<-data.frame()
ap_no_ter<-c
dim(ap_no_ter)

write.csv(ap_no_ter, 'ap_512_no_ter.csv')

#MORGAN fingerprint generator--512-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'circular', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

morgan_1<-cbind(db_id[,1], finger_hits[,1])
morgan_1<-data.frame(morgan_1, header=TRUE)
morgan_1<-morgan_1[,-3]
#morgan_1 is a table of DB IDs with an adjointed column of circular fingerprint hits (as a string).

morgan_1<-as.matrix(morgan_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(morgan_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='morgan finger 1,024-bit.csv')

#PUBCHEM fingerprint generator--881-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'pubchem', fp.mode= 'bit', size=881, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

pubchem_1<-cbind(db_id[,1], finger_hits[,1])
pubchem_1<-data.frame(pubchem_1, header=TRUE)
pubchem_1<-pubchem_1[,-3]
#pubchem_1 is a table of DB IDs with an adjointed column of PubChem fingerprint hits (as a string).

pubchem_1<-as.matrix(pubchem_1)

one_hot<-matrix(nrow=611, ncol=881)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:881)
  {
    pullout<-c(pullout,as.numeric(word(pubchem_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:881)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='pubchem finger 881-bit.csv')

#STANDARD fingerprint generator--1,024-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'standard', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

standard_1<-cbind(db_id[,1], finger_hits[,1])
standard_1<-data.frame(standard_1, header=TRUE)
standard_1<-standard_1[,-3]
#standard_1 is a table of DB IDs with an adjointed column of Standard fingerprint hits (as a string).

standard_1<-as.matrix(standard_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(standard_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='standard finger 1024-bit.csv')

#EXTENDED fingerprint generator--1,024-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'extended', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

extended_1<-cbind(db_id[,1], finger_hits[,1])
extended_1<-data.frame(extended_1, header=TRUE)
extended_1<-extended_1[,-3]
#extended_1 is a table of DB IDs with an adjointed column of Extended fingerprint hits (as a string).

extended_1<-as.matrix(extended_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(extended_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='extended finger 1024-bit.csv')

#GRAPH fingerprint generator--1,024-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'graph', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

graph_1<-cbind(db_id[,1], finger_hits[,1])
graph_1<-data.frame(extended_1, header=TRUE)
graph_1<-extended_1[,-3]
#graph_1 is a table of DB IDs with an adjointed column of Graph fingerprint hits (as a string).

graph_1<-as.matrix(graph_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(graph_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='graph finger 1024-bit.csv')

#HYBRIDIZATION fingerprint generator--1,024-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'hybridization', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

hybridization_1<-cbind(db_id[,1], finger_hits[,1])
hybridization_1<-data.frame(hybridization_1, header=TRUE)
hybridization_1<-hybidization_1[,-3]
#hybridization_1 is a table of DB IDs with an adjointed column of Hybridization fingerprint hits (as a string).

hybridization_1<-as.matrix(hybridization_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(hybridization_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='hybridization finger 1024-bit.csv')

#MACCS fingerprint generator--166-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'maccs', fp.mode= 'bit', size=166, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

maccs_1<-cbind(db_id[,1], finger_hits[,1])
maccs_1<-data.frame(maccs_1, header=TRUE)
maccs_1<-maccs_1[,-3]
#maccs_1 is a table of DB IDs with an adjointed column of MACCS fingerprint hits (as a string).

maccs_1<-as.matrix(maccs_1)

one_hot<-matrix(nrow=611, ncol=166)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:166)
  {
    pullout<-c(pullout,as.numeric(word(extended_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:166)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='maccs finger 166-bit.csv')

#ESTATE fingerprint generator--79-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'estate', fp.mode= 'bit', size=79, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

estate_1<-cbind(db_id[,1], finger_hits[,1])
estate_1<-data.frame(extended_1, header=TRUE)
estate_1<-estate_1[,-3]
#estate_1 is a table of DB IDs with an adjointed column of Estate fingerprint hits (as a string).

estate_1<-as.matrix(estate_1)

one_hot<-matrix(nrow=611, ncol=79)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:79)
  {
    pullout<-c(pullout,as.numeric(word(estate_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:79)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='estate finger 79-bit.csv')

#KR fingerprint generator--4,860-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'kr', fp.mode= 'bit', size=4860, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

kr_1<-cbind(db_id[,1], finger_hits[,1])
kr_1<-data.frame(kr_1, header=TRUE)
kr_1<-kr_1[,-3]
#kr_1 is a table of DB IDs with an adjointed column of KR fingerprint hits (as a string).

kr_1<-as.matrix(kr_1)

one_hot<-matrix(nrow=611, ncol=4860)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:4860)
  {
    pullout<-c(pullout,as.numeric(word(kr_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:4860)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='kr finger 4860-bit.csv')

#Signature fingerprint generator--1,024-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'signature', fp.mode= 'bit', size=1024, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

signature_1<-cbind(db_id[,1], finger_hits[,1])
signature_1<-data.frame(signature_1, header=TRUE)
signature_1<-signature_1[,-3]
#signature_1 is a table of DB IDs with an adjointed column of Signature fingerprint hits (as a string).

signature_1<-as.matrix(signature_1)

one_hot<-matrix(nrow=611, ncol=1024)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:1024)
  {
    pullout<-c(pullout,as.numeric(word(signature_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:1024)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='signature finger 1024-bit.csv')

#SHORTEST PATH fingerprint generator--512-bit represenations (using package rcdk)
drugbank_ids_with_structure_with_SMILES_with_teratogenicity <- read_excel("Desktop/HMS Drug Development/Archive/drugbank ids with structure with SMILES with teratogenicity.xlsx")
smiles_entry<-drugbank_ids_with_structure_with_SMILES_with_teratogenicity

smiles<-c()
for(i in 1:nrow(smiles_entry))
{
  smiles<-c(smiles,smiles_entry$'SMILES'[i])
}
sp=get.smiles.parser()
ac_smiles<-parse.smiles(smiles)
# parse.smiles accepts a vector of SMILES strings and returns a list of type AtomContainer, containing items of type IAtomContainer
# IAtomContainer is the only object type acceptable by the get.fingerprint() command.

finger<-c()
for(i in 1:611)
{
  finger<-c(finger, get.fingerprint(ac_smiles[[i]], type= 'shortestpath', fp.mode= 'bit', size=512, verbose=TRUE))
}
length(finger)

db_id<-data.frame()
for(i in 1:611)
{
  db_id[i,1]=smiles_entry$`DrugBank IDs`[i]
}

finger_hits<-data.frame()
for(i in 1:611)
{
  finger_hits[i,1]= toString(finger[[i]]@bits)
}

shortestpath_1<-cbind(db_id[,1], finger_hits[,1])
shortestpath_1<-data.frame(shortestpath_1, header=TRUE)
shortestpath_1<-shortestpath_1[,-3]
#shortestpath_1 is a table of DB IDs with an adjointed column of Shortest Path fingerprint hits (as a string).

shortestpath_1<-as.matrix(shortestpath_1)

one_hot<-matrix(nrow=611, ncol=512)

pullout<-c()
for(i in 1:611)
{  
  for(k in 1:512)
  {
    pullout<-c(pullout,as.numeric(word(kr_1[,2][i], k, sep = fixed(','))))
    print(i)
    print(k)
  }
  for(j in 1:512)
  {
    if(j %in% pullout)
    {
      one_hot[i,j]=1
    }
    else
    {
      one_hot[i,j]=0
    }
  }
  pullout<-c()
}
one_hot[1:15,1:15]

one_hot<-data.frame(one_hot)
one_hot<-cbind(db_id,one_hot)

setwd("~/Desktop/HMS Drug Development/Features")

write.csv(one_hot, file='shortestpath finger 512-bit.csv')