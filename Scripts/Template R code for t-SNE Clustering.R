library(Rtsne)
library(ggplot2)

##Attempting t-SNE clustering on 611 molecules with known teratogenicity scores

#Data processing to obtain atom pair fingerprints with binary metric of availability for teratogenicity score
ap <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap<-ap[,-1]
tsne_fit <- Rtsne(ap[,-4097], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=ap[,4097],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

##Clustering all 9,025 DB-curated molecules by presence of teratogenicity score

full<-read.csv('~/Desktop/HMS Drug Development/Features/full.csv', header=FALSE)
yes<-read.csv('~/Desktop/HMS Drug Development/Features/yes list.csv', header=FALSE)
dim(full)
dim(yes)

knowter<-c()
for(i in 1:nrow(full))
{
  if(full$V1[i] %in% yes$V1)
  {
    knowter<-c(knowter, "YES")
  }
  else
  {
    knowter<-c(knowter, "NO")
  }
}

length(knowter)

setwd("~/Desktop")
df=data.frame(knowter)
write.csv(df, file='atom pair fingerprints binary_yes_no_list.csv')

#t-SNE on 9,025 molecules on basis of availability of teratogenicity score

ap_bin<- read.csv('~/Desktop/HMS Drug Development/Features/atom pair fingerprints binary.csv', skip=1, header = FALSE)
ap_bin<-ap_bin[-c(9026:9099),-1]
tsne_fit <- Rtsne(ap_bin[,-4097], verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=ap_bin[,4097],D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

##Attempting t-SNE clustering on 611 molecules with known teratogenicity scores and 512-bit AP

ap <- read.csv('~/Desktop/HMS Drug Development/Archive/ap_512_no_ter.csv', skip=1, header = FALSE)
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap<-ap[,-1:2]
ter_full<-ter_full[,4098]
tsne_fit <- Rtsne(ap, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

##Data visualization for 611 molecules of known teratogenicity (with MORGAN fingerprints)
#Tuning: perplexity=40
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-ter_full[,4098]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

#Tuning: perplexity=10
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-ter_full[,4098]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=10, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

#Tuning: perplexity=80
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-ter_full[,4098]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=80, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

p<-ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()
p + scale_color_grey() + theme_classic()

##Attempting t-SNE for binary classification (C,D,X=YES; A,B=NO) ("kitchen sink" approach)
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.numeric(ter_full[,4098])
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

#1. Atom-pair fingerprints (4,096-bit)
ap <- read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ap<-ap[,-1]
ap<-ap[,-4097]
tsne_fit <- Rtsne(ap, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#2. Atom-pair fingerprints (512-bit)
ap<-read.csv('~/Desktop/HMS Drug Development/Archive/ap_512_no_ter.csv', skip=1, header = FALSE)
ap<-ap[,-c(1,2)]
tsne_fit <- Rtsne(ap, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df <- data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#3. Morgan fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#4. PubChem fingerprints (881-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/pubchem finger 881-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#5. Standard fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/standard finger 1024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#6. Extended fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/extended finger 1024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#7. Graph fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/graph finger 1024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#8. Hybridization fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/hybridization finger 1024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#9. MACCS fingerprints (166-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/maccs finger 166-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#10. Estate fingerprints (79-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/estate finger 79-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#11. KR fingerprints (4,860-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/kr finger 4860-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#12. Shortest Path fingerprints (512-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/shortestpath finger 512-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#PERPLEXITY TOO LARGE TO RUN NO MATTER INPUTTED VALUE; TRAINING ALGORITHM DOES NOT WORK ON FINGERPRINT DATASET
#13. Signature fingerprints (1,024-bit)
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/signature finger 1024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

#Analyzing binary t-SNE clusters (YES/NO) from Morgan fingerprints
ap_bin<-read.csv('~/Desktop/HMS Drug Development/Features/morgan finger 1,024-bit.csv', header = TRUE)
ap_bin<-ap_bin[,-c(1:2)]
ter_full<-read.csv('~/Desktop/HMS Drug Development/Features/611 viable drugs atom pair fingerprints.csv', skip=1, header = FALSE)
ter_full<-as.numeric(ter_full[,4098])
for(i in 1:length(ter_full)){
  if(ter_full[i]==3){
    ter_full[i]='YES'
  }
  if(ter_full[i]==4){
    ter_full[i]='YES'
  }
  if(ter_full[i]==5){
    ter_full[i]='YES'
  }
  if(ter_full[i]==1){
    ter_full[i]='NO'
  }
  if(ter_full[i]==2){
    ter_full[i]='NO'
  }
}

tsne_fit<-Rtsne(ap_bin, verbose = TRUE, check_duplicates = FALSE, perplexity=40, max_iter= 2000)
tsne_df<-data.frame(Class=ter_full,D1=tsne_fit$Y[,1],D2=tsne_fit$Y[,2])

ggplot(tsne_df,aes(x=D1, y=D2, color=Class)) + geom_point()

coordinates<-tsne_fit$'Y'
coordinates<-data.frame(coordinates)
setwd("~/Desktop/HMS Drug Development/Archive")
write.csv(coordinates, file='t-SNE Coordinates.csv')