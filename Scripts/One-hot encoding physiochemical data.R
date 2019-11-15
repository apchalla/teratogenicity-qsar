#Anup Challa
#IRTA Fellow, NIH-NCATS
#08 July 2019

#This program establishes cutoffs for quantum chemistry descriptors for use in the physiochemical layer of the teratogenicity QSAR. 

options(stringsAsFactors = FALSE)
library(pROC)

#Loading in the appropriate feature and label sets, and restricting the label set to only 609 drugs for which physiochemical information was obtained.ß
p_chem_descriptors <- read.csv('609 Sampled Drugs Physiochemical Calculations.csv')
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

hms_drugs <- read.csv('HMS QSAR Sampled Drugs.csv')
ncats_drugs <- read.csv('NCATS 609 QSAR Sampled Drugs.csv')

row_remove <- NULL
for(i in 1:nrow(hms_drugs))
{
  if(hms_drugs$db_ids[i] %in% ncats_drugs$db_ids_ncats_qsar)
  {
    
  }
  else
  {
    row_remove <- c(i, row_remove)
  }
}
ter_factors <- ter_factors[-row_remove]
ter_factors <- data.frame(ter_factors)
write.csv(ter_factors, 'NCATS 609 Sampled Drugs Teratogenicity Categories.csv')

#Enabling ROC on descriptors of interest to establish cutoffs
#For potential of an electron loss from the HOMO:
rows_with_homo <- NULL
for(i in 1:nrow(p_chem_descriptors))
{
  if(is.na(p_chem_descriptors$AM1_HOMO[i]))
  {
    rows_with_homo <- c(rows_with_homo, i)
  }
}
p_chem_descriptors <- p_chem_descriptors[-rows_with_homo,]
ter_factors <- ter_factors[-rows_with_homo,]
#Eliminating sparsity before running ROC analysis.
roc_homo <- multiclass.roc(ter_factors, p_chem_descriptors$AM1_HOMO)
roc_homo_stats <- roc_homo$rocs[[1]]
#Running ROC analysis and obtatining summary statistics for the call
roc_homo_table <- cbind(roc_homo_stats$thresholds, roc_homo_stats$sensitivities, roc_homo_stats$specificities)
#Creating a table of HOMO energy cutoffs and their associated predictive sensitivities and specificities

#Based on the resulting ROC curve, there is no "great" cutoff for HOMO energy to employ in binary categorization,
#as the ROC plot tends towards the line y = x (AUC = 0.55).
#However, the cutoff of HOMO energy at -9.570580 eV appears to optimize sensitivity (~52%) and specificity (~52%).

#For potential of an electron gain to the LUMO:
#We assume that each drug with calculated HOMO energy is accompanied with a calculation of LUMO energy
roc_lumo <- multiclass.roc(ter_factors, p_chem_descriptors$AM1_LUMO)
roc_lumo_stats <- roc_lumo$rocs[[1]]
#Running ROC analysis and obtatining summary statistics for the call
roc_lumo_table <- cbind(roc_lumo_stats$thresholds, roc_lumo_stats$sensitivities, roc_lumo_stats$specificities)
#Creating a table of LUMO energy cutoffs and their associated predictive sensitivities and specificities

#Based on the resulting ROC curve, there is no "great" cutoff for LUMO energy to employ in binary categorization,
#as the ROC plot tends towards the line y = x (AUC = 0.52).
#However, the cutoff of LUMO energy at -5.196165 eV appears to optimize sensitivity (~55%) and specificity (~53%).