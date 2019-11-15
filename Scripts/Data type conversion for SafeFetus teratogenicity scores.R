setwd("~/Desktop")
library(stringr)
options(stringsAsFactors = FALSE)
library(readxl)

drugbank.safefetus.pregrisk.pfam.table <- read.csv("~/Desktop/drugbank.safefetus.pregrisk.pfam.table.txt")
View(drugbank.safefetus.pregrisk.pfam.table)
write.csv(drugbank.safefetus.pregrisk.pfam.table, file='risk.csv')