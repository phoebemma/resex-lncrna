#This script is where we build a model that checks for coexpression of DE lncs
# and mRNA (protein_coding) genes



source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")


#Load the metadata
ct_metadata <- readRDS("data/contratrain_metadata.RDS")




#Load the two models, one that modeled volume of PRET, the other modelled presen 
#The volume model normalised with all genes
Vol_model <- readRDS("data/models/Filtered_coefs/Vol_model_all.RDS")

unique(Vol_model$coef)




#Load the training model
train_model <- readRDS("data/models/Filtered_coefs/training_model_all.RDS")

unique(train_model$coef)



#Laod the gene expression data in TPMS

#The volume and training models were built on raw counts.But for co-expression studies
#The TPM will be used


#Load the protein coding genes saved as TPM values
genes_TPM <- readRDS("data/protein_coding_genes_TPM.RDS")


