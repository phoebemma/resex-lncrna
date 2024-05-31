#This script is where we build a model that checks for coexpression of DE lncs
# and mRNA (protein_coding) genes



#source("R/archived_scripts/libraries.R")


library(dplyr)
library(trainomeHelper)
library(ggplot2)
library(ggrepel)
library(lme4)
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
genes_TPM <- readRDS("data/protein_coding_genes_TPM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))




#Load the full gene counts in TPM
full_df <- readRDS("data/Ct_genes_TPM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))



# Let us look at the DEs between the trained and untrained at post exercise
#Limiting the log fold 2 change to those above 1, or those below -1

trained_t4 <- train_model %>%
  dplyr::filter(coef == "timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)



#Visualize it
plot_volcano(trained_t4, "t")



# extract the DEs from the full gene counts


lncs_of_int <- full_df[full_df$gene_name %in% trained_t4$target,]



#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)





#initialising the arguments
args<- list(formula = y ~ counts + lncRNA + time  +(1|participant))


# Build the correlation model

cor_model <- seqwrap(fitting_fun = lme4::lmer,
                         arguments = args,
                         data = genes_TPM,
                         metadata = met_df,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         #return_models = F,
                        # subset = 1:550,
                         cores = ncores)



cor_model$su

#get model evaluation using the in-house for combinaing all model evaluations into a table
mod_eval <- model_eval(cor_model)



#get the model summary using the created function
#it takes as input the model name and number of unique coefficients
mod_sum <-  model_sum(vol_model_all, 10)






