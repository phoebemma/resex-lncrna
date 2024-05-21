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


#Load the full gene counts in TPM
full_df <- readRDS("data/Ct_genes_TPM.RDS")



# Let us look at the DEs between the trained and untrained at post exercise

trained_t4 <- train_model %>%
  dplyr::filter(coef == "training_statustrained:timet4")



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
  rename(lncRNA = gene_name)





#initialising the arguments
args<- list(formula = y ~  lncRNA + training_status +(1|participant),
            family = glmmTMB::nbinom2())


# Build the correlation model

cor_model <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args,
                         data = genes_TPM,
                         metadata = met_df,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         #return_models = F,
                         subset = 1:550,
                         cores = ncores)









