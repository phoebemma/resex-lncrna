source("R/Trainome_functions.R")


#This is to test the linear mixed effects model on a sample lncra for coexpression analyses
#Load the metadata
ct_metadata <- readRDS("data/contratrain_metadata.RDS")





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
#Limiting the log fold 2 change to those above 1, or those below -1

trained_t4 <- train_model %>%
  dplyr::filter(coef == "training_statustrained:timet4")%>%
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
  rename(lncRNA = gene_name)


#initialising the arguments
model <- lmer(lncRNA~counts + time*condition +(1|participant),
            data = met_df)


