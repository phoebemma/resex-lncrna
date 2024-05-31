source("R/Trainome_functions.R")


library(lme4)

library(glmmTMB)

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
  dplyr::rename(lncRNA = gene_name)


#initialising the arguments

model <- lmer(counts~  time+ condition +(1|participant),
            data = met_df)
testDispersion(model)
summary(model)





x <- data.frame(cbind(data.frame(coef = rownames(coef(summary(model)))),
                coef(summary(model)),
                row.names = NULL))
model





mod <- glmmTMB(counts~  time+condition +(1|participant),
               data = met_df, family = nbinom2() )
summary(mod)

coef(summary(mod))$cond




sum_fun <- function(x){
  
  cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$cond))),
                             coef(summary(x))$cond, 
                             
                             row.names = NULL)
  
  return(cond_effects)
  }