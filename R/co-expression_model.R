#This script is where we build a model that checks for coexpression of DE lncs
# and mRNA (protein_coding) genes



#source("R/archived_scripts/libraries.R")


library(dplyr)
library(ggplot2)
library(ggrepel)
library(lme4)
library(lmerTest)
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
  dplyr::filter(coef == "training_statustrained:timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)



#Visualize it
plot_volcano(trained_t4, "t")



#Look at volume model

cond3_t4 <- Vol_model %>%
  dplyr::filter(coef == "conditionset3:timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


#Visualize it
plot_volcano(cond3_t4, "t")






# extract the DEs from the full gene counts


lncs_of_int <- full_df[full_df$gene_name %in% cond3_t4$target,]



#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)



#dim(genes_TPM)

#initialising the arguments
args<- list(formula = log(y + 0.1) ~ counts + lncRNA + time + condition  + (1|participant))


# Build the correlation model

cor_model <- seqwrap(fitting_fun = lmerTest::lmer,
                         arguments = args,
                         data = genes_TPM,
                         metadata = met_df,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun_lmer,
                         eval_fun = NULL,
                         exported = list(),
                     save_models = FALSE,
                     return_models = FALSE,
                    
                        # subset = 1:100,
                         cores = ncores)


cor_model$summaries[[1]]


cor_model$errors%>%
  mutate(err = unlist(errors_fit)) %>%
  print()
  
  
  pivot_longer(cols = errors_fit:warn_eval) %>%
  
  filter(!is.null(unlist(value))) %>%
  print()

length(names(cor_model$summaries))

bind_rows(cor_model$summaries) %>%
  
  filter(coef == "counts") %>%
 # dim()
  
  
 #mutate(gene = names(cor_model$summaries)) %>%
  ggplot(aes(Pr...t..)) + geom_histogram()

length(names(cor_model$summaries))

x <- cor_model$models[[1]]

plot(x)
eval_mod(x)

cor_model$summaries[[1]]
cor_model$evaluations


#get model evaluation using the in-house for combinaing all model evaluations into a table
mod_eval <- model_eval(cor_model)



#get the model summary using the created function
#it takes as input the model name and number of unique coefficients
mod_sum <-  model_sum(cor_model, 15)


length(names(cor_model$summaries))



nrow(cor_model$summaries)

x <- bind_rows(cor_model$summaries) %>%
 
  mutate(target = rep(names(cor_model$summaries), each = 15))%>%
  subset(!coef == "(Intercept)") %>%
  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
         log2fc = Estimate/log(2),
         
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))


length(unique(x$coef))
n(x)
