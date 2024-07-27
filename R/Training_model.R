#This script is for the analyses of the trained versus untrained individuals

#Load the functions most regularly used
source("R/Trainome_functions.R")

#load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(marginaleffects)
library(seqwrap)
library(glmmTMB)
library(RColorBrewer)
library(ggrepel)
library(lme4)
library(lmerTest)
library(gridExtra)

#load lncRNAs counts

lncRNAS <- readRDS("data/lncRNA_genes.RDS")

#Load metadata file

ct_metadata <- readRDS("data/contratrain_metadata.RDS")


#Description of the variables in the metadata dataframe


## lib.size = library size lncRNAs alone

## lib.size_all = library size all biotypes

## norm.factors = normalization factors lncRNAs alone

## norm.factors_all = normalization factors all biotypes

## efflibsize_lncs = effective library size of lncRNAs alone

## efflibsize = effective library size of all biotypes




#argument for a  model that looks at trained versus untrained

args<- list(formula = y ~  efflibsize + training_status*time +(1|participant),
            family = glmmTMB::nbinom2())





#Model for the trained versus untrained 

training_model<- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args,
                         data = lncRNAS,
                         metadata = ct_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         save_models = FALSE,
                         return_models = FALSE,
                         cores = ncores)




#sneekpeek
training_model$summaries$`A1BG-AS1`

training_model$errors

training_model$evaluations$`A1BG-AS1`




bind_rows(training_model$summaries) %>%
  
  filter(coef == "training_statustrained:timet4") %>%
  # dim()
  
  
  #mutate(gene = names(cor_model$summaries)) %>%
  ggplot(aes(Pr...z..)) + geom_histogram()






#save model in data folder
#saveRDS(training_model, "data/models/seqwrap_generated_models/training_model.RDS")

#get model evaluation using the in-house for combinaing all model evaluations into a table
mod_eval <- model_eval(training_model)



#get the model summary using the created function
#it takes as input the model name and number of unique coefficients
#it also filters out the model coefficient called "intercept"
#creates adjusted p values, log2 fold change and fcld change significant threshold
mod_sum <-  model_sum(training_model, 7)





#use the filter model parameter function to filter the models that fit the following conditions
#(Pr...z.. <= 0.05 & fcthreshold == "s"  & pval.disp >= 0.05 & pval.unif >= 0.05 

#it takes as input the data containing the model summary, and that containing the model evaluation


training_model_filt <- filt_model_parameters(mod_eval, mod_sum)%>%
  dplyr::filter(coef != "efflibsize")


unique(training_model_filt$coef)

#save filtered model
#saveRDS(training_model_filt, "data/models/seqwrap_generated_models/filtered_training_model.RDS")










##This is done when loading the already run model

#training_model_filt <- readRDS("data/seqwrap_generated_models/filtered_training_model.RDS")





#Select those differentially expressed at midexercise
t3_int <- training_model_filt %>%
  dplyr::filter(coef == "training_statustrained:timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


t3 <- training_model_filt %>%
  dplyr::filter(coef == "timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_trained_vs_untrained_at_midexercise.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t3, "DE lncs between trained and untrained legs at mid exercise")
dev.off()



#Select those differentially expressed at post exercise
t4 <- training_model_filt %>%
  dplyr::filter(coef == "timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

t4_int <- training_model_filt %>%
  dplyr::filter(coef == "training_statustrained:timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


t1 <- plot_volcano(t3, "DE lncRNAs at mid exercise")
t2 <- plot_volcano(t4, "DE lncRNAs at post exercise")
t1_int <- plot_volcano(t3_int, "DE lncRNAs between the trained and untrained at midexercise")
t2_int <- plot_volcano(t4_int, "DE lncRNAs between trained and untrained at post exercise")

grid.arrange(t1, t2, t1_int, t2_int)

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_trained_vs_untrained.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t4, "DE lncs between trained and untrained at postexercise")
dev.off()




#Coexpression analyses

#This builds a model to check the protein coding genes coexpressed with the DE lncs

#The analyses below uses TPM values for both lncs and the proetin coding genes






#load the mRNA data
mRNA_genes_fpkm <- readRDS("data/protein_coding_genes_FPKM.RDS")

mRNA_genes_fpkm<- mRNA_genes_fpkm %>%
  dplyr::filter(rowSums(mRNA_genes_fpkm[,-1]) != 0) %>%
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))




#extract the lncs of interest 
lncs_of_int <- lncRNAS[lncRNAS$gene_name %in% t4$target,]



#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)




length(unique(met_df$lncRNA))

#initialise argument
args <- list(formula = y ~  counts + time + training_status  + sex + (1|participant))    


#Prepare and run a loop using Seqwrap
#This loops through each individual lncRNA in the metadata,
summary_results <- list()
evaluations <- list()

LR <- unique(met_df$lncRNA)

for(i in seq_along(LR)) {
  
  vol_cor_model <- seqwrap(fitting_fun = lmerTest::lmer,
                           arguments = args,
                           data = mRNA_genes_fpkm,
                           metadata = filter(met_df, lncRNA == LR[i]),
                           samplename = "seq_sample_id",
                           summary_fun = sum_fun_lmer,
                           eval_fun = eval_mod_lmer,
                           exported = list(),
                           save_models = FALSE,
                           return_models = FALSE,
                           
                           #subset = 1:5,
                           cores = ncores-2)
  
  # Add names for each and remove those with null in output
  
  
  geneids <- names(which(vol_cor_model$summaries != "NULL"))
  
  #remove all genes with null output
  excl <- names(which( vol_cor_model$evaluations == "NULL" ))
  
  evaluations[[i]] <- bind_rows(within(vol_cor_model$evaluations, rm(excl))) %>%
    mutate(geneid = geneids) %>%
    mutate(lncRNA = LR[i])
  
  
  summary_results[[i]] <- bind_rows(within(vol_cor_model$summaries, rm(excl))) %>%
    mutate(geneid = rep(geneids, each = 6)) %>%
    mutate(lncRNA =LR[i])
  
  
  
}


x<- bind_rows(summary_results) %>%
  inner_join(bind_rows(evaluations),by = c("geneid", "lncRNA"))

hist(x$Pr...t..)
hist(x$pval.unif)

#saveRDS(x, "data/seqwrap_generated_models/training_coexpression_models/trained_untrained_post_exercise_correlation.RDS")




#filter the dataframe to the required parameters




colnames(x)
length(unique(x$coef))

unique(x$lncRNA)
x_filt <- x %>%
  subset(coef != "(Intercept)") %>%
  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
         log2fc = Estimate/log(2),
         
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))%>%
  filter(fcthreshold == "s" & adj.p <= 0.05 & pval.unif >= 0.05 )

#saveRDS(x_filt, "data/seqwrap_generated_models/training_coexpression_models/filtered_trained_untrained_postexc_correlation.RDS")
#filter those completely dependent on the counts
x_counts <- x_filt %>%
  subset(coef == "counts")

length(unique(x_counts$geneid))


#load genes data in fpkm
#This is the full dataset
genes_fpkm <- readRDS("data/Ct_genes_FPKM.RDS")


ego_df <- enrichGO(gene = x_counts$geneid,
                   universe = genes_fpkm$gene_name,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)

cluster_summary <- data.frame(ego_df)


jpeg(filename = "./plots/15_top_mf_trained_untrained_at_postexc.jpeg",
     width = 850, height = 700, quality = 100)

dotplot(ego_df, showCategory = 15,
        
        font.size = 5, title = "15 top Molecular function of coexpressed proteins at  postexercise") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 10))
dev.off()


#get the entrezid of the unique genes
entrez_ids <- bitr(x_counts$geneid, "SYMBOL", "ENTREZID", org.Hs.eg.db)


#pathway overrepresentation analyses
kegg_df <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)

#kegg_summary <- data.frame(kegg_df)

jpeg(filename = "./plots/15_most_enriched_pways_set_6_at_postexc.jpeg",
     width = 850, height = 700, quality = 100)
barplot(kegg_df, showCategory = 15, title = "10 most enriched pathways in co-expressed protein coding genes set6 at postexercise")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 13))
dev.off()
