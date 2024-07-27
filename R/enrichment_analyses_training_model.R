#check the perform gene ontology  using enrichGO
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
#Load the functions most regularly used
source("R/Trainome_functions.R")

lncRNAS <-  lncRNAS <- readRDS("data/lncRNA_genes.RDS")

ct_metadata <- readRDS("data/contratrain_metadata.RDS")
#Load the four different models

#The volume model with set 3 as baseline
train_mid <- readRDS("data/seqwrap_generated_models/training_coexpression_models/filtered_trained_untrained_postexc_correlation.RDS")

unique(train_mid$coef)

colnames(train_mid)

length(unique(train_mid$geneid))
#how many correlated protein-coding genes


#The coefficient of interest is that related to the counts
x_counts <- train_mid %>%
  subset(coef == "counts")

length(unique(x_counts$geneid))


#save the correlation data
#saveRDS(x_counts, "data/seqwrap_generated_models/training_coexpression_models/mid_exc_correlation_counts.RDS")


#saveRDS(x_counts, "data/seqwrap_generated_models/training_coexpression_models/post_exc_correlation_counts.RDS")

#load genes data in fpkm
#This is the full dataset
genes_fpkm <- readRDS("data/Ct_genes_FPKM.RDS")


ego_df <- enrichGO(gene = x_counts$geneid,
                   universe = genes_fpkm$gene_name,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "cc",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)

cluster_summary <- data.frame(ego_df)


jpeg(filename = "./plots/15_top_cc_trained_untrained_at_midexc.jpeg",
     width = 850, height = 700, quality = 100)

dotplot(ego_df, showCategory = 15,
        
        font.size = 5, title = "15 top cellular components of coexpressed proteins at  midexercise") +
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

kegg_summary <- data.frame(kegg_df)

jpeg(filename = "./plots/15_most_enriched_pways_trained_untrained at_midexc.jpeg",
     width = 850, height = 700, quality = 100)
barplot(kegg_df, showCategory = 15, title = "10 most enriched pathways in co-expressed protein coding genes among trained indivuduals midexercise")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 13))
dev.off()










#Load the volume model that has set3 as baseline

Volume_model <- readRDS("data/seqwrap_generated_models/volume_coexpression_models/filtered_int_model_set6_midexc_with_set3_baseline.RDS")

unique(Volume_model$coef)



length(unique(Volume_model$geneid))
#how many correlated protein-coding genes


#The coefficient of interest is that related to the counts
vol_counts <- Volume_model %>%
  subset(coef == "counts")

length(unique(vol_counts$geneid))

#saveRDS(vol_counts, "data/seqwrap_generated_models/volume_coexpression_models/set6_midexc_counts.RDS")
