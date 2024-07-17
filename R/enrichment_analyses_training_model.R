#check the perform gene ontology  using enrichGO
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
#Load the functions most regularly used
source("R/Trainome_functions.R")

lncRNAS <-  lncRNAS <- readRDS("data/lncRNA_genes.RDS")

ct_metadata <- readRDS("data/contratrain_metadata.RDS")
#Load the four different models

#The volume model with set 3 as baseline
train <- readRDS("data/seqwrap_generated_models/training_coexpression_models/filtered_trained_untrained_midexc_correlation.RDS")

unique(train$coef)

colnames(train)

length(unique(train$geneid))
#how many correlated protein-coding genes


x_counts <- train %>%
  subset(coef == "counts")

length(unique(x_counts$geneid))


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
