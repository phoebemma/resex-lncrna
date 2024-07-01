#check the perform gene ontology  using enrichGO
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyverse)






#Load the filtered coef data of the model

#This is done to merge i with the model eveluation and filter only models that pass the evaluation parameters

str_cor_t3 <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/filtered_less_stringent_simple_correlation_model.RDS")



length(unique(str_cor_t3$target))

#laod the full data for universe
genes <- readRDS("data/protein_coding_genes_FPKM.RDS")



ego_df <- enrichGO(gene = unique(simp_mod_sum$target),
                   universe = genes$gene_name,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)


### It is different when the gene expression data is used as universe, versus when it isnt
cluster_summary <- data.frame(ego_df)

dotplot(ego_df, showCategory = 15,
        
        font.size = 8, title = "15 top ranked biological processes in coexpressed protein-coding genes between trained and untrained participants at midexercise") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




length(unique(str_cor_t3$target))




#get the entrezid of the unique genes
entrez_ids <- bitr(sim_cor_t3$target, "SYMBOL", "ENTREZID", org.Hs.eg.db)


#pathway overrepresentation analyses
kegg_df <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)

#kegg_summary <- data.frame(kegg_df)


barplot(kegg_df, showCategory = 30, title = "30 highest enriched pathways in co-expressed protein coding genes between trained and untrained individuals at midexercise")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 9), axis.title.x = element_text(size = 20))



#keggs module overrepresentation analyses
over_df <- enrichMKEGG(entrez_ids$ENTREZID,
                       organism = "hsa",
                       keyType = "kegg",
                       # OrgDb = org.Hs.eg.db, 
                       #ont = "MF", 
                       pAdjustMethod = "BH", 
                       qvalueCutoff = 0.05)



barplot(over_df, showCategory = 30, title = " highest enriched pathways in co-expressed protein coding genes at baseline")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 9), axis.title.x = element_text(size = 20))









#This is the coexpression model
less_str_cor_t3 <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/filtered_less_stringent_simple_correlation_model.RDS")



int_cor_t3 <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/filtered_interaction_training_correlation_model.RDS")



#t4 data
less_str_cor_t4 <- readRDS("data/models/seqwrap_generated_models/filtered_t4_simple_less_string_correlation.RDS")

# less_str_cor_t4$summaries[[1]]
# 
# 
# temp <-  less_str_cor_t4$errors%>%
#   
#   mutate(err = unlist(errors_fit)) %>%
#   
#   
#   
#   pivot_longer(cols = errors_fit:warn_eval) %>%
#   
#   filter(name == "err_sum") %>%
#   print()
# 
# unlist(temp$value)
# 
# x <- bind_rows(within(less_str_cor_t4$summaries, rm(LTA))) %>%
#   subset(!coef == "(Intercept)") %>%
#   mutate(target = rep(names(within(less_str_cor_t4$summaries, rm(LTA))), each = 37))%>%
#   
#   mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
#          log2fc = Estimate/log(2),
#          
#          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
#   filter(fcthreshold == "s" & adj.p <= 0.05 )%>%
#   print()
# 
# saveRDS(x, "data/models/seqwrap_generated_models/filtered_t4_simple_less_string_correlation.RDS")


