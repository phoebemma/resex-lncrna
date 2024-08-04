 

#Load the functions most regularly used
source("R/Trainome_functions.R")


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
#load lncRNAs counts

lncRNAS <- readRDS("data/lncRNA_genes.RDS")

#Load metadata file

ct_metadata <- readRDS("data/contratrain_metadata.RDS")

#reorder the levels of the metadata such that set3 becomes the baseline

#This is used for the model for DE lncs between sets 3 and 6

ct_metadata_reordered<- ct_metadata %>%
  mutate(condition = factor(condition, levels = c("set3", "set6", "set0")))


ct_metadata_reordered_2<- ct_metadata %>%
  mutate(condition = factor(condition, levels = c("set6", "set3", "set0")))

#argument for model that looks at the difference between conditions over time
args<- list(formula = y ~  efflibsize + condition*time  +(1|participant),
            family = glmmTMB::nbinom2())



#model
volume_model<- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                       arguments = args,
                       data = lncRNAS,
                       metadata = ct_metadata_reordered_2,
                       samplename = "seq_sample_id",
                       summary_fun = sum_fun,
                       eval_fun = eval_mod,
                       exported = list(),
                       save_models = FALSE,
                       return_models = FALSE,
                       cores = ncores-2)


#saveRDS(volume_model, "data/seqwrap_generated_models/volume_model_withset3_baseline.RDS")

#saveRDS(volume_model, "data/seqwrap_generated_models/volume_model_withset6_baseline.RDS")
mod_eval <- model_eval(volume_model)

hist(mod_eval$pval.zinfl)

volume_model$summaries[[1]]
#get the model summary using the created function
#it takes as input the model name and number of unique coefficients
#it also filters out the model coefficient called "intercept"
#creates adjusted p values, log2 fold change and fcld change significant threshold
mod_sum <-  model_sum(volume_model, 10)


#use the filter model parameter function to filter the models that fit the following conditions
#(Pr...z.. <= 0.05 & fcthreshold == "s"  & pval.disp >= 0.05 & pval.unif >= 0.05 

#it takes as input the data containing the model summary, and that containing the model evaluation
volume_model_filt <- filt_model_parameters(mod_eval, mod_sum)%>%
  dplyr::filter(coef != "efflibsize")


#saveRDS(volume_model_filt, "data/seqwrap_generated_models/filtered_vol_model_withset3_baseline.RDS")

#saveRDS(volume_model_filt, "data/seqwrap_generated_models/filtered_vol_model_withset6_baseline.RDS")


volume_model_filt <- readRDS("data/seqwrap_generated_models/filtered_volume_model.RDS")
unique(volume_model_filt$coef)




cond6_t3 <- volume_model_filt %>%
  dplyr::filter(coef == "conditionset6:timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

cond6_t4 <- volume_model_filt %>%
  dplyr::filter(coef == "conditionset6:timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_set3_with_set6_baseline_at_t3.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(cond3_t3, "DE lncs between set 6 and set3 at midexercise")
dev.off()


#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_set3_with_set6_baseline_at_postExc.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(cond3_t4, "DE lncs between set 6 and set3 at postexercise")
dev.off()



#Extract the DE lncs of interest


lncs_of_int <- lncRNAS[lncRNAS$gene_name %in% cond6_t3$target,]




set6_t3_lncs <- pivot_longer(data = lncs_of_int,
                        cols = -(gene_name),
                        names_to = "seq_sample_id",
                        values_to = "counts")






#Merge the lncs to the metadata

# set6_t3_lncs <- set6_t3_lncs %>%
#   inner_join(ct_metadata %>%
#                filter(condition != "set0")
#              , 
#                by = "seq_sample_id")

set6_t3_lncs$log_counts <- log(set6_t3_lncs$counts)

plot <- ggplot(data = set6_t3_lncs, 
             mapping = aes(x = time,
                           y = gene_name,
                           fill = log_counts
                           )) +
  geom_tile() +
  #use the PuOr color palette from colorBrewer
  scale_fill_distiller(palette = "PuOr")+
  #scale_fill_gradient()+
 facet_grid(~ condition)+
  #set a base for all fonts
  theme_grey(base_size=8) +
  #add border white colour of line thickness 0.25
  geom_tile(colour="white", size=0.25)+
  ggtitle("Differentially expressed lncRNAs between Set6 and Set3 at Postexercise")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot, filename = "./plots/heatmap_DE_set6_vs_set3_at_postexercise_.jpeg",
       width = 8.5, height = 5, quality = 100)




##Coexpression analyses





#load the mRNA data
mRNA_genes_fpkm <- readRDS("data/protein_coding_genes_FPKM.RDS")

mRNA_genes_fpkm<- mRNA_genes_fpkm %>%
  dplyr::filter(rowSums(mRNA_genes_fpkm[,-1]) != 0) %>%
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))




#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)


#initialise argument
args <- list(formula = y ~  counts + time + condition  + sex + (1|participant))    


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
                           
                           #subset = 1:10,
                           cores = ncores-2)
  
  # Add names for each and remove those with null in output
  
  
  geneids <- names(which(vol_cor_model$summaries != "NULL"))
  
  #remove all genes with null output
  excl <- names(which( vol_cor_model$evaluations == "NULL" ))
  
  evaluations[[i]] <- bind_rows(within(vol_cor_model$evaluations, rm(excl))) %>%
    mutate(geneid = geneids) %>%
    mutate(lncRNA = LR[i])
  
  
  summary_results[[i]] <- bind_rows(within(vol_cor_model$summaries, rm(excl))) %>%
    mutate(geneid = rep(geneids, each = 7)) %>%
    mutate(lncRNA =LR[i])
  
}


x<- bind_rows(summary_results) %>%
  inner_join(bind_rows(evaluations),by = c("geneid", "lncRNA"))

hist(x$Pr...t..)
hist(x$pval.unif)

#saveRDS(x, "data/seqwrap_generated_models/conditions_models/int_model_set6_postexc_with_set3_as_baseline.RDS")
#saveRDS(x, "data/seqwrap_generated_models/volume_coexpression_models/set3_midexc_zero_baseline_coexpression_model.RDS")
#saveRDS(x, "data/seqwrap_generated_models/volume_coexpression_models/set3_postexc_zero_baseline_coexpression_model.RDS")

#saveRDS(x, "data/seqwrap_generated_models/volume_coexpression_models/set6_postexc_zero_baseline_coexpression_model.RDS")
#saveRDS(x, "data/seqwrap_generated_models/volume_coexpression_models/set6_midexc_zero_baseline_coexpression_model.RDS")


#filter the dataframe to the required parameters




colnames(x)
length(unique(x$coef))

x_filt <- x %>%
  subset(coef != "(Intercept)") %>%
  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
                                          log2fc = Estimate/log(2),
                                          
            fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))%>%
  filter(fcthreshold == "s" & adj.p <= 0.05 & pval.unif >= 0.05 )

#saveRDS(x_filt, "data/seqwrap_generated_models/conditions_models/filtered_int_model_set6_postexc_with_set3_baseline.RDS")
#saveRDS(x_filt, "data/seqwrap_generated_models/volume_coexpression_models/filtered_set3_midexc_zero_baseline_coexpression_model.RDS")
#saveRDS(x_filt, "data/seqwrap_generated_models/volume_coexpression_models/filtered_set3_postexc_zero_baseline_coexpression_model.RDS")

#saveRDS(x_filt, "data/seqwrap_generated_models/volume_coexpression_models/filtered_set6_postexc_zero_baseline_coexpression_model.RDS")

#saveRDS(x_filt, "data/seqwrap_generated_models/volume_coexpression_models/filtered_set6_midexc_zero_baseline_coexpression_model.RDS")


#filter those completely dependent on the counts
x_counts <- x_filt %>%
  subset(coef == "counts")

length(unique(x_counts$geneid))
#saveRDS(x_counts, "data/seqwrap_generated_models/volume_coexpression_models/set3_midexc_zero_base_correlation_counts.RDS")
#saveRDS(x_counts, "data/seqwrap_generated_models/volume_coexpression_models/set3_postexc_zero_base_correlation_counts.RDS")
#saveRDS(x_counts, "data/seqwrap_generated_models/volume_coexpression_models/set6_postexc_zero_base_correlation_counts.RDS")

#saveRDS(x_counts, "data/seqwrap_generated_models/volume_coexpression_models/set6_midexc_zero_base_correlation_counts.RDS")



#load genes data in fpkm
#This is the full dataset
genes_fpkm <- readRDS("data/Ct_genes_FPKM.RDS")


ego_df <- enrichGO(gene = x_counts$geneid,
                   universe = genes_fpkm$gene_name,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)

cluster_summary <- data.frame(ego_df)


jpeg(filename = "./plots/15_top_BP_set_6_at_postexc.jpeg",
     width = 850, height = 700, quality = 100)

dotplot(ego_df, showCategory = 15,
        
        font.size = 5, title = "15 top biological processes coexpressed proteins at set 6 postexercise") +
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
