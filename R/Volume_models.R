#| #Load the needed libraries


#Load the functions most regularly used
source("R/Trainome_functions.R")


library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(marginaleffects)

#load lncRNAs counts

lncRNAS <- readRDS("data/lncRNA_genes.RDS")

#Load metadata file

ct_metadata <- readRDS("data/contratrain_metadata.RDS")











#argument for model that looks at the difference between conditions over time
args<- list(formula = y ~  efflibsize + condition*time +(1|participant),
            family = glmmTMB::nbinom2())



#model
volume_model<- seqwrap(fitting_fun = glmmTMB::glmmTMB,
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



#save model in data folder
#saveRDS(volume_model, "data/models/seqwrap_generated_models/volume_model.RDS")

#get model evaluation using the in-house for combinaing all model evaluations into a table
mod_eval <- model_eval(volume_model)


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


#saveRDS(volume_model_filt, "data/models/seqwrap_generated_models/filtered_volume_model.RDS")



volume_model_filt <- readRDS("data/models/seqwrap_generated_models/filtered_volume_model.RDS")



ct_metadata <- readRDS("data/contratrain_metadata.RDS")


unique(volume_model_filt$coef)

#extract those differentially expressed by condition 3 at midexercise
cond3_t3 <- volume_model_filt %>%
  dplyr::filter(coef == "conditionset6:timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


#load full dataset in fpkm
genes_fpkm <- readRDS("data/Ct_genes_FPKM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))

#load the mRNA data
mRNA_genes_fpkm <- readRDS("data/protein_coding_genes_FPKM.RDS")

mRNA_genes_fpkm<- mRNA_genes_fpkm %>%
  dplyr::filter(rowSums(mRNA_genes_fpkm[,-1]) != 0) %>%
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))


#extract the lncs of interest at mid exercise
lncs_of_int <- genes_fpkm[genes_fpkm$gene_name %in% cond3_t3$target,]

library(ggrepel)
library(lme4)
library(lmerTest)



#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)


#initialise argument
args <- list(formula = y ~  counts + time + condition  + (1|participant))    

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
                           
                           subset = 1:10,
                           cores = ncores-2)
  
  # Add names for each 
  
  geneids <- names(vol_cor_model$summaries)
  
  evaluations[[i]] <- bind_rows(vol_cor_model$evaluations) %>%
    mutate(geneid = geneids) %>%
    mutate(lncRNA = LR[i])
  
  
  summary_results[[i]] <- bind_rows(vol_cor_model$summaries) %>%
    mutate(geneid = rep(geneids, each = 6)) %>%
    mutate(lncRNA = LR[i])
  

  
}

bind_rows(summary_results)



# Build the correlation model



#saveRDS(vol_cor_model, "data/seqwrap_generated_models/conditions_models/conditions_set3_correlation.RDS")
#saveRDS(vol_cor_model, "data/seqwrap_generated_models/conditions_models/condition_set6_t3_correlation.RDS")
vol_cor_model$summaries[[1]]
vol_cor_model <- readRDS("data/seqwrap_generated_models/conditions_models/conditions_set3_t3_correlation.RDS")





#determine the models that raised errors
temp <-  vol_cor_model$errors%>%
  
  mutate(err = unlist(errors_fit)) %>%
  
  
  
  pivot_longer(cols = errors_fit:warn_eval) %>%
  
  filter(name == "err_sum") %>%
  print()

unlist(temp$value)



#Extract the models that gave null values in midexercise condition 3

#bind the model evaluations in one row, excluding those with null values
vol_mod_ev <- bind_rows(within(vol_cor_model$evaluations, rm(ADM, BARHL2, BFSP2, C1orf141, CFHR5, 
                                                             CST1, EGR4, GMPR2,  H3C12,IREB2, KRT85, 
                                                             MAPK8, MROH5,  OR1A1,  OR51M1, OR51T1, 
                                                             OR7E24, OR8G3P, PDHA2, PRAMEF1, PSG9,SERPINA9,
                                                             SH2D1A,  RNASE9, TMEM190,ZAR1L, ZIC5)))%>%
  mutate(target = names(within(vol_cor_model$evaluations, rm(ADM, BARHL2, BFSP2, C1orf141, CFHR5, 
                                                             CST1, EGR4, GMPR2,  H3C12,IREB2, KRT85, 
                                                             MAPK8, MROH5,  OR1A1,  OR51M1, OR51T1, 
                                                             OR7E24, OR8G3P, PDHA2, PRAMEF1, PSG9,SERPINA9,
                                                             SH2D1A,  RNASE9, TMEM190,ZAR1L, ZIC5))))

#hist(vol_mod_ev $pval.unif, main = "distribution of p unif values log converted dependent y")


vol_mod_SUM<- bind_rows(within(vol_cor_model$summaries, rm(ADM, BARHL2, BFSP2, C1orf141, CFHR5, 
                                                           CST1, EGR4, GMPR2,  H3C12,IREB2, KRT85, 
                                                           MAPK8, MROH5,  OR1A1,  OR51M1, OR51T1, 
                                                           OR7E24, OR8G3P, PDHA2, PRAMEF1, PSG9,SERPINA9,
                                                           SH2D1A,  RNASE9, TMEM190,ZAR1L, ZIC5))) %>%
  subset(!coef == "(Intercept)") %>%
  mutate(target = rep(names(within(vol_cor_model$summaries, rm(ADM, BARHL2, BFSP2, C1orf141, CFHR5, 
                                                               CST1, EGR4, GMPR2,  H3C12,IREB2, KRT85, 
                                                               MAPK8, MROH5,  OR1A1,  OR51M1, OR51T1, 
                                                               OR7E24, OR8G3P, PDHA2, PRAMEF1, PSG9,SERPINA9,
                                                               SH2D1A,  RNASE9, TMEM190,ZAR1L, ZIC5))), each = 17)) %>%
  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
       log2fc = Estimate/log(2),

       fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )

#saveRDS(vol_mod_SUM, "data/seqwrap_generated_models/conditions_models/filtered_conditions_set3_t3_correlation.RDS")

hist(vol_mod_SUM$adj.p, main = "adjusted p vlaues distribution in coexpressed mRNAs set6 mid exercise")

ggplot(data = vol_mod_SUM, aes(x = coef)) +
  geom_bar()+
  ggtitle("distribution of coefs of coexpressed proteins at midexercise set6")

unique(vol_mod_SUM$coef)
lnc_x <- vol_mod_SUM  %>%
  filter(coef == "lncRNALINC00310")

time3 <- vol_mod_SUM  %>%
  filter(coef == "timet3")




ego_df <- enrichGO(gene = time3$target,
                   universe = genes_fpkm$gene_name,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)


### It is different when the gene expression data is used as universe, versus when it isnt
cluster_summary <- data.frame(ego_df)

dotplot(ego_df, showCategory = 15,
        
        font.size = 5, title = "15 top biological processes in  protein-coding genes coexpressed postexercise") +
  theme(axis.text = element_text(size = 9), axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 10))




length(unique(str_cor_t3$target))




#get the entrezid of the unique genes
entrez_ids <- bitr(time3$target, "SYMBOL", "ENTREZID", org.Hs.eg.db)


#pathway overrepresentation analyses
kegg_df <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)

#kegg_summary <- data.frame(kegg_df)


barplot(kegg_df, showCategory = 30, title = " enriched pathways in  protein coding genes postexercise")+
  theme(axis.text = element_text(size = 9), axis.text.y = element_text(size = 7), axis.title.x = element_text(size = 10))






















































# #Load the raw model
# Vol_all <- readRDS("data/models/vol_model_all.RDS") 
# 
# 
# #get the model evaluation
# mod_eval <-  model_eval(Vol_all)
#  
#  
#  
#  #get the model summaries
#  
#  mod_sum <- model_sum(Vol_all, 10)
#  
 #colnames(model_sum)
 
#  mod_eval %>%
#    ggplot(aes(pval.unif)) + geom_histogram() +
#    ggtitle("contratrain gene level model evaluation (normalised at all genes level)")+
#    theme(axis.text = element_text(size = 15), text = element_text(size = 15),
#          plot.title = element_text(hjust = 0.5))
#  
# mod_eval%>%
#    ggplot(aes(pval.disp)) + geom_histogram() +
#    ggtitle("contratrain gene level model evaluation (normalised at all genes level)")+
#    theme(axis.text = element_text(size = 15), text = element_text(size = 15),
#          plot.title = element_text(hjust = 0.5))
# # 
# # 
#  mod_eval %>%
#    ggplot(aes(pval.zinfl)) + geom_histogram() +
#    ggtitle("contratrain gene level model evaluation (normalised at all genes level)")+
#    theme(axis.text = element_text(size = 15), text = element_text(size = 15),
#          plot.title = element_text(hjust = 0.5))
#  
 
 
 
 
 
 
 #merge the model summaries and model evaluation into one dataframe
Vol_all <- filt_model_parameters(mod_eval, mod_sum) %>%
  dplyr::filter(coef != "efflibsize")
 
 unique(Vol_all$coef)
 
 
 #save the filtered coefs 
 #saveRDS(Vol_all, file = "./data/models/Filtered_coefs/Vol_model_all.RDS")
 
 
 
 
 
 
 
 
 
 
 
 
 
 #Load the model that was normalised based on only lncRNAs
 Vol_lncs <- readRDS("data/models/vol_model_lncRNAs.RDS")
 
 
 #get the model evaluation
 mod_eval_lncs <-  model_eval(Vol_lncs)
 
 
 
 #get the model summaries
 
 mod_sum_lncs<- model_sum(Vol_lncs, 10)
 
 
 
 Vol_lncs <- filt_model_parameters(mod_eval_lncs, mod_sum_lncs) %>%
   dplyr::filter(coef != "efflibsize_lncs")
 
 unique(Vol_lncs$coef)
 
 #save the filtered coefs 
  #saveRDS(Vol_lncs, file = "./data/models/Filtered_coefs/Vol_model_lncs.RDS")
 
 
 
 
 
 
 
 
 
 
 
 
 #Select those differentially expressed at time t4 for vol_all
 t4_all <- Vol_all %>%
   dplyr::filter(coef == "timet4")
 
 
 t4_all_plot <- plot_volcano(t4_all, "DE lncs at time t4 when all genes were used in normalisation")
 
 #select t4 for lncs 
 t4_lncs <- Vol_lncs %>%
   dplyr::filter(coef == "timet4")
 
 t4_lncs_plot <- plot_volcano(t4_lncs, "DE at time t4 when lncs used in normalisation")
 
 grid.arrange(t4_all_plot, t4_lncs_plot, ncol = 2)
 
 #make a volcano plot using the plot_volcano function
 jpeg(filename = "./plots/DE_all_model_at_t4.jpeg",
      width = 850, height = 500, quality = 100)
 plot_volcano(t4, "DE lncs at time t4 when all genes were used in normalisation")
 dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 #Select those differentially expressed at time t3
 t3 <- Vol_all %>%
   dplyr::filter(coef == "timet3")
 
 #make a volcano plot using the plot_volcano function
 jpeg(filename = "./plots/DE_all_model_at_t3.jpeg",
      width = 850, height = 500, quality = 100)
 plot_volcano(t3, "DE lncs at time t3 when all genes were used in normalisation")
 dev.off()
 
 
 
 
 
 #Select those differentially expressed at set6
 set6 <- Vol_all %>%
   dplyr::filter(coef == "conditionset6")
 
 #make a volcano plot using the plot_volcano function
 jpeg(filename = "./plots/DE_all_model_at_set6.jpeg",
      width = 850, height = 500, quality = 100)
 plot_volcano(set6, "DE lncs of set 6 when all genes were used in normalisation")
 dev.off()
 
 
 #Select those differentially expressed at set3
 set3 <- Vol_all %>%
   dplyr::filter(coef == "conditionset3")
 
 #make a volcano plot using the plot_volcano function
 jpeg(filename = "./plots/DE_all_model_at_set3.jpeg",
      width = 850, height = 500, quality = 100)
 plot_volcano(set3, "DE lncs of set 3 when all genes were used in normalisation")
 dev.off()
 
 