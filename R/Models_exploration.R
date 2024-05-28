# This script contains exploratory analyses of the different models
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
#Load the functions most regularly used
source("R/Trainome_functions.R")



ct_metadata <- readRDS("data/contratrain_metadata.RDS")
#Load the four different models

#The volume model normalised with all genes
Vol_all <- readRDS("data/models/Filtered_coefs/Vol_model_all.RDS")

unique(Vol_all$coef)

#The volume model normalised with only lncs

Vol_lncs <- readRDS("data/models/Filtered_coefs/Vol_model_lncs.RDS")



#The training model normalised with all genes

train_all <- readRDS("data/models/Filtered_coefs/training_model_all.RDS")

unique(train_all$coef)
#The training model normalised with only lncs

train_lncs <- readRDS("data/models/Filtered_coefs/training_model_lncs.RDS")






#Compare DE lncs at time t4 in all 4 models
t4_vol_all <- Vol_all %>%
  dplyr::filter(coef == "timet4") %>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)



#Plot it using the  plot_volcano function
t4_vol_all_plot <- plot_volcano(t4_vol_all, "DE lncRNAs post exercise when normalised using all genes")


#Reapeat same for volume_lncs

t4_vol_lncs <- Vol_lncs %>%
  dplyr::filter(coef == "timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

t4_vol_lncs_plot <- plot_volcano(t4_vol_lncs, "DE lncRNAs post exercise when normalised using only lncRNAs")


grid.arrange(t4_vol_all_plot, t4_vol_lncs_plot, ncol = 2)



#Repeat for the training models

t4_train_all <- train_all %>%
  dplyr::filter(coef == "timet4")


#Plot it using the  plot_volcano function
t4_train_all_plot <- plot_volcano(t4_train_all, "DE lncRNAs post exercise when normalised using all genes (training model)")




#plot t4 for both volume and training model
grid.arrange(t4_vol_all_plot, t4_train_all_plot, ncol = 2)


t4_train_lncs <- train_lncs %>%
  dplyr::filter(coef == "timet4")


#Plot it using the  plot_volcano function
t4_train_lnc_plot <- plot_volcano(t4_train_lncs, "DE lncRNAs post exercise when normalised using lncRNAs (training model)")


jpeg(filename = "./plots/DE_at t4.jpeg", 
     width = 1500, height = 1000, quality = 100)
grid.arrange(t4_train_all_plot, t4_train_lnc_plot,t4_vol_all_plot, t4_vol_lncs_plot)
dev.off()


#Look at trained versus untrained in the training models at t4

trained_at_t4_all <- train_all %>%
  dplyr::filter(coef == "training_statustrained:timet4")

trained_t4_all_plot <- plot_volcano(trained_at_t4_all, "DE lncRNAs trained individuals at post exercise (normalised with all genes)")


trained_t4_lncs <- train_lncs %>%
  dplyr::filter(coef == "training_statustrained:timet4")

trained_t4_lncs_plot <- plot_volcano(trained_t4_lncs, "DE lncRNAs trained individuals at post exercise (normalised with only lncs)")

grid.arrange(trained_t4_all_plot, trained_t4_lncs_plot, ncol = 2)







#Look at trained versus untrained in the training models at t3

trained_at_t3_all <- train_all %>%
  dplyr::filter(coef == "training_statustrained:timet3")

trained_t3_all_plot <- plot_volcano(trained_at_t3_all, "DE lncRNAs trained individuals at mid exercise (normalised with all genes)")


trained_t3_lncs <- train_lncs %>%
  dplyr::filter(coef == "training_statustrained:timet3")

trained_t3_lncs_plot <- plot_volcano(trained_t3_lncs, "DE lncRNAs trained individuals at mid exercise (normalised with only lncs)")


#Plot t3 and t4 in both models
jpeg(filename = "./plots/DE_trained_mid_and_post_exercise.jpeg", 
     width = 1500, height = 1000, quality = 100)
grid.arrange(trained_t4_all_plot, trained_t4_lncs_plot, trained_t3_all_plot, trained_t3_lncs_plot)
dev.off()


#Extract the DEs at post exercise for the set 6 individuals
Vol_all_set6_t4 <- Vol_all %>%
  dplyr::filter(coef == "conditionset6:timet4")



Vol_lncs_set6_t4 <- Vol_lncs %>%
  dplyr::filter(coef == "conditionset6:timet4")


#Extract the DEs at post exercise for the set 3 individuals
Vol_all_set3_t4 <- Vol_all %>%
  dplyr::filter(coef == "conditionset3:timet4")



Vol_lncs_set3_t4 <- Vol_lncs %>%
  dplyr::filter(coef == "conditionset3:timet4")




Vol_lncs_set3_t3<- Vol_lncs %>%
  dplyr::filter(coef == "conditionset3:timet3")



Vol_lncs_set6_t3 <- Vol_lncs %>%
  dplyr::filter(coef == "conditionset6:timet3")



Vol_all_set6_t4_plot <- plot_volcano(Vol_all_set6_t4, "DE lncRNAs set6 at post exercise (normalised with all genes)")

Vol_lncs_set6_t4_plot <- plot_volcano(Vol_lncs_set6_t4, "DE lncRNAs set6 at post exercise (normalised with lncs only)")

Vol_all_set3_t4_plot <- plot_volcano(Vol_all_set3_t4, "DE lncRNAs set3 at post exercise (normalised with all genes)")

Vol_lncs_set3_t4_plot <- plot_volcano(Vol_lncs_set3_t4, "DE lncRNAs set3 at post exercise (normalised with lncs)")

Vol_lncs_set3_t3_plot <- plot_volcano(Vol_lncs_set3_t3, "DE lncRNAs set3 at mid exercise (normalised with lncs)")

Vol_lncs_set6_t3_plot <- plot_volcano(Vol_lncs_set6_t3, "DE lncRNAs set6 at mid exercise (normalised with lncs only)")


#Plot t3 and t4 in both models
jpeg(filename = "./plots/DE_Volume_post_and_mid_exercise.jpeg", 
     width = 1500, height = 1000, quality = 100)
grid.arrange( Vol_lncs_set6_t4_plot, Vol_lncs_set3_t4_plot,Vol_lncs_set6_t3_plot, Vol_lncs_set3_t3_plot )
dev.off()





#Load the lncRNA dataset

#lncs <- readRDS("data/lncRNA_genes.RDS")

# lncs <-readRDS("data/Ct_genes_TPM.RDS")
# #Extrat the DE lncs at mid exercise for set 6
# 
# lncs_of_int <- lncs[lncs$gene_name %in% Vol_lncs_set3_t4$target,]
# 
# 
# 
# 
# 
# # lncs_ <- as.matrix(lncs_of_int[, -1])
# # rownames(lncs_) <- lncs_of_int$gene_name
# # 
# # lncs_ <- as.data.frame(lncs_)
# # 
# # cem <- cemitool(lncs_)
# 
# 
# 
# #Load protein cpding genes
# prot_genes <- readRDS("data/protein_coding_genes.RDS")
# 
# #Load the proteins of interest from the TPM data
# 
# mRNAs <- lncs[lncs$gene_name %in% prot_genes$gene_name, ]
# 
# 
# #saveRDS(mRNAs, file = "data/protein_coding_genes_TPM.RDS")
# 
# 
# 
# 
# #extract the metadata and merge to the lncs of interest
# met_df <- lncs_of_int %>%
#   pivot_longer(cols = -("gene_name"),
#                names_to = "seq_sample_id",
#                values_to = "counts") %>%
#   inner_join(ct_metadata, by = "seq_sample_id") %>%
#   #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
#    rename(lncRNA = gene_name)
# 
# 
# 
# #Model building to check for coexpression between the lncs and mRNA genes
# 
# #initialising the arguments
# args<- list(formula = y ~  lncRNA + condition*time +(1|participant),
#             family = glmmTMB::nbinom2())
# 
# 
# # Build the correlation model
# 
# cor_model <- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
#                          arguments = args,
#                          data = mRNAs,
#                          metadata = met_df,
#                          samplename = "seq_sample_id",
#                          summary_fun = sum_fun,
#                          eval_fun = eval_mod,
#                          exported = list(),
#                          #return_models = F,
#                          subset = 1:550,
#                          cores = ncores)
# 
# 
# 
# 
# 
# 
# #get the model evaluation
# mod_eval <-  model_eval(cor_model)
# 
# 
# unique(cor_model$model_summarises)
# 
# 
# mod_sum <- model_sum(cor_model, 43)
# 
# 
# 
# 
# 
# 
# 
# #merge the model summaries and model evaluation into one dataframe
# cor_model <- filt_model_parameters(mod_eval, mod_sum) 
# 
# unique(cor_model$coef)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Correlation test
# cor_results_j <-list()
# cor_results_i <- list()
# 
# for (i in 1:length(lncs_of_int$gene_name)){
#   
#   
#   for (j in 1:length(mRNAs$gene_name)){
#     
#     lncRNA_Vector <- as.numeric(lncs_of_int[i,-1])
#     protein_coding_vector <- as.numeric(mRNAs[j,-1 ])
#     
#     cor_test_result <- cor.test(lncRNA_Vector, protein_coding_vector, method = "spearman")
#     # cor_results <- do.call(rbind.data.frame(cor_test_result))
#     
#     cor_results_j[[j]] <- data.frame(
#       lncRNA = as.character(lncs_of_int[i,1]),
#       protein_coding_gene = as.character(mRNAs[j, 1]),
#       correlation_coefficient = cor_test_result$estimate,
#       p_value = cor_test_result$p.value, 
#       row.names = NULL)
#     
#     
#     
#   } 
#   
#   cor_results_i[[i]] <-  bind_rows(cor_results_j)
#   
# }
# 
# correlation_results <- bind_rows(cor_results_i) %>%
# #ggplot(aes(correlation_coefficient)) + geom_histogram()
# dplyr::filter(p_value <= 0.05) %>%
#   round_df(3)
# 
# 
# 
# 
# 
# #check the perform gene ontology  using enrichGO
# 
# ego_df <- enrichGO(gene = x$gene_name,
#                    #universe = unique(prot_genes$gene_name),
#                    keyType = "SYMBOL",
#                    OrgDb = org.Hs.eg.db, 
#                    ont = "BP", 
#                    pAdjustMethod = "BH", 
#                    qvalueCutoff = 0.05, 
#                    readable = T)
# 
# 
# 
# ## Output results from GO analysis to a table
# cluster_summary <- data.frame(ego_df)
# 
# dotplot(ego_df, showCategory = 15,
#         
#         font.size = 8, title = "15 top ranked biological processes in coexpressed protein-coding genes set6 participants at t4") +
#   theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))
