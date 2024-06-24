
#read<- readr::read_csv("data/Contratrain_rsem_genes/1-subj1sample1.genes.results")

#This script builds the model based on training status of participants' legs

#It also filters the model parameters

#Load the functions most regularly used
source("R/Trainome_functions.R")



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




#argument for a second model that looks at trained versus untrained

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


library(marginaleffects)

#predictions() function calculates the regression-adjusted predicted values for every observation
pre <- predictions(training_model$models$`A1BG-AS1`, type = "link")




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


#training_model_filt <- readRDS("data/models/seqwrap_generated_models/filtered_training_model.RDS")





#Select those differentially expressed at midexercise
t3 <- training_model_filt %>%
  dplyr::filter(coef == "training_statustrained:timet3")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_trained_vs_untrained_at_midexercise.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t3, "DE lncs between trained and untrained legs at mid exercise")
dev.off()



#Select those differentially expressed at post exercise
t4 <- training_model_filt %>%
  dplyr::filter(coef == "training_statustrained:timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_trained_vs_untrained_postexercise.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t4, "DE lncs between trained and untrained at postexercise")
dev.off()




#Coexpression analyses

#This builds a model to check the protein coding genes coexpressed with the DE lncs

#The analyses below uses TPM values for both lncs and the proetin coding genes

library(ggrepel)
library(lme4)
library(lmerTest)





#Load the protein coding genes saved as TPM values
genes_TPM <- readRDS("data/protein_coding_genes_TPM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))




#Load the full gene counts in TPM
full_df <- readRDS("data/Ct_genes_TPM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))

#loda full counts in FPKM

genes_fpkm <- readRDS("data/Ct_genes_FPKM.RDS")%>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))


mRNA_genes_fpkm <- genes_fpkm[genes_fpkm$gene_name %in% genes_TPM$gene_name,]

mRNA_genes_fpkm<- mRNA_genes_fpkm %>%
  dplyr::filter(rowSums(mRNA_genes_fpkm[,-1]) != 0)


#extract the lncs of interest at mid exercise
lncs_of_int <- genes_fpkm[genes_fpkm$gene_name %in% t4$target,]



#extract the metadata and merge to the lncs of interest
met_df <- lncs_of_int %>%
  pivot_longer(cols = -("gene_name"),
               names_to = "seq_sample_id",
               values_to = "counts") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  #rename the gene_name to lncRNA to avoid mixing up with the mRNA genenames
  dplyr::rename(lncRNA = gene_name)




#initialising the arguments
args<- list(formula = y ~  lncRNA + lncRNA:counts + lncRNA:time + lncRNA:condition + (1|participant)) 
                                                                            
                                                                            
args1<- list(formula = y ~  counts + lncRNA + time + condition  + (1|participant))                                                                  

# Build the correlation model

cor_model <- seqwrap(fitting_fun = lmerTest::lmer,
                     arguments = args1,
                     data = mRNA_genes_fpkm,
                     metadata = met_df,
                     samplename = "seq_sample_id",
                     summary_fun = sum_fun_lmer,
                     eval_fun = eval_mod_lmer,
                     exported = list(),
                     save_models = FALSE,
                     return_models = FALSE,
                     
                     #subset = 1:10,
                     cores = ncores-2)

#saveRDS(cor_model, "data/models/seqwrap_generated_models/simpler_training_correlation_model.RDS")

length(unique(met_df$lncRNA))
#the one below was obtained from the less stringent filtering based on logfc changes
#saveRDS(cor_model, "data/models/seqwrap_generated_models/less_stringent_simpler_training_correlation_model.RDS")


#saveRDS(cor_model, "data/models/seqwrap_generated_models/t4_simple_less_stringent_correlation_model.RDS")
#cor_model <- readRDS("data/models/seqwrap_generated_models/simpler_training_correlation_model.RDS")

cor_model2 <- seqwrap(fitting_fun = lmerTest::lmer,
                     arguments = args,
                     data = mRNA_genes_fpkm,
                     metadata = met_df,
                     samplename = "seq_sample_id",
                     summary_fun = sum_fun_lmer,
                     eval_fun = eval_mod_lmer,
                     exported = list(),
                     save_models = FALSE,
                     return_models = FALSE,
                     
                     #subset = 1:10,
                     cores = ncores-2)



cor_model$summaries[[1]]

cor_model2$summaries[[1]]

#saveRDS(cor_model2, "data/models/seqwrap_generated_models/interaction_training_model.RDS")
#saveRDS(cor_model2, "data/models/seqwrap_generated_models/t4_interaction_less_stringent_training_model.RDS")
#cor_model$errors$err_sum
### Plotiing one gene vs lncRNA

genes_fpkm %>%
  filter(gene_name == "TSPAN6") %>%
  pivot_longer(cols = starts_with("X"), names_to = "seq_sample_id" ) %>%
  inner_join(met_df %>%
               filter(lncRNA == "LANCL1-AS1")) %>%
  ggplot(aes(value, counts)) + geom_point()
  print()



temp <-  cor_model$errors%>%
   
    mutate(err = unlist(errors_fit)) %>%

  
  
  pivot_longer(cols = errors_fit:warn_eval) %>%
    
    filter(name == "err_sum") %>%
    print()

unlist(temp$value)
print(x) 
bind_cols(x)

bind_rows(cor_model2$summaries) %>%
  
  filter(coef == "lncRNASMIM2-IT1") %>%
  # dim()
  
  
  #mutate(gene = names(cor_model$summaries)) %>%
  ggplot(aes(Pr...t..)) + geom_histogram()



length(names(cor_model$summaries))

cor_model$evaluations



#get model summary using the in-house for combinaing all model summaries into a table
#this takes as input the model name and number of coefficients


#This is the simpler model that takes the count, lncRNA, time and condition as input
simpler_model <- model_sum_lmer(cor_model, 61)

#saveRDS(simpler_model, "data/models/seqwrap_generated_models/filtered_simpler_training_model.RDS")

interaction_filled_model <- model_sum_lmer(cor_model2, 90)

cor



#GABRGI is a gene throwing errors in the model
#this removes it from the model summaries
  x <- bind_rows(within(cor_model$summaries, rm(GAL, RNASET2))) %>%
  subset(!coef == "(Intercept)") %>%
  mutate(target = rep(names(within(cor_model$summaries, rm(GAL, RNASET2))), each = 60))%>%

  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
         log2fc = Estimate/log(2),

         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )%>%
    print()

  
  unique(x$coef)
  
  #save this model
#saveRDS(x, "data/models/seqwrap_generated_models/filtered_interaction_training_correlation_model.RDS")
  

#saveRDS(x, "data/models/seqwrap_generated_models/trained_untrained_midexercise/filtered_less_stringent_simple_correlation_model.RDS")
# 
# x <- unlist(cor_model2$summaries$GABRG1)
# 
# length(rep(names(cor_model$summaries)))
# 
# names(cor_model$summaries)


x <- bind_rows(cor_model$summaries) %>%
  subset(!coef == "(Intercept)")



unique(simpler_model$coef)


unique(simpler_model$target)




#check the perform gene ontology  using enrichGO
library(clusterProfiler)
ego_df <- enrichGO(gene = x$target,
                   #universe = unique(prot_genes$gene_name),
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = T)


cluster_summary <- data.frame(ego_df)

dotplot(ego_df, showCategory = 15,

        font.size = 8, title = "15 top ranked biological processes in coexpressed protein-coding genes set6 participants at t4") +
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20))




length(unique(x$target))


library(org.Hs.eg.db)

#get the entrezid of the unique genes
entrez_ids <- bitr(x$target, "SYMBOL", "ENTREZID", org.Hs.eg.db)


#pathway overrepresentation analyses
kegg_df <- enrichKEGG(gene = entrez_ids$ENTREZID,
                      organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)

kegg_summary <- data.frame(kegg_df)


barplot(kegg_df, showCategory = 30, title = "30 highest enriched pathways in co-expressed protein coding genes at baseline")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 9), axis.title.x = element_text(size = 20))



#keggs module overrepresentation analyses
over_df <- enrichMKEGG(entrez_ids$ENTREZID,
                             organism = "hsa",
                      keyType = "kegg",
                      # OrgDb = org.Hs.eg.db, 
                      #ont = "MF", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05)



barplot(over_df, showCategory = 30, title = "30 highest enriched pathways in co-expressed protein coding genes at baseline")+
  theme(axis.text = element_text(size = 15), axis.text.y = element_text(size = 9), axis.title.x = element_text(size = 20))








# #Load the raw model
# train_all <- readRDS("data/models/training_model_all.RDS")
# 
# 
# 
# #get model evaluation
# mod_eval <- model_eval(train_all)
# 
# #get the model summaries
# mod_sum <- model_sum(train_all, 7)
# 
# 
# #merge the model summaries and model evaluation into one dataframe
# 
# 
# train_all <- filt_model_parameters(mod_sum, mod_eval) %>%
#   dplyr::filter(coef != "efflibsize")
# 
# unique(train_all$coef)
# 
# 
# #save the filtered coefs
# #saveRDS(train_all, file = "./data/models/Filtered_coefs/training_model_all.RDS")
# 
# 
# 
# 
# #Load the training model nrmalised with lncs alone
# 
# train_lncs <- readRDS("data/models/training_model_lncRNA.RDS")
# 
# 
# #get model evaluation for lnsc
# mod_eval_lncs <- model_eval(train_lncs)
# 
# #get model summaries
# mod_sum_lncs <- model_sum(train_lncs, 7)
# 
# 
# train_lncs <- filt_model_parameters(mod_sum_lncs, mod_eval_lncs) %>%
#   dplyr::filter(coef != "efflibsize_lncs")
# 
# #saveRDS(train_lncs, file = "./data/models/Filtered_coefs/training_model_lncs.RDS")
# 




















#Select those differentially expressed at time t4
t4 <- Vol_all %>%
  dplyr::filter(coef == "timet4")

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_training_model_lnscs_model_at_t4.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t4, "DE lncs at time t4 when lncs were used in normalisation")
dev.off()



#Select those differentially expressed at time t3
t3 <- Vol_all %>%
  dplyr::filter(coef == "timet3")

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_lncs_model_at_t3.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(t3, "DE lncs at time t3 when lncs were used in normalisation")
dev.off()





#Select those differentially expressed at trained
set6 <- Vol_all %>%
  dplyr::filter(coef == "training_statustrained")

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_lncs_model_at_trained.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(set6, "DE lncs of strained all were used in normalisation")
dev.off()


#Select those differentially expressed at set3 time t4
set3_time4 <- Vol_all %>%
  dplyr::filter(coef == "training_statustrained:timet3")

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_all_model_at_trained_timet3.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(set3_time4, "DE lncs of trained at time t4  when all  were used in normalisation")
dev.off()


#Select those differentially expressed at set6 time t4
set6_time4 <- Vol_all %>%
  dplyr::filter(coef == "conditionset6:timet4")

#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_lncs_model_at_set6timet4.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(set6_time4, "DE lncs of set 6 at time t4  when lncs were used in normalisation")
dev.off()
