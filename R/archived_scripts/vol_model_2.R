# The first volume model is based on condition set 0 as baseline
#This alternative model will use set3 as baseline and compare the other conditions to it. 

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

#load lncRNAs counts

lncRNAS <- readRDS("data/lncRNA_genes.RDS")

#Load metadata file

ct_metadata <- readRDS("data/contratrain_metadata.RDS")

#reorder the levels of the metadata such that set3 becomes the baseline

ct_metadata_reordered<- ct_metadata %>%
  mutate(condition = factor(condition, levels = c("set3", "set6", "set0")))

#argument for model that looks at the difference between conditions over time
args<- list(formula = y ~  efflibsize + condition*time +(1|participant),
            family = glmmTMB::nbinom2())



#model
volume_model<- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                       arguments = args,
                       data = lncRNAS,
                       metadata = ct_metadata_reordered,
                       samplename = "seq_sample_id",
                       summary_fun = sum_fun,
                       eval_fun = eval_mod,
                       exported = list(),
                       save_models = FALSE,
                       return_models = FALSE,
                       cores = ncores-2)


#saveRDS(volume_model, "data/seqwrap_generated_models/volume_model_withset3_baseline.RDS")
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

unique(volume_model_filt$coef)

cond6_t3 <- volume_model_filt %>%
  dplyr::filter(coef == "conditionset6:timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)


#make a volcano plot using the plot_volcano function
jpeg(filename = "./plots/DE_set6_with_set3_baseline_at_t4.jpeg",
     width = 850, height = 500, quality = 100)
plot_volcano(cond6_t3, "DE lncs between set 6 and set3 at postexercise")
dev.off()




#Extrat the DE lncs at mid exercise 

set6_t3_lncs <- lncRNAS[lncRNAS$gene_name %in% cond6_t3$target,]




set6_t3_lncs <- pivot_longer(data = set6_t3_lncs,
                        cols = -(gene_name),
                        names_to = "seq_sample_id",
                        values_to = "counts")






#Merge the lncs to the metadata

set6_t3_lncs <- set6_t3_lncs %>%
  inner_join(ct_metadata %>%
               filter(condition != "set0"), by = "seq_sample_id")

set6_t3_lncs$log_counts <- log(set6_t3_lncs$counts)

ggplot(data = set6_t3_lncs, 
             mapping = aes(x = time,
                           y = gene_name,
                           fill = log_counts
                           )) +
  geom_tile() +
 facet_grid(~ condition)


