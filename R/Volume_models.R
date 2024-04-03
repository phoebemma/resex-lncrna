#| #Load the needed libraries
source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")


#Load the raw model
Vol_all <- readRDS("data/models/vol_model_all.RDS") 


#get the model evaluation
mod_eval <-  model_eval(Vol_all)
 
 
 
 #get the model summaries
 
 mod_sum <- model_sum(Vol_all, 10)
 
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
 
 