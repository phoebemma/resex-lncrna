#| #Load the needed libraries
source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")


#Load the raw model
Vol_all <- readRDS("data/models/vol_model_all.RDS") 

 mod_eval <-  bind_rows(Vol_all$model_evaluations) %>%
  mutate(target = names(Vol_all$model_evaluations)) %>%
  print()
#get the model evaluation

 
 
 
 #get the model summaries
 
 model_sum <- bind_rows(Vol_all$model_summarises) %>%
   mutate(target = rep(names(Vol_all$model_summarises), each = 10))%>%
   subset(!coef == "(Intercept)") %>%
   mutate(adj.p = p.adjust(Pr...z.., method = "fdr"),
          log2fc = Estimate/log(2),
          
          fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))
 
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
Vol_all <- inner_join(model_sum, mod_eval, by = "target") %>%
 dplyr::filter(Pr...z.. <= 0.05 & fcthreshold == "s"  & pval.disp >= 0.05 
               & pval.unif >= 0.05 & coef != "efflibsize" ) %>%
  
   print()
 
 unique(Vol_all$coef)
 
 
 #save the filtered coefs 
 #saveRDS(Vol_all, file = "./data/models/Filtered_coefs/Vol_model_all_coefs.RDS")
 
 
 
 #Select those differentially expressed at time t4
 t4 <- Vol_all %>%
   dplyr::filter(coef == "timet4")
 
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
 
 