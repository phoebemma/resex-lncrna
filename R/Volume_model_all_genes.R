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
Vol_all <- inner_join(model_sum, mod_eval, by = "target") #%>%
 dplyr::filter(Pr...z.. <= 0.05 & fcthreshold == "s"  & pval.disp >= 0.05 
               & pval.unif >= 0.05) %>%
   print()
 
 unique(Vol_all$coef)
 