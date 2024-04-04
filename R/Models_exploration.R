# This script contains exploratory analyses of the different models
source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")


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
  dplyr::filter(coef == "timet4")


#Plot it using the  plot_volcano function
t4_vol_all_plot <- plot_volcano(t4_vol_all, "DE lncRNAs post exercise when normalised using all genes")


#Reapeat same for volume_lncs

t4_vol_lncs <- Vol_lncs %>%
  dplyr::filter(coef == "timet4")

t4_vol_lncs_plot <- plot_volcano(t4_vol_lncs, "DE lncRNAs post exercise when normalised using only lncRNAs")


grid.arrange(t4_vol_all_plot, t4_vol_lncs_plot, ncol = 2)



#Repeat for the training models

t4_train_all <- train_all %>%
  dplyr::filter(coef == "timet4")


#Plot it using the  plot_volcano function
t4_train_all_plot <- plot_volcano(t4_train_all, "DE lncRNAs post exercise when normalised using all genes (training model)")



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



