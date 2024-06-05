


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






















#Load the raw model
train_all <- readRDS("data/models/training_model_all.RDS") 



#get model evaluation
mod_eval <- model_eval(train_all)

#get the model summaries
mod_sum <- model_sum(train_all, 7)


#merge the model summaries and model evaluation into one dataframe


train_all <- filt_model_parameters(mod_sum, mod_eval) %>%
  dplyr::filter(coef != "efflibsize")

unique(train_all$coef)


#save the filtered coefs 
#saveRDS(train_all, file = "./data/models/Filtered_coefs/training_model_all.RDS")




#Load the training model nrmalised with lncs alone

train_lncs <- readRDS("data/models/training_model_lncRNA.RDS")


#get model evaluation for lnsc
mod_eval_lncs <- model_eval(train_lncs)

#get model summaries
mod_sum_lncs <- model_sum(train_lncs, 7)


train_lncs <- filt_model_parameters(mod_sum_lncs, mod_eval_lncs) %>%
  dplyr::filter(coef != "efflibsize_lncs")

#saveRDS(train_lncs, file = "./data/models/Filtered_coefs/training_model_lncs.RDS")





















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
