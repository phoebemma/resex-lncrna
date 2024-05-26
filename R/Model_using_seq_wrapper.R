#Load the needed libraries

library(glmmTMB)
library(seqwrap)
#Load the functions most regularly used
source("R/Trainome_functions.R")


#load counts containing all biotypes
#all_biotypes_df <- readRDS("data/filtered_all_gene_counts.RDS")

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













#Define the model parametes


#Volume dependednt model parameters for all the biotypes
args<- list(formula = y ~  efflibsize + condition*time +(1|participant),
            family = glmmTMB::nbinom2())


#argument for a second model that looks at trained versus untrained in all biotypes

args_2<- list(formula = y ~  efflibsize + training_status*time +(1|participant),
              family = glmmTMB::nbinom2())


#arguements for a volume dependent model in lncs alone
args_3 <- list(formula = y ~  efflibsize_lncs + condition*time +(1|participant),
               family = glmmTMB::nbinom2())


#arguement for trained versus untraind in lncs alone

args_4<- list(formula = y ~  efflibsize_lncs + training_status*time +(1|participant),
              family = glmmTMB::nbinom2())




#Volume_dependent model
vol_model_all<- seqwrap(fitting_fun = glmmTMB::glmmTMB,
                       arguments = args,
                       data = lncRNAS,
                       metadata = ct_metadata,
                       samplename = "seq_sample_id",
                       summary_fun = sum_fun,
                       eval_fun = eval_mod,
                       exported = list(),
                       #return_models = F,
                       #subset = NULL,
                       cores = ncores)


#saveRDS(vol_model_all, file = "./data/models/vol_model_all.RDS")

#get model evaluation
mod_eval <- model_eval(vol_model_all)


#Model for the trained versus untrained 

training_model_all<- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                         arguments = args_2,
                         data = lncRNAS,
                         metadata = ct_metadata,
                         samplename = "seq_sample_id",
                         summary_fun = sum_fun,
                         eval_fun = eval_mod,
                         exported = list(),
                         #return_models = F,
                         #subset = NULL,
                         cores = ncores)


#saveRDS(training_model_all, file = "./data/models/training_model_all.RDS")



#Volume_dependent model lncs alone
vol_model_lncRNAs<- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                            arguments = args_3,
                            data = lncRNAS,
                            metadata = ct_metadata,
                            samplename = "seq_sample_id",
                            summary_fun = sum_fun,
                            eval_fun = eval_mod,
                            exported = list(),
                            #return_models = F,
                            #subset = NULL,
                            cores = ncores)

#saveRDS(vol_model_lncRNAs, file = "./data/models/vol_model_lncRNAs.RDS")





#Model for the trained versus untrained in lncRNAs

training_model_lncRNA<- seq_wrapper(fitting_fun = glmmTMB::glmmTMB,
                                 arguments = args_4,
                                 data = lncRNAS,
                                 metadata = ct_metadata,
                                 samplename = "seq_sample_id",
                                 summary_fun = sum_fun,
                                 eval_fun = eval_mod,
                                 exported = list(),
                                 #return_models = F,
                                 #subset = NULL,
                                 cores = ncores)


#saveRDS(training_model_lncRNA, file = "./data/models/training_model_lncRNA.RDS")
