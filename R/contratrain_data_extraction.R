#Load the needed libraries
#source("R/libraries.R")

library(dplyr)
library(trainomeMetaData)
library(trainomeHelper)



#Load the functions most regularly used
source("R/Trainome_functions.R")


# Download metadata from  the TrainomeMetadata package
data("ct_participants")
data("ct_samples")  
data("ct_strength")
data("ct_thickness")


#rename the values in the time column to match those in the other dataframes
ct_strength["time"][ct_strength["time"] == "pre"] <- "t1"
ct_strength["time"][ct_strength["time"] == "mid"] <- "t3"
ct_strength["time"][ct_strength["time"] == "post"] <- "t4"

ct_metadata <- ct_samples %>% inner_join(ct_participants, by = c("study", "participant", "sex")) %>%
  dplyr::select(study, participant, sex, condition,leg, time, seq_sample_id) %>%
  inner_join(ct_thickness, by = c("study", "participant", "leg", "condition", "time") ) %>%
  #drop the rows containing missing data
  drop_na()%>%
  #drop columns not needed for lncRNA analyses
  dplyr::select(-c("study", "leg", "thickness")) %>%
  #inner_join(ct_strength, by = c("study", "participant", "group", "leg", "condition", "time"))
  mutate(time = factor(time, levels = c("t1", "t3",   "t4")), 
         condition = factor(condition, levels = c("set0", "set3", "set6")), 
         sex = factor(sex, levels = c("male", "female")),
         training_status = if_else(condition == "set0", "untrained", "trained"),
         training_status = factor(training_status,
                                  levels = c("untrained", "trained"))) 



#Visualize metadata

#save the resulting image in jpeg format
jpeg(filename = "plots/lncRNA_metadata_distribution.jpeg",
     width = 850, height = 500, quality = 100)
ct_metadata %>%
  ggplot(aes(condition, color = condition, fill = condition))+
  geom_bar( alpha = 0.5) +
  facet_grid(time ~ group) +
  ggtitle("Distribution of RNA sequence data in the Contratrain study (Group by time, by condition)")+
  theme(axis.text = element_text(size = 10), text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
dev.off()


#Load the RSEM TPM counts

TPM <- extract_rsem_gene_counts("data/Contratrain_rsem_genes/")%>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  dplyr::select(-(gene_id))
 
TPM[,-1] <-  round(TPM[,-1], 0)

#saveRDS(TPM, file = "data/Ct_genes_TPM.RDS")

# download the contratrain gene counts
download_ome(download = "ct_gene_rsem")


ct_genes <- read_csv("ome-data/ct_gene_rsem.csv") %>%
  separate(gene_id, c("gene_id", "gene_name"), sep = "_", extra = "merge") %>%
  #drop gene_id, select gene_name and any of the sample names that match sample name in metadata
  dplyr::select(gene_name, any_of(ct_metadata$seq_sample_id))


# # Keep nonzero rows
nonzero <- ct_genes %>%
  dplyr::filter(rowSums(ct_genes[,-1]) != 0)
 
 
# # Keep filtered genes (based on group)
filtered <- nonzero %>%
  dplyr::filter(filterByExpr(nonzero[,-1], group = ct_metadata$time)) 


#saveRDS(filtered, file = "./data/filtered_all_gene_counts.RDS")




# Create dge lists and calculate norm factors
dge_ct   <- DGEList(filtered[,-1])


dge_ct  <- calcNormFactors(dge_ct  )



## Add effective library size to 
ct_metadata <- dge_ct$samples %>%
  rownames_to_column(var = "seq_sample_id") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  
  mutate(efflibsize = (lib.size * norm.factors) / median(lib.size * norm.factors)) %>%
  #Remove the group.x column as it was used for filtering
  dplyr::select(-("group.x")) %>%
  #rename the group.y column as it it the original grouping of the data
  #rename the library size and the normalization factors to reflect that all biotypes are accounted for
  rename("group" = "group.y", 
         "lib.size_all" = "lib.size",
         "norm.factors_all" = "norm.factors") %>%
  
  print()

#Save metadata for further use

#saveRDS(ct_metadata, file = "data/contratrain_metadata.RDS")







#Extract only the lncRNAs from the data



#Use the ensemble database to get the annotation of the transcripts
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)
attributes <- listAttributes(ensembl)

# #extract GENE biotypes and  names
annotation<- getBM(attributes = c("gene_biotype", "external_gene_name"), values =filtered$gene_name, mart = ensembl )
annotation <- inner_join(filtered, annotation, by= c("gene_name" = "external_gene_name"))%>%
  #Round the counts to zero decimal place
  mutate_at(2:137, ~ as.integer(., 0))


#check the unique biotypes
unique(annotation$gene_biotype)

#Extract only the those annotated as protein coding and lncRNA"
#This will be used for exploratory data analyses
mRNAs_and_lncs <- annotation %>%
  dplyr::filter(gene_biotype == "lncRNA" | gene_biotype == "protein_coding")

unique(mRNAs_and_lncs$gene_biotype)

#saveRDS(mRNAs_and_lncs, file = "data/protein_coding_and_lncRNA_genes_combined.RDS")


#xtract only the lncRNAs and save in the data file


lncRNAS <- annotation %>%
  dplyr::filter(gene_biotype == "lncRNA") %>%
  dplyr::select(-gene_biotype)


#save the lncRNA count file
#saveRDS(lncRNAS, file = "data/lncRNA_genes.RDS")



#Extract only the protein-coding genes
#Needed for correlation testing

mRNAs <- annotation %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::select(-gene_biotype)


#Save the file
#saveRDS(mRNAs, file = "data/protein_coding_genes.RDS")












#Extracting the library size based only on the lncRNAs

# # Keep nonzero rows
nonzero_2 <- lncRNAS %>%
  dplyr::filter(rowSums(lncRNAS[,-1]) != 0)


# # Keep filtered lncs (based on group)
filtered_2 <- nonzero_2 %>%
  dplyr::filter(filterByExpr(nonzero_2[,-1], group = ct_metadata$time)) 

#save the filtered lncRNA count file
#saveRDS(filtered_2, file = "data/lncRNA_genes.RDS")




# Create dge lists and calculate norm factors
dge_ct_2   <- DGEList(filtered_2[,-1])


dge_ct_2  <- calcNormFactors(dge_ct_2 )



## Add effective library size to 
ct_metadata <- dge_ct_2$samples %>%
  rownames_to_column(var = "seq_sample_id") %>%
  inner_join(ct_metadata, by = "seq_sample_id") %>%
  
  mutate(efflibsize_lncs = (lib.size * norm.factors) / median(lib.size * norm.factors)) %>%
  #Remove the group.x column as it was used for filtering
  dplyr::select(-("group.x")) %>%
  #rename the group.y column as it it the original grouping of the data
  rename("group" = "group.y") %>%
  
  print()

#saveRDS(ct_metadata, file = "data/contratrain_metadata.RDS")


