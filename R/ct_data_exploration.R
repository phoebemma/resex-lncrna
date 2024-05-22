#Load the needed libraries


#Load the functions most regularly used
source("R/Trainome_functions.R")


#Load the meatadata 
meta_df <- readRDS("data/contratrain_metadata.RDS")


#Load the file containing protein-coding genes and lncs


genes <- readRDS("data/protein_coding_and_lncRNA_genes_combined.RDS")

# mRNA_genes <-readRDS("data/protein_coding_genes.RDS")
# 


genes%>%
  group_by(gene_name, gene_biotype) %>%
  ggplot(aes(x = gene_biotype))+
  geom_bar()+
  theme(axis.text.x = element_text( vjust = 1, hjust = 0.5))+
  geom_text(stat = "count",  aes(label= after_stat(count)), vjust = -0.5)

# 
# trained_df <- meta_df %>%
#   dplyr::filter(training_status == "trained")
# 
# 
# trained_genes <- genes %>%
#   dplyr::select(gene_name, gene_biotype, any_of(trained_df$seq_sample_id)) 
# 
# trained_genes %>%
#   group_by(gene_name, gene_biotype) %>%
#   ggplot(aes(x = gene_biotype))+
#   geom_bar()+
#   theme(axis.text.x = element_text( vjust = 1, hjust = 0.5))+
#   geom_text(stat = "count",  aes(label= after_stat(count)), vjust = -0.5)
# 



# 
# #convert the genes into  long formats

long_genes <- genes %>%
  pivot_longer(names_to = "sample_id",
               values_to = "count",
               cols = 2 : 137)


#merge the metadata with the genes long data

merged_df <- inner_join(meta_df, long_genes, by = c("seq_sample_id" = "sample_id"))


merged_df %>%
  group_by(seq_sample_id, gene_biotype, condition) %>%
  ggplot(aes(x = gene_biotype)) +
  geom_bar()+
  facet_grid(training_status~time)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5))+
  geom_text(stat = "count",  aes(label=..count..), vjust = -0.5)
