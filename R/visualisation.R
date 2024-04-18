#CEMiTool analyses
#| #Load the needed libraries
source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")



#Load the lncRNA dataset

# lncs <- readRDS("data/lncRNA_genes.RDS")
# 
# lncs_ <- as.matrix(lncs[, -1])
# rownames(lncs_) <- lncs$gene_name
# 
# lncs_ <- as.data.frame(lncs_)
# 
# cem <- cemitool(lncs_)
# 
# 
# 
# #load metadata
# met <- readRDS("data/contratrain_metadata.RDS")
# 
# annot <- met %>%
#   dplyr::select(seq_sample_id, condition)
# 
# 
# # r ename the column names as stated on https://cemitool.sysbio.tools/tutorial
# colnames(annot)[colnames(annot) == "seq_sample_id"] <- "SampleName"
# 
# 
# colnames(annot)[colnames(annot) == "condition"] <- "Class"
# 
# 
# 
# ## Read example pathways file
# gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
# gmt_in <- read_gmt(gmt_fname)
# ## Get example interactions file
# int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))
# ## Run CEMiTool
# 
# 
# cem <- cemitool(lncs_, annot = annot, gmt = gmt_in, interactions = int_df)
# 
# 
# 
# 
# head(lncs)
# 
# cem <- cemitool(lncs)



#Load the volume model

Vol_all <- readRDS("data/models/Filtered_coefs/Vol_model_all.RDS")
#check the unique coefficients
unique(Vol_all$coef)

# Vol_df <- Vol_all %>%
#   dplyr::select(target, coef)

#select those at post exercise

t4_vol_all <- Vol_all %>%
  dplyr::filter(coef == "timet4")%>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)

#select those at mid exercise

t3_vol <-  Vol_all %>%
  dplyr::filter(coef == "timet3" ) %>%
  dplyr::filter(log2fc >= 1 | log2fc <= -1)



#Load the metadata
meta_df <- readRDS("data/contratrain_metadata.RDS") %>%
  dplyr::select(seq_sample_id, time, condition, group)


#Load lncs data

lncs_df <- readRDS("data/lncRNA_genes.RDS")


#Extrat the DE lncs at mid exercise and post exercise

t4_lncs <- lncs_df[lncs_df$gene_name %in% t4_vol_all$target,]

t3_lncs <- lncs_df[lncs_df$gene_name %in% t3_vol$target,]


t4_lncs <- pivot_longer(data = t4_lncs,
                        cols = -(gene_name),
                        names_to = "seq_sample_id",
                        values_to = "counts")






#Merge the t4 lncs to the metadata

t4_df <- t4_lncs %>%
  inner_join(meta_df, by = "seq_sample_id")

t4_df$log_counts <- log(t4_df$counts)

df <- ggplot(data = t4_df, 
             mapping = aes(x = time,
                           y = gene_name,
                           fill = log_counts)) +
  geom_tile() +
  facet_grid(~ condition)
df




#Merge the t4 lncs to the metadata





#Repeat for mid exercise


t3_lncs <- pivot_longer(data = t3_lncs,
                        cols = -(gene_name),
                        names_to = "seq_sample_id",
                        values_to = "counts")
t3_df <- t3_lncs %>%
  inner_join(meta_df, by = "seq_sample_id")

t3_df$log_counts <- log(t3_df$counts)

df <- ggplot(data = t3_df, 
             mapping = aes(x = time,
                           y = gene_name,
                           fill = log_counts)) +
  geom_tile() +
  facet_grid(~ condition)
df