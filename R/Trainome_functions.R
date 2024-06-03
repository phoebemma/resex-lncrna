
 source("R/archived_scripts/libraries.R")

 
 
 
 
 
 
 ncores <- parallel::detectCores()








#Function to extract the model evaluations into one dataframe



model_eval <-  function(x){ 
   bind_rows(x$evaluations) %>%
    mutate(target = names(x$evaluations))
  
}
  
 
#Function to extract model summary. Takes s input the model file and the number of coefficients
model_sum <- function(x, y){
  bind_rows(x$summaries) %>%
       mutate(target = rep(names(x$summaries), each = y))%>%
       subset(!coef == "(Intercept)") %>%
       mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
              log2fc = Estimate/log(2),
              
              fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"))
}




#Function to put the model summaries and evaluation in one dataframe
#Takes as input the model summary and model evaluation dataframes

filt_model_parameters <- function(x, y){
  inner_join(x, y, by = "target") %>%
    dplyr::filter(Pr...z.. <= 0.05 & fcthreshold == "s"  & pval.disp >= 0.05 
                  & pval.unif >= 0.05  )
  
}
 



#The plot_volcano function accepts as input the results of the model summary and evaluation
#And the string with which to name the volcano plot

#It returns a volcano plot

plot_volcano <- function(x, y){
  
  # add a column of NAs
  x$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  x$diffexpressed[x$log2fc > 0.6 & x$Pr...z.. < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  x$diffexpressed[x$log2fc < -0.6 & x$Pr...z..< 0.05] <- "DOWN"
  
  
  
  # Re-plot but this time color the points with "diffexpressed"
  p <- ggplot(data=x, aes(x=log2fc, y=-log10(Pr...z..), col=diffexpressed)) + geom_point() + theme_minimal()
  
  # Add lines as before...
  p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red")
  
  
  
  
  ## Change point color 
  
 
  p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
  
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  p3 <- p2 + scale_colour_manual(values = mycolors)
  
  
  
  
  # Now write down the name of genes beside the points...
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  x$delabel <- NA
  x$delabel[x$diffexpressed != "NO"] <- x$target[x$diffexpressed != "NO"]
  
  chart <- ggplot(data=x, aes(x=log2fc, y=-log10(Pr...z..), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    ggtitle(y) +
    theme(axis.text = element_text(size = 15),text = element_text(size = 15))+
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red")) +
    ylab("-log10(pvalue)")
  
  return(chart)
}











## Function for converting tibble to matrix and transform all values to integers
# Copied from Daniel's codes
make_matrix <- function(df,rownames = NULL){
  # Convert to matrix
  my_matrix <-  as.matrix(df)
  # Round all values
  my_matrix <- round(my_matrix, 0)
  
  # Change mode to integers
  mode(my_matrix) <- "integer"
  
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}


sum_fun <- function(x){
  
  cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$cond))),
                             coef(summary(x))$cond, 
                             
                             row.names = NULL)
  
  return(cond_effects)
  
}




sum_fun_lmer <- function(x){
  
  cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))))),
                                   coef(summary(x)), 
                                   
                                   row.names = NULL)
  
  return(cond_effects)
}




eval_mod <- function(x) {
  
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  
   disp <- DHARMa::testDispersion(sim, plot = FALSE) #tests if the simulated dispersion is equal to the observed dispersion
  unif <- DHARMa::testUniformity(sim, plot = FALSE) # tests if the overall distribution conforms to expectations
   zinfl <- DHARMa::testZeroInflation(sim, plot = FALSE) #tests if there are more zeros in the data than expected from the simulations
  
  results <- data.frame(pval.disp = disp$p.value, 
     pval.unif = unif$p.value, 
   pval.zinfl = zinfl$p.value)
  
  return(results)
}


#Evaluation model when using lmer
eval_mod_lmer <- function(x) {
  
  sim <- DHARMa::simulateResiduals(x, n = 1000)
  
#  disp <- DHARMa::testDispersion(sim, plot = FALSE) #tests if the simulated dispersion is equal to the observed dispersion
  unif <- DHARMa::testUniformity(sim, plot = FALSE) # tests if the overall distribution conforms to expectations
#  zinfl <- DHARMa::testZeroInflation(sim, plot = FALSE) #tests if there are more zeros in the data than expected from the simulations
  
  results <- data.frame(#pval.disp = disp$p.value, 
                        pval.unif = unif$p.value) #, 
                        #pval.zinfl = zinfl$p.value)
  
  return(results)
}





#function to read Kallisto files
read_kallisto_output <- function(file){
  df <- readr::read_table(file)
  df <- separate(df, target_id, c("transcript_ID", "gene_ID", "Havana_gene_ID",
                                  "Havana_transcript_ID", "transcript_name",
                                  "gene_name", "sequence_length", 
                                  "transcript_biotype"), sep = "\\|")
  df <- df %>% select(transcript_ID, gene_ID, transcript_name, gene_name,
                      transcript_biotype, length, est_counts, tpm)
  #extract all transcripts with est_count below 1
  # df <- filter(df, est_counts > 1)
  return(df)
}





## a function to extract all the lncRNAs based on EBI's definition in the link
#https://www.ensembl.org/info/genome/genebuild/biotypes.html
read_biotype_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- df %>% filter(transcript_biotype %in% c("processed_transcript", "lncRNA",
                                                "lincRNA", "3prime_overlapping_ncrna",
                                                "antisense", "non_coding", "sense_intronic",
                                                "sense_overlapping", "TEC", "known_ncrna",
                                                "bidirectional_promoter_lncrna",
                                                "macro_lncRNA"))
  return(df)
}
#function to extract those that are lncRNAs among the lncRNA biotype
read_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- filter(df, transcript_biotype == "lncRNA")
  return(df)
}

#A function to load load non-coding RNAs and protein coding RNAs
read_coding_and_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- df %>% filter(transcript_biotype %in% c("processed_transcript","protein_coding", "lncRNA",
                                                "lincRNA", "3prime_overlapping_ncrna",
                                                "antisense", "non_coding", "sense_intronic",
                                                "sense_overlapping", "TEC", "known_ncrna",
                                                "bidirectional_promoter_lncrna",
                                                "macro_lncRNA"))
  return(df)
}
## Using data table with large number of rows instead...
# function to remove all columns except TPM and transcript name
# the function also combines all files into a data fram

####The subfunction which extracts the desired biotype is interchangable

####extract all transcripts
extract_all_transcripts <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  transcripts <- list()
  
  for(i in 1:length(files)){
    
    transcripts[[i]] <- read_kallisto_output(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
  }
  
  comd.df <- data.table::rbindlist(transcripts)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
}



extract_coding_lncRNA <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  imports <- list()
  
  for(i in 1:length(files)){
    
    imports[[i]] <- read_coding_and_lncRNAs(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(imports)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
  
  
}



extract_lncRNAs <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  imports1 <- list()
  
  for(i in 1:length(files)){
    
    imports1[[i]] <- read_biotype_lncRNAs(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub("^[^-]*-", "", gsub(".tsv", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_ID, file_id, tpm)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(imports1)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID ~ file_id, value.var = "tpm")
  
  return(data.frame(comb.df))
  
}






#function to read Rsem .gene.results files
read_Rsem_genes <- function(file){
  df <- readr::read_delim(file)
  df <- df %>% dplyr::select(gene_id, length, TPM, effective_length, expected_count)
  return(df)
}

#function to read Rsem isoform.results file
read_Rsem_isoforms <- function(file){
  df <- readr::read_delim(file)
  df <- df %>% dplyr::select(transcript_id, effective_length, length, expected_count,)
  return(df)
}



extract_rsem_isoform_counts <- function(folder){
  
  files <- list.files(folder, pattern =".isoforms.results") 
  
  isoforms <- list()
  
  for(i in 1:length(files)){
    
    isoforms[[i]] <- read_Rsem_isoforms(paste0(folder,"/", files[i])) %>% 
      mutate(file_id =gsub("^[^_]*-",  "",gsub(".isoforms.results", "", files[i]))) %>% ## This removes the file number
      dplyr::select(transcript_id, file_id, expected_count, effective_length, length )
  }
  
  comd.df <- data.table::rbindlist(isoforms)
  
  comb.df <- data.table::dcast(comd.df, transcript_id ~file_id, value.var = "expected_count")
  
  return(data.frame(comb.df))
  
}




#Extracts the TPM 

extract_rsem_gene_counts <- function(folder){
  
  files <- list.files(folder, pattern =".genes.results") 
  
  genes <- list()
  
  for(i in 1:length(files)){
    
    genes[[i]] <- read_Rsem_genes(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub(".genes.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(gene_id, file_id, length, effective_length, expected_count, TPM)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(genes)

  
  comb.df <- data.table::dcast(comd.df, gene_id ~ file_id, value.var = "TPM")
  
  return(data.frame(comb.df))
}




#extract the length of the genes. This extracts the column called "effective length"
extract_rsem_gene_lengths <- function(folder){
  
  files <- list.files(folder, pattern =".genes.results") 
  
  genes <- list()
  
  for(i in 1:length(files)){
    
    genes[[i]] <- read_Rsem_genes(paste0(folder,"/", files[i])) %>% 
      mutate(file_id = gsub(".genes.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(gene_id, file_id, length, effective_length, expected_count)
    
    
    
  }
  
  comd.df <- data.table::rbindlist(genes)
  
  
  comb.df <- data.table::dcast(comd.df, gene_id ~ file_id, value.var = "effective_length")
  
  return(data.frame(comb.df))
}




extract_rsem_isoform_lengths <- function(folder){
  
  files <- list.files(folder, pattern =".isoforms.results") 
  
  isoforms <- list()
  
  for(i in 1:length(files)){
    
    isoforms[[i]] <- read_Rsem_isoforms(paste0(folder,"/", files[i])) %>% 
      mutate(file_id =gsub(".isoforms.results", "", files[i])) %>% ## This removes the file number
      dplyr::select(transcript_id, file_id, expected_count, effective_length, length )
  }
  
  comd.df <- data.table::rbindlist(isoforms)
  
  comb.df <- data.table::dcast(comd.df, transcript_id ~file_id, value.var = "effective_length")
  
  return(data.frame(comb.df))
  
}


#function to read Splice_q results file
read_Splice_Q <- function(file){
  df <- readr::read_tsv(file)
  df <- df %>% select(transcript_ID, intron_ID, score)
  return(df)
}

extract_splice_q <- function(folder){
  
  files <- list.files(folder, pattern =".tsv") 
  
  splice_list <- list()
  
  for(i in 1:length(files)){
    
    splice_list[[i]] <- read_Splice_Q(paste0(folder, files[i])) %>% 
      mutate(file_id = gsub(".tsv", "", files[i])) %>% ## This removes the file number
      dplyr::select(transcript_ID, intron_ID, file_id, score)
  }
  
  comd.df <- data.table::rbindlist(splice_list)
  
  comb.df <- data.table::dcast(comd.df, transcript_ID + intron_ID   ~ file_id, value.var = "score")
  
  return(data.frame(comb.df))
}



#Function to round all floating numbers to two decimal places

round_df <- function(x, digits) {
 
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


