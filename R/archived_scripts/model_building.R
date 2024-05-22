source("R/libraries.R")

#Load the functions most regularly used
source("R/Trainome_functions.R")


#Load the metadata
meta_df <- readRDS("data/contratrain_metadata.RDS")


#Load the lncRNA data
lncRNA_df <- readRDS("data/lncRNA_genes.RDS")


#Define the arguments. Will be looking at the interaction effect of condition
#and time on the expression of lncRNAs

args <- list(formula = y ~ efflibsize + time * condition + (1|participant),
                     
                     family = glmmTMB::nbinom2())


sum_fun_df <- function(x) {
  
  
  # Create a reference grid from the model
  em <- emmeans::emmeans(x, ~ time + condition, cov.reduce = FALSE)
  
  # Create a custom contrasts comparing conditions: set0 vs the average of set3 and set6 at t3 and t4
  # This should be double checked in additional analyses
  contr <- emmeans::contrast(em,  
                             list("t3, train - ctrl" = c( 0, -1, 0, 0, 0.5, 0, 0, 0.5, 0) - c( -1, 0, 0, 0.5, 0, 0, 0.5, 0, 0), 
                                  "t4, train - ctrl" = c( 0, 0, -1, 0, 0, 0.5, 0, 0, 0.5) - c( -1, 0, 0, 0.5, 0, 0, 0.5, 0, 0)),
                             adjust = "none")
  
  # Extract conditional effects
  results <-   tibble::as_tibble(coef(summary(x))$cond, 
                                 rownames = "coef") |>
    dplyr::select(coef, 
                  estimate = Estimate, 
                  se = 'Std. Error', 
                  z = 'z value', 
                  p = 'Pr(>|z|)') |>
    # Bind regression model to the custom contrasts
    rbind(data.frame(contr) |>
            dplyr::select(coef = contrast,
                          estimate, 
                          se = SE, 
                          z = z.ratio, 
                          p = p.value)) 
  
  ## Return results
  results
  
  
  
}



#Train the model
ct_models <- seqwrap(fitting_fun = glmmTMB::glmmTMB, 
                          arguments = args, 
                          data = lncRNA_df, 
                          metadata = meta_df, 
                          summary_fun = sum_fun_df, 
                          eval_fun = eval_mod,
                          #  subset = 1:5,
                          return_models = FALSE,
                          cores = ncores)


#emm1 <- emmeans(ct_models, specs = pairwise ~ time:condition)


err_warn_ct_gene <- ct_models$errors %>%
  dplyr::filter(map_lgl(warnings_fit, ~ !is.null(unlist(.x)))) %>%
  print()


keep_models_zi_gene <- ct_models$summaries[!(names(ct_models$summaries) %in% pull(err_warn_ct_gene, target))] 

gene_level_fig <- bind_rows(keep_models_zi_gene) %>%
  mutate(target = rep(names(keep_models_zi_gene), each = 12)) %>%
 dplyr:: filter(coef %in% c("t3, train - ctrl", "t4, train - ctrl")) #%>%
  
  ## Filter transcripts with p-values == 1 and abs estimates > 0
  # filter(!(p > 0.95 & abs(estimate) > 0.5)) %>%
  
   gene_level_fig %>% ggplot(aes(estimate, -log10(p))) + 
           
           geom_point(shape = 21)+ facet_wrap(~ coef) + 
           
           theme_gray()  + 
           labs(x = "Estimate (log-fold change)", 
                y = "-log10(P-value)") +
           
           labs(subtitle = "Gene level data") +
           
           theme_gray() + theme(legend.title = element_blank())
         
