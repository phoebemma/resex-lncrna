log_simp <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/log_simple_stringent_training_model.RDS")
log_simp_ev <- bind_rows(within(log_simp$evaluations, rm(NUDT12, PPIL6, RAG2)))%>%
  mutate(target = names(within(log_simp$evaluations, rm(NUDT12, PPIL6, RAG2))))

hist(log_simp_ev$pval.unif, main = "distribution of p unif values log converted dependent y")


log_simp_summ <- bind_rows(within(log_simp$summaries, rm(NUDT12, PPIL6, RAG2))) %>%
  subset(!coef == "(Intercept)") %>%
   mutate(target = rep(names(within(log_simp$summaries, rm(NUDT12, PPIL6, RAG2))), each = 19)) #%>%
  # mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
  #        log2fc = Estimate/log(2),
  #        
  #        fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  # filter(fcthreshold == "s" & adj.p <= 0.05 )

hist(log_simp_summ$Pr...t.., main = "distribution of p values log converted model")
length(unique(log_simp_summ$coef))



simp_string <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/simpler_training_correlation_model.RDS")
simp_string_ev <- model_eval(simp_string)

hist(simp_string_ev$pval.unif, main = "punif distribution in simple model")


simp_string_sum <- bind_rows(simp_string$summaries) %>%
  subset(!coef == "(Intercept)") %>%
  mutate(target = rep(names(simp_string$summaries), each = 19)) #%>%
  # mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
  #        log2fc = Estimate/log(2),
  #        
  #        fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  # filter(fcthreshold == "s" & adj.p <= 0.05 )

hist(simp_string_sum$Pr...t.., main = "distribution p values in simple model")

length(unique(simp_string_sum$coef))



int_model <- readRDS("data/models/seqwrap_generated_models/trained_untrained_midexercise/interaction_training_correlation_model.RDS")
int_mod_ev <- bind_rows(within(int_model$evaluations, rm(GABRG1)))%>%
  mutate(target = names(within(int_model$evaluations, rm(GABRG1))))
  
  hist(int_mod_ev$pval.unif, main = "p unif distribution in interaction model")
  
  
int_mod_summ <- bind_rows(within(int_model$summaries, rm(GABRG1))) %>%
  subset(!coef == "(Intercept)") %>%
  mutate(target = rep(names(within(int_model$summaries, rm(GABRG1))), each = 89))%>%
  mutate(adj.p = p.adjust(Pr...t.., method = "fdr"),
         log2fc = Estimate/log(2),
         
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>%
  filter(fcthreshold == "s" & adj.p <= 0.05 )

hist(int_mod_summ$Pr...t.., main = "distribution p values in interaction model")

length(unique(int_mod_summ$coef))
