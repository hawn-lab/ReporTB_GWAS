library(tidyverse)
signif <- 5E-8

#### covariate as outcome ####
cov <- read_csv("result/model_final/model_signif_cov.csv") %>% 
  filter(term=="value") %>% 
  select(model, snpID, pval)

cov_signif <- cov %>% 
  filter(pval<signif) %>% 
  group_by(snpID) %>% 
  summarise(note = paste("signif predictor of",
                         paste0(unique(model), collapse = ", "))) %>% 
  ungroup()

write_csv(cov_signif, file = "result/model_final/ReporTB_signif_snp_notes.csv")

#### covariates in full model ####
signif_full <- read_csv("result/model_final/model_full.csv") %>% 
  filter(pval<signif & term=="value") %>% 
  pull(snpID)

cov_signif2 <- read_csv("result/model_final/model_full.csv") %>% 
  filter(snpID %in% signif_full) %>% 
  filter(term %in% c("sexF","bl_age",
                     "smokhx2Y","bl_hivY","diabetesY")) %>% 
  select(snpID, term, pval) %>% 
  mutate(pval=case_when(pval<signif~term)) %>% 
  pivot_wider(names_from = term, values_from = pval) %>% 
  mutate(note2 = paste(sexF, bl_age, smokhx2Y, bl_hivY, diabetesY, sep=","),
         note2 = gsub(",NA|NA$","",note2),
         note2 = paste(note2, "signif in model")) %>% 
  select(snpID, note2)

write_csv(cov_signif2, file = "result/model_final/ReporTB_signif_cov_notes.csv")
