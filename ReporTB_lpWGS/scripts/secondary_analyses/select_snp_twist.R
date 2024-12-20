library(tidyverse)

# Pass filter SNP
imp_maf5 <- read_tsv("result/plink/reportTB_filter_imp_maf5.bim",
                     col_names = FALSE) %>% pull(X2)

# Concordance
imp_maf5_concord <- read_tsv("result/snp_imp_maf5_concord.txt",
                             col_names=FALSE) %>% pull(X1)

# GWAS results
full <- read_csv("result/model_final/model_full.csv") %>% 
  filter(term=="value") %>% 
  select(snpID, pval)
full2 <- read_tsv("result/model_final/glmm_score_full.txt") %>%
  filter(! SNP %in% full$snpID) %>% 
  rename(snpID=SNP,pval=PVAL) %>% 
  select(snpID, pval)

full_all <- bind_rows(full, full2) %>% 
  #MAF > 5%
  filter(snpID %in% imp_maf5) %>% 
  #positions
  separate(snpID, into = c("CHR","POS","REF","ALT"), sep=":", remove=FALSE)

rm(full, full2)
gc()

full_sugg <- full_all %>% 
  filter(pval<5E-5) %>% 
  mutate(priority = case_when(
    pval<5E-8 ~"high",
    pval<5E-7 & snpID %in% imp_maf5_concord~"high-med",
    pval<5E-6 & snpID %in% imp_maf5_concord~"medium",
    pval<5E-5 & snpID %in% imp_maf5_concord~"med-low",
    TRUE~"low"),
    priority = factor(priority, 
                      levels=c("high","high-med",
                               "medium","med-low","low"))) %>% 
  arrange(priority,pval)

table(full_sugg$priority)

full_sugg %>% 
  filter(priority!="low") %>% 
  mutate(CHR = paste0("chr",CHR)) %>% 
  mutate(start=as.numeric(POS), end=as.numeric(POS)+1) %>% 
  select(CHR, start, end, snpID, priority) %>% 
  write_csv("result/other/snp_twist.csv")


#ZNF snps
full_sugg %>%  filter(priority=="high" & CHR=="3") %>% 
  arrange(as.numeric(POS))
