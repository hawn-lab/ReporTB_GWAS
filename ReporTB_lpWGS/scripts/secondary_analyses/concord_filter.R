library(tidyverse)

# 1 and 5% minor allele frequency
maf5 <- read_tsv("result/plink/reportTB_filter_imp_maf5.bim", col_names = FALSE) %>% 
  pull(X2) 

maf1 <- read_tsv("result/plink/reportTB_filter_imp_maf1.bim", col_names = FALSE) %>% 
  pull(X2)

#Concordant SNPs
concord <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz")

# Save lists
concord_pass <- concord %>% 
  filter(pct>=60)

concord_pass %>% 
  filter(snpLow %in% maf5) %>% 
  select(snpLow) %>% 
  write_tsv("result/snp_imp_maf0_concord.txt", col_names = FALSE)
concord_pass %>% 
  filter(snpLow %in% maf5) %>% 
  select(snpLow) %>% 
  write_tsv("result/snp_imp_maf5_concord.txt", col_names = FALSE)
concord_pass %>% 
  filter(snpLow %in% maf1) %>% 
  select(snpLow) %>% 
  write_tsv("result/snp_imp_maf1_concord.txt", col_names = FALSE)
