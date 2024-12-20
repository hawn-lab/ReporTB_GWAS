library(tidyverse)

# SNP modeled (filters = MAF 5%, HWE, imp > 0.5)
impPF <- read_tsv("../ReporTB_lpWGS/result/snp_impute_pf05.txt", 
                  col_names=FALSE) %>% pull(X1) %>% sort()
maf5 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_maf5.bim",
                 col_names=FALSE) %>% pull(X2) %>% sort()
snp_maf5 <- intersect(impPF, maf5)
length(snp_maf5)

# SNP modeled and concord (filters = MAF 5%, HWE, imp > 0.5, concord>60)
maf5_concord <- read_csv("../ReporTB_lpWGS/result/snp_maf5_concord.txt", 
                         col_names=FALSE) %>% pull(X1) %>% sort()
snp_final <- intersect(impPF, maf5_concord)
length(snp_final)

# SNP 
