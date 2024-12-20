library(tidyverse)

#### MAGMA data ####
maf5_concord <- read_tsv("result/snp_imp_maf5_concord.txt", col_names=FALSE) %>% 
  pull(X1)

fdr <- read_csv("result/model_final/model_full.csv") %>% 
  filter(snpID %in% maf5_concord) %>% 
  filter(term == "value") %>% 
  rename(SNP=snpID, P=pval) %>% 
  select(SNP, P)

scores <- read_tsv("result/model_final/glmm_score_full.txt") %>% 
  filter(SNP %in% maf5_concord) %>% 
  filter(!SNP %in% fdr$SNP) %>% 
  rename(P=PVAL) %>% 
  select(SNP, P)

#Combine and save
bind_rows(fdr,scores) %>% 
  mutate(NOBS=2754) %>% 
  select(SNP, P, NOBS) %>% 
  write.table(file="result/magma/ReportTB_magma_pval.txt", 
              sep = " ", row.names = FALSE, col.names = TRUE, 
              quote = FALSE)

#### MSigDB anno databases for magma ####
library(msigdbr)

db_h <- msigdbr("human", "H") %>% 
  select(gs_name, human_entrez_gene) 

db_c2 <- msigdbr("human", "C2") %>% 
  dplyr::filter(gs_subcat %in% c("CP:BIOCARTA","CP:KEGG","CP:REACTOME")) %>% 
  dplyr::select(gs_name, human_entrez_gene) 

db_c5 <- msigdbr("human", "C5") %>% 
  dplyr::filter(gs_subcat == "GO:BP") %>% 
  dplyr::select(gs_name, human_entrez_gene) 

bind_rows(db_h, db_c2, db_c5)%>% 
  write_delim(file = "/Users/kadm/Applications/magma_v1.10_mac/ref/MSigDB/msigdb_7.5.1_H.C2.C5BP.txt",
              col_names = FALSE, quote = "none")
