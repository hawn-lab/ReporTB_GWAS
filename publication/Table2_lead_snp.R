library(tidyverse)

#GWAS results
anno <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  select_if(~sum(!is.na(.)) > 0)

bio_cols <- colnames(anno)[grepl("biotype",colnames(anno))]
anno_pc <- anno %>% 
  filter(if_any(all_of(bio_cols), ~ . == "protein_coding"))

#Concordance
concord <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz") %>% 
  filter(snpLow %in% anno_pc$snpID)

#Impute score
imp <- read_tsv("../ReporTB_lpWGS/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv.gz") %>% 
  filter(snpID %in% anno_pc$snpID)

#Combine
anno_format <- anno %>% 
  select(snpID, rsID, full_pval, full_OR, cc_pval, cc_OR,
         MAF_control, MAF_case,
         contains("_gene_"), contains("_biotype_"), contains("_distance_")) %>%
  mutate(across(contains("_distance_"), ~as.character(.))) %>% 
  pivot_longer(c(intragenic_gene_1:cis_distance_9)) %>% 
  drop_na(value) %>% 
  separate(name, into=c("cis","name","ID"), sep="_") %>% 
  pivot_wider() %>% 
  filter(biotype=="protein_coding") %>% 
  mutate(distance = ifelse(is.na(distance),0, distance)) %>% 
  select(cis, snpID, rsID, full_pval, full_OR, cc_pval, cc_OR,
         gene, distance, MAF_control, MAF_case) %>% 
  mutate(across(c(MAF_control, MAF_case), ~.*100)) %>% 
  mutate(across(c(full_pval, full_OR, cc_pval, cc_OR, MAF_control, MAF_case),
                ~signif(., digits=2))) %>% 
  mutate(distance = signif(as.numeric(distance)/1000, digits=3)) %>% 
  arrange(desc(cis), snpID, gene) %>% 
  left_join(concord %>% select(snpLow, pct) %>% 
              rename(snpID=snpLow, pct_concord=pct)) %>% 
  left_join(imp %>% select(snpID, INFO) %>% rename(imp_score=INFO)) %>% 
  select(snpID:MAF_case,imp_score,pct_concord, everything())

# 
load("../ReporTB_sceQTL/result/sceQTL/sceQTL_df.RData")
model_noCov %>% 
  filter(FDR < 0.05 & variable_group == "SNP") %>% 
  select(genotype, cell, gene, condition) %>% 
  arrange(genotype)
table(anno_format$snpID)

write_csv(anno_format, "Table2_lead_snp.csv", na="")
