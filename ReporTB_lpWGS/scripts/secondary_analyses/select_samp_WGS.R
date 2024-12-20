# Select samples for 30X whole genome sequencing
# Focus on signif SNP

library(tidyverse)
library(SNPRelate)

#### DATA ####
#Signif SNP
signif <- read_csv("result/model_final/ReportTB_signif_snp_anno_allcis.csv") %>% 
  mutate(OR_group = ifelse(full_OR<1, "ref","alt")) %>% 
  select(snpID, full_group, cc_group, MAF_overall, full_OR, cc_OR, OR_group)

#SNP data
GDS_file1 <- "result/gds/reportTB_filter_maf1.gds"
genofile1 <- snpgdsOpen(GDS_file1)
#List SNP in this gds file
snp.id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))
snp.id1 <- snp.id1[snp.id1 %in% unique(signif$snpID)]

geno.num <- snpgdsGetGeno(genofile1, snp.id=snp.id1)

## Add row and column names
colnames(geno.num) <- snp.id1

sample.id1 <- read.gdsn(index.gdsn(genofile1, "sample.id"))
rownames(geno.num) <- sapply(strsplit(sample.id1,"_"), `[`, 1)
snpgdsClose(genofile1)

#Metadat
load("result/reportTB_meta_filter.RData")

#### SELECT SAMPLES #####
#Equalize case and control in subset
meta_FULL %>% count(tb)

set.seed(555)
meta.even <- meta_FULL %>% 
  group_by(tb) %>% 
  slice_sample(n=25)
meta.even %>% count(tb)

set.seed(666)
meta.ex <- meta_FULL %>% 
  filter(!sample_id %in% meta.even$sample_id) %>% 
  group_by(tb) %>% 
  slice_sample(n=5) 
meta.ex %>% count(tb)

#Total people with SNPs
snp.top <- as.data.frame(geno.num) %>% 
  rownames_to_column("ptID") %>% 
  pivot_longer(-ptID, names_to = "snpID") %>% 
  left_join(signif) %>% 
  filter(ptID %in% meta.even$sample_id) %>% 
  # mutate(value = ifelse(value==2,1,value)) %>% 
  group_by(snpID, MAF_overall, full_OR) %>% 
  # summarise(sumSNP = sum(value, na.rm=TRUE)) %>% 
  summarise(allele_0=sum(value==0),
            allele_1=sum(value==1),
            allele_2=sum(value==2)) %>% 
  ungroup() %>% 
  arrange(full_OR)
View(snp.top)

write_csv(snp.top, "result/other/ReportTB_snp_30X.csv")

samp1 <- meta.even %>% 
  ungroup() %>% 
  select(sample_id, tb) %>% 
  mutate(key = rep(1:25, 2)) %>% 
  pivot_wider(names_from = tb, values_from = sample_id) %>% 
  select(-key)

sampEx <- meta.ex %>% 
  ungroup() %>% 
  select(sample_id, tb) %>% 
  mutate(key = rep(1:5, 2)) %>% 
  pivot_wider(names_from = tb, values_from = sample_id, 
              names_prefix = "extra_") %>% 
  bind_rows(data.frame(matrix(nrow = 20))) %>% 
  select(-key, -`matrix.nrow...20.`)


bind_cols(samp1, sampEx) %>% 
  write_csv("result/other/ReportTB_samp_30X.csv", 
            na = "")
