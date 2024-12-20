library(tidyverse)
load("../metadata/ReportTB_metadata_clean.RData")

key <- read_tsv("result/plink/reportTB_filter_imp_maf5.fam", col_names = FALSE) %>% 
  select(X1,X2) %>% 
  separate(X2, into=c("sample_id"), sep="_", remove=FALSE, extra = "drop")

#ALL FILES first two columns as family id and individual id followed by variable columns. Pull from key
#### TB outcome ####
#TB.phen case-control class
meta_FULL %>% 
  select(tb, sample_id) %>% 
  mutate(tb=as.numeric(as.factor(tb))) %>% 
  inner_join(key, by="sample_id") %>% 
  select(X1,X2,tb) %>% 
  write_tsv("result/heritability/TB.phen", col_names = FALSE)

#### Categorical ####
meta_FULL %>% 
  mutate(smokhx2 = fct_recode(smokhx, 
                              "curr_past"="current",
                              "curr_past"="past")) %>% 
  select(sample_id, sex, smokhx2, bl_hiv, diabetes) %>% 
  inner_join(key, by="sample_id") %>% 
  select(X1,X2,sex:diabetes) %>% 
  write_tsv("result/heritability/covar_cat.txt", col_names = FALSE)

# No TB risk factors
meta_FULL %>% 
  mutate(smokhx2 = fct_recode(smokhx, 
                              "curr_past"="current",
                              "curr_past"="past")) %>% 
  select(sample_id, sex) %>% 
  inner_join(key, by="sample_id") %>% 
  select(X1,X2,sex) %>% 
  write_tsv("result/heritability/covar_cat_noTB.txt", col_names = FALSE)

#### Continuous ####
attach("result/kinship/reportTB_PCair.RData")
PC.dat <- as.data.frame(PC.dat) %>% 
  rownames_to_column("X2") %>% 
  select(X2,PC1,PC2)

meta_FULL %>% 
  mutate(smokhx2 = fct_recode(smokhx, 
                              "curr_past"="current",
                              "curr_past"="past")) %>% 
  select(sample_id, bl_age) %>% 
  inner_join(key, by="sample_id") %>% 
  inner_join(PC.dat, by="X2") %>% 
  select(X1,X2,PC1,PC2,bl_age) %>% 
  write_tsv("result/heritability/covar_cont.txt", col_names = FALSE)

# PCs only
meta_FULL %>% 
  mutate(smokhx2 = fct_recode(smokhx, 
                              "curr_past"="current",
                              "curr_past"="past")) %>% 
  select(sample_id, bl_age) %>% 
  inner_join(key, by="sample_id") %>% 
  inner_join(PC.dat, by="X2") %>% 
  select(X1,X2,PC1,PC2) %>% 
  write_tsv("result/heritability/covar_PC.txt", col_names = FALSE)

#### TB prevalence in cohort ####
meta_FULL %>% 
  inner_join(key, by="sample_id") %>% 
  count(tb) %>% 
  mutate(pct=n/sum(n)) %>% 
  filter(tb=="culture_confirm") %>% 
  pull(pct)
