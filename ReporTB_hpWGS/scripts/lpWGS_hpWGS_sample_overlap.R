library(tidyverse)

hp <- read_tsv("/home/ec2-user/SEAsnake/result/7_plink_filter/reportb_wgs_all_maf0.fam",
               col_names=FALSE) %>% 
  select(X1) %>% 
  separate(X1, into=c("trash","sample_id"), sep=" ", remove=FALSE) %>% 
  select(-trash) %>% 
  mutate(sample_id = gsub("/home/ec2-user/SEAsnake/result/3_bam_filter/",
                          "", sample_id),
         sample_id = gsub("_S1_Aligned.sortedByCoord.filter.bam",
                          "", sample_id)) %>% 
  rename(hp=X1)

lp <- read_tsv("/home/ec2-user/SEAsnake/ReporTB_lpWGS/plink/reportTB_filter.fam",
               col_names=FALSE) %>% 
  select(X1,X2) %>% 
  separate(X2, into=c("sample_id"), sep="_", remove=FALSE) %>% 
  rename(lp=X2)

inner_join(hp, lp) %>% 
  select(X1,lp) %>% 
  write_tsv("/home/ec2-user/SEAsnake/result2/lpWGS/hp_samples.txt", col_names = FALSE)
