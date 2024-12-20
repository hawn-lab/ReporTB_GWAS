library(tidyverse)
library(data.table)

#Sample metadata
attach("~/project/data/metadata/ReportTB_metadata_clean.RData")

#Samples with genotype data
fam <- fread(file="~/project/result/plink/reportTB.fam",
             col.names = c("FID","IID","father","mother","sex","pheno")) %>% 
  separate(IID, into=c("sample_id"), sep="_", remove = FALSE, extra = "drop") %>% 
  #Add case control classification
  left_join(meta_FULL %>% select(sample_id, tb), by="sample_id") %>% 
  mutate(tb = case_when(tb=="culture_confirm"~2,
                        tb=="contact"~1,
                        TRUE~-9)) %>% 
  select(-pheno, -sample_id)

#Save
write_tsv(fam, "~/project/result/plink/reportTB.fam", 
          quote="none", col_names=FALSE)

#Filter file for case-control subset
# fam %>% 
#   select(IID) %>% 
#   separate(IID, into=c("sample_id"), sep="_", remove = FALSE, extra = "drop") %>% 
#   inner_join(meta_CC %>% select(sample_id, tb)) %>% 
#   select(IID) %>% 
#   write_tsv("~/project/result/other/cc_keep.txt", col_names = FALSE)
  
