library(tidyverse)
library(data.table)

#Sample metadata
attach("~/project/data/metadata/ReportTB_metadata_clean.RData")
#Samples with genotype data
samps <- fread(file="~/project/result/plink/reportTB.fam",
               col.names=c("FID","IID","father","mother","sex","pheno")) %>% 
  separate(IID, into=c("sample_id"), sep="_", remove = FALSE, extra="drop")

#Samples missing metadata
pt_missing <- samps %>% 
  filter(!sample_id %in% meta_FULL$sample_id) %>% 
  pull(IID)

#Save
write_tsv(as.data.frame(pt_missing), col_names = FALSE, 
          file = "~/project/result/other/pt_remove.txt")

######
#Duplicates and unconfirmed twins
# See twin_check.R and metadata/RePORTBR_potential_twins Notes MCF Oct27.xlsx 
# for final classifications of potential twins / duplicates

dups <- c("7212-BA-2547",
          "7212-BA-0678",
          "7212-BA-0839",
          "7212-BA-2320",
          "7212-BA-1047",
          "7212-BA-0534",
          "7212-BA-0875",
          "7212-BA-1572",
          "7212-BA-2488")

dup_remove <- samps %>% 
  filter(sample_id %in% dups) %>% 
  pull(IID)

#Save
write_tsv(as.data.frame(dup_remove), col_names = FALSE, 
          file = "~/project/result/other/dups_remove.txt")

#### Case-control subset ####
#participants to keep 
cc_keep <- samps %>% 
  filter(sample_id %in% meta_CC$sample_id &
           !sample_id %in% dups) %>% 
  pull(IID)

write_tsv(as.data.frame(cc_keep), col_names = FALSE, 
          file = "~/project/result/other/cc_keep.txt")

print("Done")
