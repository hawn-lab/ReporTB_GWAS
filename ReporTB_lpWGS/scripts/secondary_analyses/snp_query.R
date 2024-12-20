library(tidyverse)
dt <- Sys.Date() %>% gsub("-",".",.)

#### SNP lists ####
# called in lpWGS
maf0 <- read_tsv("result/plink/reportTB_filter.bim",
                 col_names=FALSE) %>% pull(X2)
# imputation score > 0.6
maf0_imp <- read_tsv("result/plink/reportTB_filter_imp.bim",
                     col_names=FALSE) %>% pull(X2)
# imputation score > 0.6, MAF > 5%, HWE < 1E-6
maf5_imp <- read_tsv("result/plink/reportTB_filter_imp_maf5.bim",
                     col_names=FALSE) %>% pull(X2)
#Imputation score
imp_score <- read_tsv("result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv.gz") %>% 
  select(snpID, INFO) %>% 
  rename(imp_score=INFO)

#concordance
concord <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz") %>%
  filter(snpLow %in% maf0) %>% 
  select(snpLow,tot_calls, pct) %>% 
  rename(snpID = snpLow)

# combined 
dat <- concord %>% 
  left_join(imp_score, by = join_by(snpID)) %>% 
  mutate(
    lpWGS_called = ifelse(snpID %in% maf0, "Y","N"),
    lpWGS_passfilter = ifelse(snpID %in% maf0_imp, "Y","N"),
    lpWGS_passfilter_maf5 = ifelse(snpID %in% maf5_imp, "Y","N"),
    lpWGS_passfilter_maf5_hpWGSconcord = ifelse(
      pct >= 60 & tot_calls > 30, "Y","N"),
    .after="snpID") %>% 
  mutate(pos = gsub(":[A-Z]:[A-Z]$","",snpID), .after="snpID")

#### Query ####
query <- c("1:247448734",
           "1:247440956",
           "2:32264782",
           "19:48234449",
           "19:48240212",
           "2:112836810",
           "11:112164735")

result <- dat %>% 
  filter(pos %in% query)
# result <- dat %>% 
#   filter(snpID %in% query)

write_csv(result, paste0("~/Desktop/",dt,"_RePORT_Brazil_SNP_query.csv"))

"
snpID Chromosome:position:reference:alternate
pos Chromosome:position
lpWGS_called Y/N if the SNP was called in lpWGS after imputation
lpWGS_passfilter Y/N if SNP passed imputation score filter (score > 0.6) in lpWGS
lpWGS_passfilter_maf5 Y/N if SNP passed imputation score filter, MAF > 5%, and HWE P < 1E-6
lpWGS_passfilter_maf5_hpWGSconcord Y/N if SNP passed imputation score filter, MAF > 5%, HWE P < 1E-6, and >= 60% concordant in at least 30 individuals in hpWGS data
tot_calls Total individuals with successful hpWGS call
pct Percent concordance between lpWGS and hpWGS
imp_score Imputation score
"
