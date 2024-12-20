library(tidyverse)
library(openxlsx)

#### SNP lists ####
# Pass-filter concord SNP
impPF <- read_tsv("../ReporTB_lpWGS/result/snp_impute_pf06.txt", 
                  col_names=FALSE) %>% pull(X1) %>% sort()
maf5 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_imp_maf5.bim", 
                 col_names = FALSE) %>% pull(X2) %>% sort()
maf1_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf1_concord.txt",
                         col_names=FALSE) %>% pull(X1) %>% sort()
maf5_concord <- read_csv("../ReporTB_lpWGS/result/snp_imp_maf5_concord.txt", 
                         col_names=FALSE) %>% pull(X1) %>% sort()
length(maf5_concord)

#### Table S1. GWAS models ####
TabS1 <- list()
TabS1[["full cohort, final model"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>%
  select(snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  filter(snpID %in% maf5)

TabS1[["case-control, final model"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_cc.csv") %>% 
  select(snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  filter(snpID %in% maf5)

TabS1[["full cohort, no TB risks"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_full_comp.csv") %>% 
  select(snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp"))%>% 
  filter(snpID %in% maf5)

TabS1[["case-control, no TB risks"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_cc_comp.csv") %>% 
  select(snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  filter(snpID %in% maf5)

TabS1[["Significant SNP, univariate"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_univar.csv") %>%
  select(model, snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  mutate(model = recode(model, "plus_hiv"="plus_HIV","plus_diabetes"="plus_DM")) %>% 
  filter(snpID %in% maf5)

TabS1[["Significant SNP, leave-one-out"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_loo.csv") %>%
  select(model, snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  mutate(model = recode(model, "minus_hiv"="minus_HIV","minus_diabetes"="minus_DM")) %>% 
  filter(snpID %in% maf5)

TabS1[["AIC model fit"]] <- read_csv("../ReporTB_lpWGS/result/model_fitting/model_pcs.csv") %>%
  bind_rows(read_csv("../ReporTB_lpWGS/result/model_fitting/model_kin.csv")) %>% 
  bind_rows(read_csv("../ReporTB_lpWGS/result/model_fitting/model_cov.csv")) %>% 
  bind_rows(read_csv("../ReporTB_lpWGS/result/model_fitting/model_cov_minus1.csv")) %>% 
  filter(term == "value") %>% 
  select(model, snpID, AIC) %>% 
  pivot_wider(names_from = model, values_from = AIC)

TabS1[["smoke past vs current"]] <- read_csv("../ReporTB_lpWGS/result/model_fitting/model_cov.csv") %>%
  filter(model == "smoke") %>% 
  mutate(model = "smoke_3level") %>% 
  filter(term %in% c("smokhxpast","smokhxnever")) %>% 
  distinct(snpID, model, term, pval) %>% 
  pivot_wider(names_from = term, values_from = pval)

TabS1[["covariate outcomes"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_signif_cov.csv") %>% 
  select(model, snpID, term, estimate, pval, AIC, OR) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = recode(term, "sexF"="sex","bl_age"="age",
                       "smokhx2Y"="smoke","bl_hivY"="HIV", "diabetesY"="DM",
                       "value"="snp")) %>% 
  separate(model, into=c("outcome","subset"),sep="_") %>% 
  mutate(outcome = recode(outcome, "hiv"="HIV","diabetes"="DM")) %>% 
  filter(snpID %in% maf5)

# Imputation score and LD
impScore <- read_tsv("../ReporTB_lpWGS/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename.tsv.gz")
signif <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv")

temp=impScore %>% mutate(snpID=paste(`#CHROM`,POS,REF,ALT, sep=":")) %>% 
  inner_join(select(signif,snpID, full_pval))

biotypes <- colnames(signif)[grepl("biotype",colnames(signif))]
signif <- signif %>% 
  filter(if_any(biotypes, ~.=="protein_coding")) %>% 
  pull(snpID)

#Full cohort
full <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value") %>% 
  select(snpID, pval)
full2 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_full.txt") %>%
  filter(! SNP %in% full$snpID) %>% 
  rename(snpID=SNP,pval=PVAL) %>% 
  select(snpID, pval)

full_all <- bind_rows(full, full2) %>% 
  #MAF > 1%
  filter(snpID %in% maf1_concord)

maf1 <- data.frame()

for(snp in signif){
  print(snp)
  CHR <- strsplit(snp, split=":")[[1]][1]
  POS <- as.numeric(strsplit(snp, split=":")[[1]][2])
  
  #LD
  LD_name <- gsub(":","_",snp)
  LD_name <- paste0("../ReporTB_lpWGS/result/LD/",LD_name,".RData")
  load(LD_name)
  temp0 <- LD.df %>% 
    rename(signifSNP=lead.snp)
  
  maf1 <- full_all %>% 
    select(snpID, pval) %>% 
    inner_join(temp0, by = join_by(snpID)) %>% 
    mutate(signifSNP=snp) %>% 
    select(signifSNP, snpID, pval, R2) %>% 
    bind_rows(maf1)
}

TabS1[["full cohort, maf1%, LD"]] <- maf1

write.xlsx(TabS1, file="TableS1_gwas_models.xlsx")

#### Table S5. signif SNP ####
TabS5 <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  select(snpID, rsID:cc_pval, full_OR:MAF_case, contains("gene")) %>% 
  select(where(~ !all(is.na(.))))

#add concordance
concord <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz") %>% 
  filter(snpLow %in% TabS5$snpID)
#add impute score
imp <- read_tsv("../ReporTB_lpWGS/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv.gz") %>% 
  filter(snpID %in% TabS5$snpID)

TabS5 <- TabS5 %>% 
  left_join(concord %>% select(snpLow, tot_calls, pct) %>% 
              rename(hpWGS_calls=tot_calls, hpWGS_concordance=pct), 
            by=c("snpID"="snpLow")) %>% 
  left_join(imp %>% select(snpID, INFO) %>% rename(impute_score=INFO), 
            by="snpID") %>% 
  select(snpID:rsID, impute_score, hpWGS_calls, hpWGS_concordance,
         everything())

write_csv(TabS5, "TableS5_signif_snp.csv", na="")

#### Table S3. concordance ####
concord <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz")

TabS3 <- concord %>% 
  select(snpLow, tot_calls, pct) %>% 
  filter(tot_calls>0) %>% 
  rename(snpID=snpLow, total_call=tot_calls, percent_concordance=pct)

write_csv(TabS3, file="TableS3_hpWGS_concordance.csv", na="")

#### Table S6. enrichment ####
load("../ReporTB_lpWGS/result/enrich/enrich_iter_signif_concord_all.RData")

TabS6 <- enrich_pc %>% 
  filter(FDR < 0.2 & k_median > 1) %>% 
  mutate(database = paste(gs_cat,gs_subcat,sep=":")) %>% 
  select(database, pathway, pvalue_median, FDR, k_median, K, genes) %>% 
  mutate(genes=as.character(genes))

write_csv(TabS6, file="TableS6_enrichment.csv", na="")

#### Table S7. Array overlap ####
overlap <- read.csv("../ReporTB_lpWGS/result/liftover/compare_to_SNPchip.csv") %>% 
  select(-none) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  filter(value != "") %>% 
  distinct(name, value, snpID) %>% 
  mutate(name = recode(name,
                       "af_limaa"="Affymetrix LIMAArray",
                       "il_ominiZhongHua"="Illumina Omni ZhongHua",
                       "il_omni2.5"="Illumina Omni 2.5",
                       "il_omni5"="Illumina Omni 5",
                       "il_megaex38"="Illumina MegaEx38",
                       "il_660Quad"="Illumina 660Quad",
                       "il_omni1"="Illumina Omni 1",
                       "af_100k"="Affymetrix 100K",
                       "af_500k"="Affymetrix 500K",
                       "af_6.0"="Affymetrix 6.0",
                       "af_axiompel"="Affymetrix AxiomPEL",
                       "af_chb1"="Affymetrix CHB1",
                       "af_chb2"="Affymetrix CHB2",
                       "il_1M"="Illumina 1M",
                       "il_610Quad"="Illumina 610Quad",
                       "il_core_exome"="Illumina Core Exome",
                       "il_omniExpress0"="Illumina Omni Express 0",
                       "il_omniExpress1"="Illumina Omni Express 1",
                       "il_omniExpress2"="Illumina Omni Express 2")) %>% 
  arrange(name, snpID) %>% 
  rename("SNP array"=name, "SNP (this study)"=snpID,
         "SNP (on array)"=value) 

length(unique(overlap$`SNP array`))
length(unique(overlap$`SNP (this study)`))

overlap %>% 
  drop_na() %>% 
  count(`SNP array`) %>% 
  arrange(-n)

write_csv(overlap, "temp/TableS7_array_overlap.csv", na="")

#### Table S2. sceQTL ####
attach("../ReporTB_sceQTL/result/sceQTL/eQTL_model_fitting.RData")
# attach("../ReporTB_sceQTL/result/sceQTL/sceQTL_noCov_df.RData")
# attach("../ReporTB_sceQTL/result/sceQTL/sceQTL_noCov_no2_df.RData")
attach("../ReporTB_sceQTL/result/sceQTL/sceQTL_df.RData")
TabS2 <- list()

#snpID	term	estimate	pval	AIC	OR
TabS2[["final model"]] <- model_noCov %>% 
  full_join(model_noCov_fit %>% 
              select(cell, condition, genotype, gene, AIC)) %>% 
  mutate(snpID = gsub("snp[.]","",genotype),
         snpID = gsub("_",":",snpID)) %>% 
  select(condition, cell, gene, snpID, variable_group, 
         estimate, pval, FDR, AIC) %>% 
  rename(term=variable_group) %>% 
  mutate(term = recode(term, "SNP"="snp")) %>% 
  mutate(cell = gsub("_"," ",cell))

TabS2[["with covariates"]] <- model_lmerel %>%
  full_join(model_fit %>% 
              select(cell, condition, genotype, gene, AIC)) %>% 
  mutate(snpID = gsub("snp[.]","",genotype),
         snpID = gsub("_",":",snpID)) %>% 
  select(condition, cell, gene, snpID, variable_group, 
         estimate, pval, FDR, AIC) %>% 
  rename(term=variable_group) %>% 
  mutate(term = recode(term, "bl_age"="age",
                       "smokhx2"="smoke","diabetes"="DM",
                       "SNP"="snp")) %>% 
  mutate(cell = gsub("_"," ",cell))

#snpID	model	AIC
TabS2[["univariate fitting"]] <- models_univar %>% 
  mutate(snpID = gsub("snp[.]","",genotype),
         snpID = gsub("_",":",snpID)) %>% 
  mutate(condition="media") %>% 
  select(condition, cell, gene, snpID, model, AIC) %>% 
  mutate(cell = gsub("_"," ",cell))

TabS2[["leave-one-out fitting"]] <- models_loo %>% 
  mutate(snpID = gsub("snp[.]","",genotype),
         snpID = gsub("_",":",snpID)) %>% 
  mutate(condition="media") %>% 
  select(condition, cell, gene, snpID, model, AIC)  %>% 
  mutate(cell = gsub("_"," ",cell))

TabS2[["name key"]] <- read_csv("../ReporTB_sceQTL/result/sc_key.csv") %>% 
  mutate(anno_long = gsub("cluster ","",anno_long))

write.xlsx(TabS2, file="TableS2_sceQTL_models.xlsx")

#### Table 2 Signif sceQTL ####
model_noCov %>% 
  #signif
  filter(FDR < 0.05 & variable_group=="SNP") %>% 
  mutate(snpID = gsub("snp[.]","",genotype),
         snpID = gsub("_",":",snpID)) %>% 
  select(condition, cell, gene, snpID, 
         estimate, FDR) %>% 
  arrange(snpID) %>% 
  write_csv("temp/Table2_sceQTL.csv")
