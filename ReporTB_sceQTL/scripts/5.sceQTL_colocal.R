#Single cell eQTL analysis
library(tidyverse)
library(SNPRelate)
library(limma)
library(kimma)
# library(BIGpicture)
# library(patchwork)

#### Signif eQTL ####
attach("result/sceQTL/sceQTL_df.RData")

signif_media <- model_noCov %>% 
  filter(variable_group == "SNP" & FDR < 0.05 &
           condition=="media") %>% 
  distinct(cell,genotype,gene) %>% 
  separate(genotype, into=c("trash","chr","pos"), sep="[.]|_", 
           remove = FALSE, extra="drop")

#### Find surrounding SNP ####
GDS_file <- "result/gds/reportTB_sc_filter_imp_maf1.gds"
#SNP data
genofile1 <- snpgdsOpen(GDS_file)
##List SNP in this gds file
snp.id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))

#Chr in signif
signif_chr <- unique(signif_media$chr) %>% paste0("^",.,":") %>% 
  paste(., collapse="|")

# +/- 10 kb from signif snp
signif_span <- signif_media %>% 
  mutate(min.POS = as.numeric(pos)-10E3,
         max.POS = as.numeric(pos)+10E3) %>% 
  distinct(genotype, chr, max.POS, min.POS) %>% 
  rename(lead.snp=genotype)

maf1_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf1_concord.txt",
                         col_names=FALSE) %>% pull(X1)

snp.regions.signif <- data.frame(genotype=snp.id1) %>% 
  filter(genotype %in% maf1_concord) %>% 
  filter(grepl(signif_chr, genotype)) %>% 
  separate(genotype, into=c("chr","pos"), sep=":", 
           remove = FALSE, extra="drop") %>% 
  mutate(pos=as.numeric(pos)) %>% 
  inner_join(signif_span,
             by="chr", relationship = "many-to-many") %>% 
  filter(pos >= min.POS & pos <= max.POS)

#Recode everything to leade ZNF SNP
snp.regions.signif <- snp.regions.signif %>% 
  mutate(lead.snp = ifelse(grepl("^snp.3_", lead.snp) &
                             lead.snp != "snp.3_75816280_A_T", 
                           "snp.3_75816280_A_T",
                           lead.snp)) %>% 
  distinct(genotype, chr, pos, lead.snp)

map <- signif_media %>% 
  select(genotype, gene)

#Fix names in gds
signif_media_snp <- gsub("snp.","", unique(snp.regions.signif$genotype)) %>% 
  gsub("_",":",.)

#### Meta data ####
#Study metadata
load("../metadata/ReportTB_metadata_clean.RData")

#### SNP data ####
##Select signif SNP and surrounding region
snp.id1 <- snp.id1[snp.id1 %in% signif_media_snp]
geno.num <- snpgdsGetGeno(genofile1, snp.id=snp.id1)
## Add row and column names
colnames(geno.num) <- snp.id1
sample.id1 <- read.gdsn(index.gdsn(genofile1, "sample.id"))
rownames(geno.num) <- sapply(strsplit(sample.id1,"_"), `[`, 1)
snpgdsClose(genofile1)

geno.num <- as.data.frame(geno.num) %>% 
  rownames_to_column("sample_id")

#### Single cell RNAseq data ####
load("data_clean/ReporTB_sc_voom.RData")

#### Check donor overlap ####
voom_B$targets %>%
  rownames_to_column() %>% 
  distinct(sample_id, TB) %>% 
  mutate(scRNAseq = "Y") %>% 
  mutate(SNP = ifelse(sample_id %in% geno.num$sample_id, "Y", "N")) %>% 
  arrange(TB, sample_id) %>% 
  count(TB, scRNAseq, SNP)

#All donors present in both datasets

#### Rename SNP ####
#Avoid ":" because kmFit uses that to split interaction variables
geno.num.rename <- geno.num %>% 
  rename_with(~gsub(":","_", .), everything()) %>% 
  rename_with(~paste0("snp.",.), !contains("sample_id")) 

#### Kinship data ####
kin.raw <- read_csv("result/kinship/reportTB_sc_kin.csv")

kin.matrix <- kin.raw %>% 
  #order by geno data
  mutate(rowname = factor(rowname, levels=geno.num.rename$sample_id)) %>% 
  arrange(rowname) %>% 
  column_to_rownames() %>% 
  select(all_of(geno.num.rename$sample_id)) %>% 
  as.matrix()

#### NO COVARIATES ####
colocal_media_noCov_ls <- list()
for(voom in ls(pattern = "voom_")){
  #get snp with signif in this cell type
  cellOI <- gsub("voom_","",voom)
  #Deal with multiple SNPs in ZNF
  if(voom %in% c("voom_CD14_Mono-1","voom_CD4_Treg")){
    to_model <- signif_media %>% 
      filter(gene != "ZNF717" | genotype=="snp.3_75816280_A_T") %>% 
      filter(cell==cellOI)
  } else {
    to_model <- signif_media %>% 
      filter(cell==cellOI)
  }
  if(nrow(to_model>0)){
    print(voom)
    voom.temp <- get(voom)
    
    #Media only samples
    ##meta
    voom.temp$targets <- voom.temp$targets %>%
      filter(condition=="Media")
    ##weights
    colnames(voom.temp$weights) <- colnames(voom.temp$E)
    rownames(voom.temp$weights) <- rownames(voom.temp$E)
    voom.temp$weights <- as.data.frame(voom.temp$weights) %>%
      select(all_of(voom.temp$targets$libID)) %>%
      as.matrix()
    ##expression
    voom.temp$E <- as.data.frame(voom.temp$E) %>%
      select(all_of(voom.temp$targets$libID)) %>%
      as.matrix()
    
    for(i in 1:nrow(to_model)){
      to_model_temp <- to_model[i,]
      
      genotemp <- to_model_temp$genotype
      if(genotemp %in% c("snp.3_75816342_A_C","snp.3_75816310_C_A")){
        snp.regions.signif.temp <- snp.regions.signif %>% 
          mutate(lead.snp = recode(lead.snp,
                                   "snp.3_75816280_A_T"=genotemp))
      } else {
        snp.regions.signif.temp <- snp.regions.signif
      }
      #Format map
      map.rename.temp <- snp.regions.signif.temp %>% 
        inner_join(to_model_temp %>% select(genotype, gene),
                   by=c("lead.snp"="genotype")) %>% 
        mutate(snpID = paste("snp",genotype,sep="."),
               snpID = gsub(":","_",snpID)) %>% 
        distinct(snpID, gene)
      
      #Rm genotypes with only 1 level and without matching gene
      geno.num.rename.temp <- geno.num.rename %>%
        filter(sample_id %in% voom.temp$targets$sample_id) %>% 
        select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
        select_if(function(col) length(unique(col))>1) 
      
      #Run models
      colocal_media_noCov_ls[[paste0(voom,".media",".",to_model_temp$gene)]] <-
        kmFit_eQTL(
          dat_snp = geno.num.rename.temp, 
          dat = voom.temp,
          dat_map = map.rename.temp,
          kin = kin.matrix,
          geneID="gene", genotypeID="snpID",
          patientID="sample_id",
          model="~genotype + PC1 + PC2 + (1|sample_id)",
          use_weights = TRUE, run_lmerel=TRUE, run_contrast = FALSE,
          metrics=TRUE, processors=1)
    }
  }
  
}

#Format to df
colocal_noCov <- data.frame()

for(m in names(colocal_media_noCov_ls)){
  model.temp <- colocal_media_noCov_ls[[m]]
  
  colocal_noCov <- model.temp$lmerel %>% 
    mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
    separate(cell, into=c("cell","condition","gene"), sep="[.]") %>% 
    #Fill pairs with only 1 SNP level
    mutate(estimate = ifelse(estimate=="seeContrasts",NA, estimate),
           statistic = ifelse(statistic=="seeContrasts",NA, statistic)) %>% 
    mutate(estimate = as.numeric(estimate),
           statistic = as.numeric(statistic)) %>% 
    bind_rows(colocal_noCov)
}

#Recalc FDR
colocal_noCov <- colocal_noCov %>% 
  rename(FDR_orig=FDR) %>% 
  mutate(variable_group = case_when(grepl("^snp", variable)~"SNP",
                                    TRUE~variable)) %>% 
  group_by(condition, cell, variable_group, genotype) %>% 
  mutate(FDR = p.adjust(pval, method="BH"), .before=FDR_orig) %>% 
  ungroup()

#### Save ####
save(colocal_noCov, colocal_media_noCov_ls, file="result/sceQTL/colocal.RData")

beepr::beep()
