#Single cell eQTL analysis
library(tidyverse)
library(SNPRelate)
library(limma)
library(kimma)
# library(BIGpicture)
# library(patchwork)

#### Meta data ####
#Study metadata
load("../metadata/ReportTB_metadata_clean.RData")

#### SNP data ####
#Signif SNP
signif <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv")
#Filter for at least 1 protein coding annotation

biotypes <- colnames(signif)[grepl("biotype",colnames(signif))]
signif_pc_concord <- signif %>% 
  filter(if_any(biotypes, ~.=="protein_coding"))

# Total SNP
nrow(signif_pc_concord)
# Total genes
all_pc <- signif_pc_concord %>% 
  select(snpID,contains("biotype"),contains("gene")) %>% 
  pivot_longer(-snpID) %>% 
  separate(name, into=c("cis","name","ID"), sep="_") %>% 
  pivot_wider() %>% 
  filter(biotype == "protein_coding") %>% 
  pull(gene) %>% unique()

signif_pc_concord %>% 
  select(contains("gene")) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  filter(value %in% all_pc) %>% 
  drop_na() %>% pull(value) %>% unique() %>% length()

GDS_file <- "result/gds/reportTB_sc_filter_imp_maf1.gds"
#SNP data
genofile1 <- snpgdsOpen(GDS_file)
##List SNP in this gds file
snp.id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))
##Select signif SNP
snp.id1 <- snp.id1[snp.id1 %in% unique(signif_pc_concord$snpID)]
geno.num <- snpgdsGetGeno(genofile1, snp.id=snp.id1)
## Add row and column names
colnames(geno.num) <- snp.id1
sample.id1 <- read.gdsn(index.gdsn(genofile1, "sample.id"))
rownames(geno.num) <- sapply(strsplit(sample.id1,"_"), `[`, 1)
snpgdsClose(genofile1)

geno.num <- as.data.frame(geno.num) %>% 
  rownames_to_column("sample_id")

# Save
write_csv(geno.num, file="data_clean/geno.num.subset.csv")

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

##### Define SNP-gene pairs ####
temp <- signif_pc_concord %>% 
  select(snpID, intragenic_gene_1:cis_distance_30) %>% 
  select(!contains("distance") & !contains("ensembl")) %>% 
  pivot_longer(-snpID) %>% 
  separate(name, into=c("anno_type","name","nameID"), sep="_")

to_keep <- temp %>%
  filter(value=="protein_coding") %>% 
  select(snpID, anno_type, nameID)

all_genes <- c()
for(voom in ls(pattern = "voom_")){
  genes.temp <- get(voom)
  genes.temp <- genes.temp$genes$hgnc_symbol
  all_genes <- c(all_genes,genes.temp) %>% unique()
}

map <- temp %>% 
  inner_join(to_keep) %>% 
  filter(name=="gene") %>% 
  select(snpID,value) %>% 
  rename(hgnc_symbol=value) %>% 
  #keep only genes and snps in datasets
  filter(snpID %in% colnames(geno.num)) %>% 
  filter(hgnc_symbol %in% all_genes)

#SNP be gene
map %>% 
  write_csv("result/sceQTL/sceQTL_map.csv")

map %>% 
  count(hgnc_symbol) %>% 
  write_csv("result/sceQTL/sceQTL_snp_per_gene.csv")

map %>% 
  count(snpID) %>% 
  arrange(snpID) %>% 
  write_csv("result/sceQTL/sceQTL_gene_per_snp.csv")

nrow(map)
map %>% pull(hgnc_symbol) %>% unique() %>% length()
map %>% pull(snpID) %>% unique() %>% length()

#### Rename SNP ####
#Avoid ":" because kmFit uses that to split interaction variables
geno.num.rename <- geno.num %>% 
  select(sample_id, all_of(map$snpID)) %>%
  rename_with(~gsub(":","_", .), everything()) %>% 
  rename_with(~paste0("snp.",.), !contains("sample_id")) 

map.rename <- map %>% 
  mutate(snpID = gsub(":","_",snpID),
         snpID = paste0("snp.",snpID))

#### Kinship data ####
kin.raw <- read_csv("result/kinship/reportTB_sc_kin.csv")

kin.matrix <- kin.raw %>% 
  #order by geno data
  mutate(rowname = factor(rowname, levels=geno.num.rename$sample_id)) %>% 
  arrange(rowname) %>% 
  column_to_rownames() %>% 
  select(all_of(geno.num.rename$sample_id)) %>% 
  as.matrix()

#### ALL COVARIATES ####
#All HIV-
#### Media only model ####
model_media_ls <- list()
for(voom in ls(pattern = "voom_")){
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
  
  #Rm genes not in expression. Differs for indiv cell types
  map.rename.temp <- map.rename %>% 
    filter(hgnc_symbol %in% rownames(voom.temp$E))
  
  #Rm genotypes with only 1 level and without matching gene
  geno.num.rename.temp <- geno.num.rename %>%
    filter(sample_id %in% voom.temp$targets$sample_id) %>% 
    select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
    select_if(function(col) length(unique(col))>1) 
  
  #Run models
  model_media_ls[[paste0(voom,".media")]] <- kmFit_eQTL(
    dat_snp = geno.num.rename.temp, 
    dat = voom.temp,
    dat_map = map.rename.temp,
    kin = kin.matrix,
    geneID="hgnc_symbol", genotypeID="snpID",
    patientID="sample_id",
    model="~genotype + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes + (1|sample_id)",
    use_weights = TRUE, run_lmerel=TRUE, run_contrast = FALSE,
    metrics=TRUE, processors=1)
}

# Format to df
model_lmerel <- data.frame()
model_fit <- data.frame()

for(m in names(model_media_ls)){
  model.temp <- model_media_ls[[m]]
  
  model_lmerel <- model.temp$lmerel %>% 
    mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
    separate(cell, into=c("cell","condition"), sep="[.]", extra="merge") %>% 
    #Fill pairs with only 1 SNP level
    mutate(estimate = ifelse(estimate=="seeContrasts",NA, estimate),
           statistic = ifelse(statistic=="seeContrasts",NA, statistic)) %>% 
    mutate(estimate = as.numeric(estimate),
           statistic = as.numeric(statistic)) %>% 
    bind_rows(model_lmerel)
  
  model_fit <- model.temp$lmerel.fit %>% 
    mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
    separate(cell, into=c("cell","condition"), sep="[.]", extra="merge") %>% 
    bind_rows(model_fit)
}

#Recalc FDR
model_lmerel <- model_lmerel %>% 
  rename(FDR_orig=FDR) %>% 
  mutate(variable_group = case_when(grepl("^snp", variable)~"SNP",
                                    TRUE~variable)) %>% 
  group_by(condition, cell, variable_group, genotype) %>% 
  mutate(FDR = p.adjust(pval, method="BH"), .before=FDR_orig) %>% 
  ungroup()

#### NO COVARIATES ####
model_media_noCov_ls <- list()
for(voom in ls(pattern = "voom_")){
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
  
  #Rm genes not in expression. Differs for indiv cell types
  map.rename.temp <- map.rename %>% 
    filter(hgnc_symbol %in% rownames(voom.temp$E))
  
  #Rm genotypes with only 1 level and without matching gene
  geno.num.rename.temp <- geno.num.rename %>%
    filter(sample_id %in% voom.temp$targets$sample_id) %>% 
    select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
    select_if(function(col) length(unique(col))>1) 
  
  #Run models
  model_media_noCov_ls[[paste0(voom,".media")]] <- kmFit_eQTL(
    dat_snp = geno.num.rename.temp, 
    dat = voom.temp,
    dat_map = map.rename.temp,
    kin = kin.matrix,
    geneID="hgnc_symbol", genotypeID="snpID",
    patientID="sample_id",
    model="~genotype + PC1 + PC2 + (1|sample_id)",
    use_weights = TRUE, run_lmerel=TRUE, run_contrast = FALSE,
    metrics=TRUE, processors=1)
}

#Format to df
model_noCov <- data.frame()
model_noCov_fit <- data.frame()

for(m in names(model_media_noCov_ls)){
  model.temp <- model_media_noCov_ls[[m]]
  
  model_noCov <- model.temp$lmerel %>% 
    mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
    separate(cell, into=c("cell","condition"), sep="[.]", extra="merge") %>% 
    #Fill pairs with only 1 SNP level
    mutate(estimate = ifelse(estimate=="seeContrasts",NA, estimate),
           statistic = ifelse(statistic=="seeContrasts",NA, statistic)) %>% 
    mutate(estimate = as.numeric(estimate),
           statistic = as.numeric(statistic)) %>% 
    bind_rows(model_noCov)
  
  model_noCov_fit <- model.temp$lmerel.fit %>% 
    mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
    separate(cell, into=c("cell","condition"), sep="[.]", extra="merge") %>% 
    bind_rows(model_noCov_fit)
}

#Recalc FDR
model_noCov <- model_noCov %>% 
  rename(FDR_orig=FDR) %>% 
  mutate(variable_group = case_when(grepl("^snp", variable)~"SNP",
                                    TRUE~variable)) %>% 
  group_by(condition, cell, variable_group, genotype) %>% 
  mutate(FDR = p.adjust(pval, method="BH"), .before=FDR_orig) %>% 
  ungroup()

#### Save ####
save(model_media_ls, model_media_noCov_ls,
     file="result/sceQTL/sceQTL_lists.RData")

save(model_lmerel, model_fit,
     model_noCov, model_noCov_fit,
     file="result/sceQTL/sceQTL_df.RData")

beepr::beep()

#### Confirm hits without "2" genotypes ####
#Not applicable in file data set
# model_media_noCov_ls2 <- list()
# for(voom in ls(pattern = "voom_")){
#   print(voom)
#   voom.temp <- get(voom)
#   #Media only samples
#   ##meta
#   voom.temp$targets <- voom.temp$targets %>%
#     filter(condition=="Media")
#   ##weights
#   colnames(voom.temp$weights) <- colnames(voom.temp$E)
#   rownames(voom.temp$weights) <- rownames(voom.temp$E)
#   voom.temp$weights <- as.data.frame(voom.temp$weights) %>%
#     select(all_of(voom.temp$targets$libID)) %>%
#     as.matrix()
#   ##expression
#   voom.temp$E <- as.data.frame(voom.temp$E) %>%
#     select(all_of(voom.temp$targets$libID)) %>%
#     as.matrix()
#   
#   #Rm genes not in expression. Differs for indiv cell types
#   map.rename.temp <- map.rename %>% 
#     filter(hgnc_symbol %in% rownames(voom.temp$E))
#   
#   #Rm genotypes with only 1 level and without matching gene
#   #Select genotypes with at least one 2 
#   #Remove indiv with 2 genotype
#   geno.num.rename.temp <- geno.num.rename %>%
#     filter(sample_id %in% voom.temp$targets$sample_id) %>% 
#     select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
#     select_if(function(col) length(unique(col))>1) %>% 
#     column_to_rownames("sample_id") %>% 
#     select_if(function(col) max(col)==2) %>% 
#     rownames_to_column("sample_id") %>% 
#     mutate(across(where(is.numeric), ~ifelse(.==2,NA, .)))
#   #Alternatively make 2->1
#   geno.num.rename.temp2 <- geno.num.rename %>%
#     filter(sample_id %in% voom.temp$targets$sample_id) %>% 
#     select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
#     select_if(function(col) length(unique(col))>1) %>% 
#     column_to_rownames("sample_id") %>% 
#     select_if(function(col) max(col)==2) %>% 
#     rownames_to_column("sample_id") %>% 
#     mutate(across(where(is.numeric), ~ifelse(.==2,1, .)))
#   
#   #Run models
#   if(ncol(geno.num.rename.temp)>1){
#     model_media_noCov_ls2[[paste0(voom,".rm2.media")]] <- kmFit_eQTL(
#       dat_snp = geno.num.rename.temp, 
#       dat = voom.temp,
#       dat_map = map.rename.temp,
#       kin = kin.matrix,
#       geneID="hgnc_symbol", genotypeID="snpID",
#       patientID="sample_id",
#       model="~genotype + PC1 + PC2 + (1|sample_id)",
#       use_weights = TRUE, run_lmerel=TRUE, run_contrast = FALSE,
#       metrics=TRUE, processors=1)
#     
#     model_media_noCov_ls2[[paste0(voom,".dom.media")]] <- kmFit_eQTL(
#       dat_snp = geno.num.rename.temp2, 
#       dat = voom.temp,
#       dat_map = map.rename.temp,
#       kin = kin.matrix,
#       geneID="hgnc_symbol", genotypeID="snpID",
#       patientID="sample_id",
#       model="~genotype + PC1 + PC2 + (1|sample_id)",
#       use_weights = TRUE, run_lmerel=TRUE, run_contrast = FALSE,
#       metrics=TRUE, processors=1)
#     
#     }
# }
# 
# #Format to df
# model_noCov2 <- data.frame()
# model_fit2 <- data.frame()
# 
# for(m in names(model_media_noCov_ls2)){
#   model.temp <- model_media_noCov_ls2[[m]]
#   
#   model_noCov2 <- model.temp$lmerel %>% 
#     mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
#     separate(cell, into=c("cell","model","condition"), sep="[.]", extra="merge") %>% 
#     #Fill pairs with only 1 SNP level
#     mutate(estimate = ifelse(estimate=="seeContrasts",NA, estimate),
#            statistic = ifelse(statistic=="seeContrasts",NA, statistic)) %>% 
#     mutate(estimate = as.numeric(estimate),
#            statistic = as.numeric(statistic)) %>% 
#     bind_rows(model_noCov2)
#   
#   model_fit2 <- model.temp$lmerel.fit %>% 
#     mutate(cell=gsub("voom_","",m), .before="genotype") %>% 
#     separate(cell, into=c("cell","model","condition"), sep="[.]", extra="merge") %>% 
#     bind_rows(model_fit2)
# }
# 
# #Recalc FDR
# model_noCov2 <- model_noCov2 %>% 
#   rename(FDR_orig=FDR) %>% 
#   mutate(variable_group = case_when(grepl("^snp", variable)~"SNP",
#                                     TRUE~variable)) %>% 
#   group_by(model, condition, cell, variable_group, genotype) %>% 
#   mutate(FDR = p.adjust(pval, method="BH"), .before=FDR_orig) %>% 
#   ungroup()

#### Save ####
save(model_media_ls, model_media_noCov_ls, #model_media_noCov_ls2,
     file="result/sceQTL/sceQTL_lists.RData")

save(model_lmerel, model_fit,
     model_noCov, model_noCov_fit,
     #model_noCov2, model_fit2,
     file="result/sceQTL/sceQTL_df.RData")

beepr::beep()
