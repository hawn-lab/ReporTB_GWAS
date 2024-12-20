#Single cell eQTL model fitting
library(tidyverse)
library(SNPRelate)
library(limma)
library(kimma)
source("scripts/models/reportb_eqtl.R")

#### Meta data ####
#Study metadata
load("../metadata/ReportTB_metadata_clean.RData")

#### SNP data ####
#Signif SNP
signif <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv")
#Filter for at least 1 protein coding annotation

biotypes <- colnames(signif)[grepl("biotype",colnames(signif))]
signif_pc <- signif %>% 
  filter(if_any(biotypes, ~.=="protein_coding"))

# Total SNP
nrow(signif_pc)
# Total genes
signif_pc %>% 
  select(contains("gene")) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  drop_na() %>% pull(value) %>% unique() %>% length()

GDS_file <- "result/gds/reportTB_sc_filter_imp_maf1.gds"
#Check if file exists
test.files <- list.files(path="result/gds/", 
                         pattern=".gds")

if(! "reportTB_sc_filter_imp_maf1.gds" %in% test.files){
  SNPRelate::snpgdsBED2GDS(bed.fn="result/plink/reportTB_sc_filter_imp_maf1.bed", 
                           fam.fn="result/plink/reportTB_sc_filter_imp_maf1.fam", 
                           bim.fn="result/plink/reportTB_sc_filter_imp_maf1.bim",
                           out.gdsfn=GDS_file)
}

#SNP data
genofile1 <- snpgdsOpen(GDS_file)
##List SNP in this gds file
snp.id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))
##Select signif SNP
snp.id1 <- snp.id1[snp.id1 %in% unique(signif_pc$snpID)]
geno.num <- snpgdsGetGeno(genofile1, snp.id=snp.id1)
## Add row and column names
colnames(geno.num) <- snp.id1
sample.id1 <- read.gdsn(index.gdsn(genofile1, "sample.id"))
rownames(geno.num) <- sapply(strsplit(sample.id1,"_"), `[`, 1)
snpgdsClose(genofile1)
rm(genofile1)

geno.num <- as.data.frame(geno.num) %>% 
  rownames_to_column("sample_id")

#### Single cell RNAseq data ####
load("data_clean/ReporTB_sc_voom.RData")

##### Define SNP-gene pairs ####
temp <- signif_pc %>% 
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
rm(genes.temp)

map <- temp %>% 
  inner_join(to_keep) %>% 
  filter(name=="gene") %>% 
  select(snpID,value) %>% 
  rename(hgnc_symbol=value) %>% 
  #keep only genes and snps in datasets
  filter(snpID %in% colnames(geno.num)) %>% 
  filter(hgnc_symbol %in% all_genes)
rm(temp)

#Pairs
nrow(map)
#SNPs
length(unique(map$snpID))
#Genes
length(unique(map$hgnc_symbol))
#Cells
length(ls(pattern = "voom_"))

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

#####
rm(geno.num, kin.raw)
gc()

#### Univariate models ####
#All HIV-
vooms <- ls(pattern = "voom_")
mods1 <- c("~genotype + PC1 + PC2 + (1|sample_id)",
          "~genotype + PC1 + PC2 + sex + (1|sample_id)",
          "~genotype + PC1 + PC2 + bl_age + (1|sample_id)",
          "~genotype + PC1 + PC2 + smokhx2 + (1|sample_id)",
          "~genotype + PC1 + PC2 + diabetes + (1|sample_id)")

models_univar <- reportb_eqtl(dat = vooms, models = mods1,
                              group = "univariate", processors = 6)

#### Leave one out models ####
mods2 <- c("~genotype + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes + (1|sample_id)",
          "~genotype + PC1 + PC2 + bl_age + smokhx2 + diabetes + (1|sample_id)",
          "~genotype + PC1 + PC2 + sex + smokhx2 + diabetes + (1|sample_id)",
          "~genotype + PC1 + PC2 + sex + bl_age + diabetes + (1|sample_id)",
          "~genotype + PC1 + PC2 + sex + bl_age + smokhx2 + (1|sample_id)")

models_loo <- reportb_eqtl(dat = vooms, models = mods2, 
                           group = "loo", processors = 6)

#Fix age naming error b/c column name is bl_age
models_loo <- models_loo %>% 
  mutate(model = recode(model, age="minus.age"))

#### Combine and save ####
save(models_univar, models_loo, file="result/sceQTL/eqtl_model_fitting.RData")

bind_rows(models_loo,models_univar) %>% 
  mutate(cell = gsub("_minus.bl$","", cell),
         cell = gsub(".bl$","", cell),
         cell = gsub("_bl$","", cell)) %>% 
  write_csv(file = "result/sceQTL/eqtl_model_fitting.csv")

beepr::beep()
