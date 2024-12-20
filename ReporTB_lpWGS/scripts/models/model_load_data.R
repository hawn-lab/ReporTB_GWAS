library(tidyverse)
library(SNPRelate)
library(GMMAT)

#### Metadata ####
attach("~/project/data/metadata/ReportTB_metadata_clean.RData")
attach("~/project/result/kinship/reportTB_PCair.RData")

pc <- as.data.frame(PC.dat) %>%
  rownames_to_column("sample_id") %>% 
  separate(sample_id, into=c("sample_id"), sep="_", extra="drop")

#Add PCs and format
meta_FULL_PC <- meta_FULL %>% 
  inner_join(pc) %>% 
  #recode tb
  mutate(tb_recode = recode(tb, "contact"="0", "culture_confirm"="1"),
         tb_recode = as.numeric(tb_recode)) %>% 
  arrange(sample_id) %>% 
  mutate(smokhx2 = case_when(smokhx %in% c("current","past")~"Y",
                             smokhx == "never"~"N",
                             TRUE~smokhx),
         smokhx2 = factor(smokhx2, levels=c("N","Y")))

meta_CC_PC <- meta_CC %>% 
  inner_join(pc) %>% 
  #recode tb
  mutate(tb_recode = recode(tb, "contact"="0", "culture_confirm"="1"),
         tb_recode = as.numeric(tb_recode)) %>% 
  arrange(sample_id) %>% 
  mutate(smokhx2 = case_when(smokhx %in% c("current","past")~"Y",
                             smokhx == "never"~"N",
                             TRUE~smokhx),
         smokhx2 = factor(smokhx2, levels=c("N","Y")))

#### Genotypes ####
if(MAF==5){
  # MAF 5% and HWE 1E-6 filtered SNPs
  GDS_file <- "~/project/result/gds/reportTB_filter_imp_maf5.gds"
  
  #Check if file exists
  test.files <- list.files(path="~/project/result/gds/", 
                           pattern=".gds")
  
  if(! "reportTB_filter_imp_maf5.gds" %in% test.files){
    SNPRelate::snpgdsBED2GDS(bed.fn="~/project/result/plink/reportTB_filter_imp_maf5.bed", 
                             fam.fn="~/project/result/plink/reportTB_filter_imp_maf5.fam", 
                             bim.fn="~/project/result/plink/reportTB_filter_imp_maf5.bim",
                             out.gdsfn=GDS_file)
  }
} else if(MAF==1){
  # MAF 1% and HWE 1E-6 filtered SNPs
  GDS_file <- "~/project/result/gds/reportTB_filter_imp_maf1.gds"
  
  #Check if file exists
  test.files <- list.files(path="~/project/result/gds/", 
                           pattern=".gds")
  
  if(! "reportTB_filter_imp_maf1.gds" %in% test.files){
    SNPRelate::snpgdsBED2GDS(bed.fn="~/project/result/plink/reportTB_filter_imp_maf1.bed", 
                             fam.fn="~/project/result/plink/reportTB_filter_imp_maf1.fam", 
                             bim.fn="~/project/result/plink/reportTB_filter_imp_maf1.bim",
                             out.gdsfn=GDS_file)
  }
}

# Extract and format genotypes
genofile <- snpgdsOpen(GDS_file)
geno.num <- snpgdsGetGeno(genofile)

## Add row and column names
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
colnames(geno.num) <- snp.id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
# rownames(geno.num) <- sapply(strsplit(sample.id,"_"), `[`, 1)
rownames(geno.num) <- sample.id

## Add SNP library name to metadata
meta_FULL_PC <- data.frame(libID = sample.id) %>% 
  separate(libID, into = c("sample_id", "dragen"), sep="_", remove = FALSE) %>% 
  inner_join(meta_FULL_PC)

meta_CC_PC <- data.frame(libID = sample.id) %>% 
  separate(libID, into = c("sample_id", "dragen"), sep="_", remove = FALSE) %>% 
  inner_join(meta_CC_PC)

# Runs out of RAM, Transpose after filtering by score instead
# geno.num.t <- t(geno.num)
# geno.num.t <- WGCNA::transposeBigData(geno.num, blocksize = 100)

snpgdsClose(genofile)

#### Kinship data ####
kin.raw <- read_csv("~/project/result/kinship/reportTB_kin.csv")

kin.matrix <- kin.raw %>% 
  #rename to match SNP data
  inner_join(meta_FULL_PC %>% select(sample_id, libID), 
             by = c("rowname"="sample_id")) %>% 
  select(-rowname) %>% 
  dplyr::rename(rowname=libID) %>% 
  pivot_longer(-rowname) %>% 
  inner_join(meta_FULL_PC %>% select(sample_id, libID),
             by = c("name"="sample_id")) %>% 
  select(-name) %>% 
  pivot_wider(names_from = libID) %>% 
  #order by geno data
  mutate(rowname = factor(rowname, levels=rownames(geno.num))) %>% 
  arrange(rowname) %>% 
  column_to_rownames() %>% 
  select(all_of(rownames(geno.num))) %>% 
  as.matrix()

check1 <- identical(rownames(kin.matrix),colnames(kin.matrix))
check2 <- identical(rownames(kin.matrix), rownames(geno.num))

if(!check1 | !check2){
  stop("Kinship and genotypes do not contain the same samples in the same order.")
}

#### Case-control subset ####
geno.num.cc <- geno.num[meta_CC_PC$libID,]

kin.matrix_CC <- kin.raw %>% 
  #rename to match SNP data
  inner_join(meta_CC_PC %>% select(sample_id, libID), 
             by = c("rowname"="sample_id")) %>% 
  select(-rowname) %>% 
  dplyr::rename(rowname=libID) %>% 
  pivot_longer(-rowname) %>% 
  inner_join(meta_CC_PC %>% select(sample_id, libID), 
             by = c("name"="sample_id")) %>% 
  select(-name) %>% 
  pivot_wider(names_from = libID) %>% 
  #order by geno data
  mutate(rowname = factor(rowname, levels=rownames(geno.num.cc))) %>% 
  arrange(rowname) %>% 
  column_to_rownames() %>% 
  select(all_of(rownames(geno.num.cc))) %>% 
  as.matrix()

check3 <- identical(rownames(kin.matrix_CC),colnames(kin.matrix_CC))
check4 <- identical(rownames(kin.matrix_CC), rownames(geno.num.cc))

if(!check3 | !check4){
  stop("Case-control kinship and genotypes do not contain the same samples in the same order.")
}

print("Data load complete.")
