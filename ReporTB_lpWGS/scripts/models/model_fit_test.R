library(tidyverse)
library(SNPRelate)
# library(SNPassoc)
library(GMMAT)
library(doParallel)
library(broom)
library(foreach)

dir.create("~/project/result/model_fitting/", showWarnings = FALSE)

print("Data")

#### Fxns ####
source("~/project/scripts/models/glmmkin_summ.R")
source("~/project/scripts/models/glm_summ.R")

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
         tb_recode = as.numeric(tb_recode),
         smokhx2 = case_when(smokhx %in% c("current","past")~"Y",
                             smokhx == "never"~"N",
                             TRUE~smokhx),
         smokhx2 = factor(smokhx2, levels=c("N","Y"))) %>% 
  arrange(sample_id)

#### Genotypes ####
# LD filtered
# MAF 1% 
GDS_file <- "~/project/result/gds/reportTB_filter_imp_LDprune.gds"
#Check if file exists
test.files <- list.files(path="~/project/result/gds/", 
                         pattern=".gds")

if(! "reportTB_filter_imp_LDprune.gds" %in% test.files){
  SNPRelate::snpgdsBED2GDS(bed.fn="~/project/result/plink/reportTB_filter_imp_LDprune.bed", 
                           fam.fn="~/project/result/plink/reportTB_filter_imp_LDprune.fam", 
                           bim.fn="~/project/result/plink/reportTB_filter_imp_LDprune.bim",
                           out.gdsfn=GDS_file)
}

genofile <- snpgdsOpen(GDS_file)

# Filter SNPs to random 10K
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
set.seed(42)
to_keep <- sample(1:length(snp.id), size = 10000)

geno.num <- snpgdsGetGeno(genofile, snp.id = snp.id[to_keep])
colnames(geno.num) <- snp.id[to_keep]
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
rownames(geno.num) <- sapply(strsplit(sample.id,"_"), `[`, 1)

#Keep samples with metadata
## Doesn't remove any
to_keep2 <- rownames(geno.num) %in% meta_FULL_PC$sample_id
geno.num.sub<-t(geno.num[to_keep2,]) %>% as.data.frame()

snpgdsClose(genofile)
rm(geno.num, geno.num.cc)
gc()

#### Kinship data ####
kin.matrix <- read_csv("~/project/result/kinship/reportTB_kin.csv") %>% 
  #keep samples with metadata
  #order by geno data
  filter(rowname %in% meta_FULL_PC$sample_id) %>%
  mutate(rowname = factor(rowname, levels=colnames(geno.num.sub))) %>% 
  arrange(rowname) %>% 
  column_to_rownames() %>% 
  select(all_of(colnames(geno.num.sub))) %>% 
  as.matrix()

check1 <- identical(rownames(kin.matrix),colnames(kin.matrix))
check2 <- identical(rownames(kin.matrix), colnames(geno.num.sub))

if(!check1 | !check2){
  stop("Kinship and genotypes do not contain the same samples in the same order.")
}

#### Run PC models ####
#Set up parallel processors
print("PC models")
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 8, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "Beta", "SE", "pval", "sigma", "AIC")
write_csv(temp, "~/project/result/model_fitting/model_pcs.csv")

#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit0 <- fit1a <- fit1b <- fit1c <- fit1d <- fit1e <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "sample_id")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit base glm ####
  fit0 <- glm(tb_recode ~ value, family = "binomial", data = geno.num.row)
  results0 <- glm_summ(fit0, snpID) %>% mutate(model = "base", .before = 0) 
  
  ####fit base glm with PCs####
  fit1a <- glm(tb_recode ~ value + PC1, family = "binomial", data = geno.num.row)
  results1a <- glm_summ(fit1a, snpID) %>% mutate(model = "PC1", .before = 0) 
  
  fit1b <- glm(tb_recode ~ value + PC1 + PC2, family = "binomial", data = geno.num.row)
  results1b <- glm_summ(fit1b, snpID) %>% mutate(model = "PC2", .before = 0) 
  
  fit1c <- glm(tb_recode ~ value + PC1 + PC2 + PC3, family = "binomial", data = geno.num.row)
  results1c <- glm_summ(fit1c, snpID) %>% mutate(model = "PC3", .before = 0) 
  
  fit1d <- glm(tb_recode ~ value + PC1 + PC2 + PC3 + PC4, family = "binomial", data = geno.num.row)
  results1d <- glm_summ(fit1d, snpID) %>% mutate(model = "PC4", .before = 0) 
  
  fit1e <- glm(tb_recode ~ value + PC1 + PC2 + PC3 + PC4 + PC5, family = "binomial", data = geno.num.row)
  results1e <- glm_summ(fit1e, snpID) %>% mutate(model = "PC5", .before = 0) 
  
  #### Save ####
  bind_rows(results0, results1a, results1b, results1c, results1d, results1e) %>% 
    write_csv(file = "~/project/result/model_fitting/model_pcs.csv", append = TRUE)
  
  return(NULL)
}

#### Run kinship model ####
print("Kinship model")
#Set up parallel processors
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 8, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "Beta", "SE", "pval", "sigma", "AIC")
write_csv(temp, "~/project/result/model_fitting/model_kin.csv")

#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit2 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "sample_id")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit kinship model####
  fit2 <- glmmkin(tb_recode ~ value + PC1 + PC2, 
                  data=geno.num.row, 
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  results2 <- glmmkin_summ(fit2, snpID) %>% mutate(model = "PC2_kinship", .before = 0) 
  
  #### Save ####
  write_csv(results2, file = "~/project/result/model_fitting/model_kin.csv", append = TRUE)
  
  return(NULL)
}

##### Run cov models ####
#320,000 lines
#Set up parallel processors
print("Cov models")
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 8, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "estimate", "std_error", "pval", "sigma", "AIC")
write_csv(temp, "~/project/result/model_fitting/model_cov.csv")

#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit3 <- fit4 <- fit5 <- fit5b <- fit6 <- fit7 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "sample_id")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit age####
  fit3 <- glmmkin(tb_recode ~ value + PC1 + PC2 + bl_age,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit3)){
    results3 <- glmmkin_summ(fit3, snpID) %>% mutate(model = "age", .before = 0)
  } else {
    results3 <- NULL
  }
  
  ####fit sex####
  fit4 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit4)){
    results4 <- glmmkin_summ(fit4, snpID) %>% mutate(model = "sex", .before = 0)
  } else {
    results4 <- NULL
  }
  
  ####fit smoke####
  fit5 <- glmmkin(tb_recode ~ value + PC1 + PC2 + smokhx,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit5)){
    results5 <- glmmkin_summ(fit5, snpID) %>% mutate(model = "smoke", .before = 0)
  } else {
    results5 <- NULL
  }
  
  fit5b <- glmmkin(tb_recode ~ value + PC1 + PC2 + smokhx2,
                   data=geno.num.row,
                   kins=kin.matrix,
                   id = "name",
                   family = binomial(link = "logit"))
  if(!is.null(fit5b)){
    results5b <- glmmkin_summ(fit5b, snpID) %>% mutate(model = "smoke2", .before = 0)
  } else {
    results5b <- NULL
  }
  
  ####fit hiv####
  fit6 <- glmmkin(tb_recode ~ value + PC1 + PC2 + bl_hiv,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit6)){
    results6 <- glmmkin_summ(fit5, snpID) %>% mutate(model = "hiv", .before = 0)
  } else {
    results6 <- NULL
  }
  
  ####fit diabetes####
  fit7 <- glmmkin(tb_recode ~ value + PC1 + PC2 + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit7)){
    results7 <- glmmkin_summ(fit7, snpID) %>% mutate(model = "diabetes", .before = 0)
  } else {
    results7 <- NULL
  }
  
  #### Save ####
  bind_rows(results3,results4,results5,results5b,results6,results7) %>%
    write_csv(file = "~/project/result/model_fitting/model_cov.csv", append = TRUE)
  
  return(NULL)
}

##### Run cov models minus 1 ####
#490,000 lines
#Set up parallel processors
print("Minus 1 models")
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 8, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "estimate", "std_error", "pval", "sigma", "AIC")
write_csv(temp, "~/project/result/model_fitting/model_cov_minus1.csv")

#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit2 <- fit3 <- fit4 <- fit5 <- fit6 <- fit7 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "sample_id")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit all cov####
  fit2 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + bl_hiv + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit2)){
    results2 <- glmmkin_summ(fit2, snpID) %>% mutate(model = "all", .before = 0)
  } else {
    results2 <- NULL
  }
  
  ####fit age####
  fit3 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + smokhx2 + bl_hiv + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit3)){
    results3 <- glmmkin_summ(fit3, snpID) %>% mutate(model = "minus_age", .before = 0)
  } else {
    results3 <- NULL
  }
  
  ####fit sex####
  fit4 <- glmmkin(tb_recode ~ value + PC1 + PC2 + bl_age + smokhx2 + bl_hiv + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit4)){
    results4 <- glmmkin_summ(fit4, snpID) %>% mutate(model = "minus_sex", .before = 0)
  } else {
    results4 <- NULL
  }
  
  ####fit smoke####
  fit5 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + bl_hiv + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit5)){
    results5 <- glmmkin_summ(fit5, snpID) %>% mutate(model = "minus_smoke", .before = 0)
  } else {
    results5 <- NULL
  }
  
  ####fit hiv####
  fit6 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit6)){
    results6 <- glmmkin_summ(fit5, snpID) %>% mutate(model = "minus_hiv", .before = 0)
  } else {
    results6 <- NULL
  }
  
  ####fit diabetes####
  fit7 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + bl_hiv,
                  data=geno.num.row,
                  kins=kin.matrix,
                  id = "name",
                  family = binomial(link = "logit"))
  if(!is.null(fit7)){
    results7 <- glmmkin_summ(fit7, snpID) %>% mutate(model = "minus_diabetes", .before = 0)
  } else {
    results7 <- NULL
  }
  
  #### Save ####
  bind_rows(results2,results3,results4,results5,results6,results7) %>%
    write_csv(file = "~/project/result/model_fitting/model_cov_minus1.csv", append = TRUE)
  
  return(NULL)
}

#### FIN ####
print("FIN")
