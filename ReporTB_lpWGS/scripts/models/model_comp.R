library(tidyverse)
library(SNPRelate)
library(GMMAT)
# library(GWASTools)
library(doParallel)
library(foreach)
# library(WGCNA)

#### Fxns ####
print("Load data")
source("~/project/scripts/models/glmmkin_summ.R")

#### Data #### 
print("Load data")
MAF <- 5
source("~/project/scripts/models/model_load_data.R")

#Check data
count(meta_FULL_PC, tb)
dim(geno.num)
count(meta_CC_PC, tb)
dim(geno.num.cc)

#### Run score tests ####
print("Score tests - FULL")
glmmkin_full <- glmmkin(tb_recode ~ PC1 + PC2 + sex + bl_age,
                        data = meta_FULL_PC,
                        kins = kin.matrix,
                        id = "libID",
                        family = binomial(link = "logit"))

glmm.score(glmmkin_full,
           infile = "/home/ec2-user/project/result/plink/reportTB_filter_imp_maf5",
           infile.header.print = c("snpID"),
           outfile = "/home/ec2-user/project/result/model_final/glmm_score_full_comp.txt",
           MAF.range = c(0.01, 0.99))

#### Run full model on significant scores ####
# Get sig snps
sig_scores_comp <- read_tsv("~/project/result/model_final/glmm_score_full_comp.txt")%>%
  filter(PVAL < 1e-4)
# Add orig model SNPs
sig_scores <- read_tsv("~/project/result/model_final/glmm_score_full.txt") %>%
  filter(PVAL < 1e-4)

sig_scores_all <- unique(c(sig_scores_comp$SNP,sig_scores$SNP))

# Filter genotype data
geno.num.sub <- t(geno.num[,colnames(geno.num) %in% sig_scores_all]) %>% 
  as.data.frame()
rm(geno.num)
gc()

#Set up parallel processors
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 9, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "estimate", "std_error", "pval", "sigma", "AIC", "OR")
write_csv(temp, "~/project/result/model_final/model_full_comp.csv")

out <- read_csv("~/project/result/model_final/model_full_comp.csv")
complete_snp <- unique(out$snpID)

print("Run models - FULL")
#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit4 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit sex####
  if(!snpID %in% complete_snp){
    fit4 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit4)){
      results4 <- glmmkin_summ(fit4, snpID) %>% mutate(model = "full_comp",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "full_comp",
                     term = names(exp(coef(fit4))),
                     snpID = snpID,
                     OR = exp(coef(fit4))),
          by = join_by(model, term, snpID))
      
    } else {
      results4 <- NULL
    }
    
    #### Save ####
    results4 %>%
      write_csv(file = "~/project/result/model_final/model_full_comp.csv", append = TRUE)
  }
  
  return(NULL)
}

#### Case-control subset ####
#### Run CC score tests ####
print("Score tests - CC")
glmmkin_cc <- glmmkin(tb_recode ~ PC1 + PC2 + sex + bl_age,
                      data = meta_CC_PC,
                      kins = kin.matrix_CC,
                      id = "libID",
                      family = binomial(link = "logit"))

glmm.score(glmmkin_cc,
           infile = "/home/ec2-user/project/result/plink/reportTB_filter_imp_maf5_cc",
           infile.header.print = c("snpID"),
           outfile = "/home/ec2-user/project/result/model_final/glmm_score_cc_comp.txt",
           MAF.range = c(0.01, 0.99))

#### Run CC model on significant scores ####
sig_scores_comp_cc <- read_tsv("~/project/result/model_final/glmm_score_cc_comp.txt") %>%
  filter(PVAL < 1e-4)
# Add orig model SNPs
sig_scores_cc <- read_tsv("~/project/result/model_final/glmm_score_cc.txt") %>%
  filter(PVAL < 1e-4)
#Add full cohort SNPs
## Already loaded above
# sig_scores_comp <- read_tsv("~/project/result/model_final/glmm_score_full_comp.txt") %>%
#   filter(PVAL < 1e-4)
# sig_scores <- read_tsv("~/project/result/model_final/glmm_score_full.txt") %>%
#   filter(PVAL < 1e-4)

sig_scores_all2 <- unique(c(sig_scores_comp_cc$SNP,sig_scores_cc$SNP,
                           sig_scores_comp$SNP,sig_scores$SNP))

# Filter genotype data
geno.num.cc.sub <- t(geno.num.cc[,colnames(geno.num.cc) %in% sig_scores_all2]) %>% 
  as.data.frame()
rm(geno.num.cc)
gc()

#create blank df to hold results
write_csv(temp, "~/project/result/model_final/model_cc_comp.csv")

out2 <- read_csv("~/project/result/model_final/model_cc_comp.csv")
complete_snp2 <- unique(out2$snpID)

print("Run models - CC")
#run models in a loop
foreach(i= 1:nrow(geno.num.cc.sub)) %dopar%{
  fit4 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.cc.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_CC_PC, by = c("name" = "libID")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####fit sex####
  if(!snpID %in% complete_snp2){
    fit4 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age,
                    data=geno.num.row,
                    kins=kin.matrix_CC,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit4)){
      results4 <- glmmkin_summ(fit4, snpID) %>% mutate(model = "cc_comp",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "cc_comp",
                     term = names(exp(coef(fit4))),
                     snpID = snpID,
                     OR = exp(coef(fit4))),
          by = join_by(model, term, snpID))
      
    } else {
      results4 <- NULL
    }
    
    #### Save ####
    results4 %>%
      write_csv(file = "~/project/result/model_final/model_cc_comp.csv", append = TRUE)
  }
  
  return(NULL)
}
print("Fin")

#### FIN ####
