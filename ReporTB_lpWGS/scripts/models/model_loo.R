# Leave-one-out models for all SNP: HIV, diabetes, smoking
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

#### Pull full model significant ####
#MAF5 significant only
fdr <- read_csv("~/project/result/model_final/model_full.csv") %>% 
  filter(term=="value" & pval < 5E-8 & snpID %in% colnames(geno.num))

#### Leave-one-out models ####
## 24*68 = 1632 lines
# Filter genotype data
geno.num.sub <- t(geno.num[,colnames(geno.num) %in% fdr$snpID]) %>% 
  as.data.frame()
rm(geno.num, geno.num.cc)
gc()

#Set up parallel processors
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 9, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "estimate", "std_error", "pval", "sigma", "AIC", "OR")
write_csv(temp, "~/project/result/model_final/model_loo.csv")

out <- read_csv("~/project/result/model_final/model_loo.csv")
complete_snp <- unique(out$snpID)

print("Run models - LOO")
#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit1 <- fit2 <- fit3 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####Minus HIV####
  if(!snpID %in% complete_snp){
    fit1 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit1)){
      results1 <- glmmkin_summ(fit1, snpID) %>% mutate(model = "minus_hiv",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "minus_hiv",
                     term = names(exp(coef(fit1))),
                     snpID = snpID,
                     OR = exp(coef(fit1))),
          by = join_by(model, term, snpID))
      
    } else {
      results1 <- NULL
    }
    
    ##### Minus diabetes ####
    fit2 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + bl_hiv,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit2)){
      results2 <- glmmkin_summ(fit2, snpID) %>% mutate(model = "minus_diabetes",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "minus_diabetes",
                     term = names(exp(coef(fit2))),
                     snpID = snpID,
                     OR = exp(coef(fit2))),
          by = join_by(model, term, snpID))
      
    } else {
      results2 <- NULL
    }
    
    #### Minus smoke ####
    fit3 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + bl_hiv + diabetes,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit3)){
      results3 <- glmmkin_summ(fit3, snpID) %>% mutate(model = "minus_smoke",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "minus_smoke",
                     term = names(exp(coef(fit3))),
                     snpID = snpID,
                     OR = exp(coef(fit3))),
          by = join_by(model, term, snpID))
      
    } else {
      results3 <- NULL
    }
    #### Save ####
    bind_rows(results1,results2,results3) %>%
      write_csv(file = "~/project/result/model_final/model_loo.csv", append = TRUE)
  }
  
  return(NULL)
}

#### Univariate models ####
## 21*68 = 1428 lines
#create blank df to hold results
write_csv(temp, "~/project/result/model_final/model_univar.csv")

out <- read_csv("~/project/result/model_final/model_univar.csv")
complete_snp <- unique(out$snpID)

print("Run models - UNIVAR")
#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit1 <- fit2 <- fit3 <- NULL
  #subset to row of interest
  geno.num.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  snpID <- unique(geno.num.row$snpID)
  
  ####Add HIV####
  if(!snpID %in% complete_snp){
    fit1 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + bl_hiv,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit1)){
      results1 <- glmmkin_summ(fit1, snpID) %>% mutate(model = "plus_hiv",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "plus_hiv",
                     term = names(exp(coef(fit1))),
                     snpID = snpID,
                     OR = exp(coef(fit1))),
          by = join_by(model, term, snpID))
      
    } else {
      results1 <- NULL
    }
    
    ##### Add diabetes ####
    fit2 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + diabetes,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit2)){
      results2 <- glmmkin_summ(fit2, snpID) %>% mutate(model = "plus_diabetes",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "plus_diabetes",
                     term = names(exp(coef(fit2))),
                     snpID = snpID,
                     OR = exp(coef(fit2))),
          by = join_by(model, term, snpID))
      
    } else {
      results2 <- NULL
    }
    
    #### Add smoke ####
    fit3 <- glmmkin(tb_recode ~ value + PC1 + PC2 + sex + bl_age + smokhx2,
                    data=geno.num.row,
                    kins=kin.matrix,
                    id = "name",
                    family = binomial(link = "logit"))
    if(!is.null(fit3)){
      results3 <- glmmkin_summ(fit3, snpID) %>% mutate(model = "plus_smoke",
                                                       .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "plus_smoke",
                     term = names(exp(coef(fit3))),
                     snpID = snpID,
                     OR = exp(coef(fit3))),
          by = join_by(model, term, snpID))
      
    } else {
      results3 <- NULL
    }
    #### Save ####
    bind_rows(results1,results2,results3) %>%
      write_csv(file = "~/project/result/model_final/model_univar.csv", append = TRUE)
  }
  
  return(NULL)
}

#### FIN ####
print("FIN")
