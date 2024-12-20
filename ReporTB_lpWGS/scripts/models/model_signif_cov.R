library(tidyverse)
library(SNPRelate)
library(GMMAT)
# library(GWASTools)
library(doParallel)
library(foreach)
# library(WGCNA)

signif <- 5E-8

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

#### Filter snp data ####
#filter to signif SNP
full.signif <- read_csv("~/project/result/model_final/model_full.csv") %>% 
  filter(term == "value" & pval < signif) %>% 
  pull(snpID)
cc.signif <- read_csv("~/project/result/model_final/model_cc.csv") %>% 
  filter(term == "value" & pval < signif) %>% 
  pull(snpID)
all.signif <- unique(c(full.signif,cc.signif))

length(full.signif)
length(cc.signif)
length(all.signif)

geno.num.sub <- geno.num[,colnames(geno.num) %in% all.signif]
dim(geno.num.sub)
geno.num.sub <- t(geno.num.sub) %>% as.data.frame()

##Split case and control
meta_case <- meta_FULL_PC %>% 
  filter(tb == "culture_confirm")
geno.num.case <- geno.num[meta_case$libID, 
                          colnames(geno.num) %in% all.signif]
geno.num.case <- t(geno.num.case) %>% as.data.frame()
kin.case <- kin.matrix[meta_case$libID, meta_case$libID]

meta_ctrl <- meta_FULL_PC %>% 
  filter(tb == "contact")
geno.num.ctrl <- geno.num[meta_ctrl$libID,
                          colnames(geno.num) %in% all.signif]
geno.num.ctrl <- t(geno.num.ctrl) %>% as.data.frame()
kin.ctrl <- kin.matrix[meta_ctrl$libID, meta_ctrl$libID]

rm(geno.num, geno.num.cc)
gc()

#### Run full model on significant SNP ####
## 79*68 = 5372 lines
#Set up parallel processors
registerDoParallel(50)

#create blank df to hold results
temp <- data.frame(matrix(ncol = 9, nrow = 0 ))
colnames(temp) <- c("model", "term", "snpID", "estimate", "std_error", "pval", "sigma", "AIC", "OR")
write_csv(temp, "~/project/result/model_final/model_signif_cov.csv")

out <- read_csv("~/project/result/model_final/model_signif_cov.csv")
complete_snp <- unique(out$snpID)

print("Run models")
#run models in a loop
foreach(i= 1:nrow(geno.num.sub)) %dopar%{
  fit1a <- fit1b <- fit2a <- fit2b <- fit3a <- fit3b <- fit4a <- fit4b <- fit5a <- fit5b <- NULL
  
  #### data ####
  #subset to row of interest
  geno.num.sub.row <- geno.num.sub[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  geno.num.row.case <- geno.num.case[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  geno.num.row.ctrl <- geno.num.ctrl[i,] %>% 
    rownames_to_column("snpID") %>% 
    pivot_longer(-snpID) %>% 
    filter(!is.na(value)) %>% 
    left_join(meta_FULL_PC, by = c("name" = "libID")) 
  snpID <- unique(geno.num.sub.row$snpID)
  snpID1 <- unique(geno.num.row.case$snpID)
  snpID2 <- unique(geno.num.row.ctrl$snpID)
  
  #### models ####
  if(!snpID %in% complete_snp){
    #### sex ####
    fit1a <- glmmkin(sex ~ value + PC1 + PC2 + bl_age + smokhx2 + bl_hiv + diabetes,
                     data=geno.num.row.case,
                     kins=kin.case,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit1a)){
      results1a <- glmmkin_summ(fit1a, snpID1) %>% 
        mutate(model = "sex_case", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "sex_case",
                     term = names(exp(coef(fit1a))),
                     snpID = snpID,
                     OR = exp(coef(fit1a))),
          by = join_by(model, term, snpID))
    } else {
      results1a <- NULL
    }
    
    fit1b <- glmmkin(sex ~ value + PC1 + PC2 + bl_age + smokhx2 + bl_hiv + diabetes,
                     data=geno.num.row.ctrl,
                     kins=kin.ctrl,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit1b)){
      results1b <- glmmkin_summ(fit1b, snpID2) %>% 
        mutate(model = "sex_ctrl", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "sex_ctrl",
                     term = names(exp(coef(fit1b))),
                     snpID = snpID,
                     OR = exp(coef(fit1b))),
          by = join_by(model, term, snpID))
    } else {
      results1b <- NULL
    }
    
    #### hiv ####
    fit2a <- glmmkin(bl_hiv ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes,
                     data=geno.num.row.case,
                     kins=kin.case,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit2a)){
      results2a <- glmmkin_summ(fit2a, snpID1) %>% 
        mutate(model = "hiv_case", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "hiv_case",
                     term = names(exp(coef(fit2a))),
                     snpID = snpID,
                     OR = exp(coef(fit2a))),
          by = join_by(model, term, snpID))
    } else {
      results2a <- NULL
    }
    
    fit2b <- glmmkin(bl_hiv ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes,
                     data=geno.num.row.ctrl,
                     kins=kin.ctrl,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit2b)){
      results2b <- glmmkin_summ(fit2b, snpID2) %>% 
        mutate(model = "hiv_ctrl", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "hiv_ctrl",
                     term = names(exp(coef(fit2b))),
                     snpID = snpID,
                     OR = exp(coef(fit2b))),
          by = join_by(model, term, snpID))
    } else {
      results2b <- NULL
    }
    
    #### diabetes ####
    fit3a <- glmmkin(diabetes ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + bl_hiv,
                     data=geno.num.row.case,
                     kins=kin.case,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit3a)){
      results3a <- glmmkin_summ(fit3a, snpID1) %>% 
        mutate(model = "diabetes_case", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "diabetes_case",
                     term = names(exp(coef(fit3a))),
                     snpID = snpID,
                     OR = exp(coef(fit3a))),
          by = join_by(model, term, snpID))
    } else {
      results3a <- NULL
    }
    
    fit3b <- glmmkin(diabetes ~ value + PC1 + PC2 + sex + bl_age + smokhx2 + bl_hiv,
                     data=geno.num.row.ctrl,
                     kins=kin.ctrl,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit3b)){
      results3b <- glmmkin_summ(fit3b, snpID2) %>% 
        mutate(model = "diabetes_ctrl", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "diabetes_ctrl",
                     term = names(exp(coef(fit3b))),
                     snpID = snpID,
                     OR = exp(coef(fit3b))),
          by = join_by(model, term, snpID))
    } else {
      results3b <- NULL
    }
    
    #### age ####
    fit4a <- glmmkin(bl_age ~ value + PC1 + PC2 + sex + smokhx2 + bl_hiv + diabetes,
                     data=geno.num.row.case,
                     kins=kin.case,
                     id = "name",
                     family = gaussian(link = "identity"))
    if(!is.null(fit4a)){
      results4a <- glmmkin_summ(fit4a, snpID1) %>% 
        mutate(model = "age_case", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "age_case",
                     term = names(exp(coef(fit4a))),
                     snpID = snpID,
                     OR = exp(coef(fit4a))),
          by = join_by(model, term, snpID))
    } else {
      results4a <- NULL
    }
    
    fit4b <- glmmkin(bl_age ~ value + PC1 + PC2 + sex + smokhx2 + bl_hiv + diabetes,
                     data=geno.num.row.ctrl,
                     kins=kin.ctrl,
                     id = "name",
                     family = gaussian(link = "identity"))
    if(!is.null(fit4b)){
      results4b <- glmmkin_summ(fit4b, snpID2) %>% 
        mutate(model = "age_ctrl", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "age_ctrl",
                     term = names(exp(coef(fit4b))),
                     snpID = snpID,
                     OR = exp(coef(fit4b))),
          by = join_by(model, term, snpID))
    } else {
      results4b <- NULL
    }
    
    #### smoke ####
    fit5a <- glmmkin(smokhx2 ~ value + PC1 + PC2 + sex + bl_age + bl_hiv + diabetes,
                     data=geno.num.row.case,
                     kins=kin.case,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit5a)){
      results5a <- glmmkin_summ(fit5a, snpID1) %>% 
        mutate(model = "smoke_case", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "smoke_case",
                     term = names(exp(coef(fit5a))),
                     snpID = snpID,
                     OR = exp(coef(fit5a))),
          by = join_by(model, term, snpID))
    } else {
      results5a <- NULL
    }
    
    fit5b <- glmmkin(smokhx2 ~ value + PC1 + PC2 + sex + bl_age + bl_hiv + diabetes,
                     data=geno.num.row.ctrl,
                     kins=kin.ctrl,
                     id = "name",
                     family = binomial(link = "logit"))
    if(!is.null(fit5b)){
      results5b <- glmmkin_summ(fit5b, snpID2) %>% 
        mutate(model = "smoke_ctrl", .before = 0) %>% 
        full_join(
          #odds-ratio
          data.frame(model = "smoke_ctrl",
                     term = names(exp(coef(fit5b))),
                     snpID = snpID,
                     OR = exp(coef(fit5b))),
          by = join_by(model, term, snpID))
    } else {
      results5b <- NULL
    }
    
    #### Save ####
    bind_rows(results1a, results1b, results2a, results2b,
              results3a, results3b, results4a, results4b,
              results5a, results5b) %>%
      write_csv(file = "~/project/result/model_final/model_signif_cov.csv", 
                append = TRUE)
  }
  
  return(NULL)
}

print("Fin")
#### FIN ####
