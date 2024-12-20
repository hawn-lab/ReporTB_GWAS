#GCTA, REML method
#See https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpz1.734
#From indiv level data

#Prepare metadata
## See scripts/secondary_analysis/heritability_covar.R

cd result

#Calculate GRM
## 3 hrs
/Applications/gcta-1.94.2/gcta64 \
  --bfile plink_concord/reportTB_filter_imp_maf5_concord \
  --autosome --make-grm-bin \
  --out heritability/reportTB_filter_imp_maf5_concord \
  --thread-num 10

# 2754 individuals 
# 2527008 SNPs from chromosome 1 to chromosome 22
# 
# Summary of the GRM:
# Mean of diagonals = 1.01767
# Variance of diagonals = 0.0266436
# Mean of off-diagonals = -0.000369676
# Variance of off-diagonals = 0.00128047

#Remove 1 from related pairs
/Applications/gcta-1.94.2/gcta64 \
  --grm heritability/reportTB_filter_imp_maf5_concord \
  --grm-cutoff 0.05 --make-grm-bin \
  --out heritability/reportTB_filter_imp_maf5_concord_0.05 \
  --thread-num 10

# 1025 individuals
head heritability/reportTB_filter_imp_maf5_concord_0.05.grm.id

#Model heritaibility 
##Correct for all covar used in GWAS
/Applications/gcta-1.94.2/gcta64 \
  --grm heritability/reportTB_filter_imp_maf5_concord_0.05 \
  --pheno heritability/TB.phen \
  --reml \
  --thread-num 10 \
  --covar heritability/covar_cat.txt \
  --qcovar heritability/covar_cont.txt \
  --prevalence 0.0005707 \
  --out heritability/reportTB_filter_imp_maf5_concord_greml 
  
#No clinical or risk cov. PCs only
/Applications/gcta-1.94.2/gcta64 \
  --grm heritability/reportTB_filter_imp_maf5_concord_0.05 \
  --pheno heritability/TB.phen \
  --reml \
  --thread-num 10 \
  --qcovar heritability/covar_PC.txt \
  --prevalence 0.0005707 \
  --out heritability/reportTB_filter_imp_maf5_concord_greml_PConly
 
#No addtl TB risk cov 
/Applications/gcta-1.94.2/gcta64 \
  --grm heritability/reportTB_filter_imp_maf5_concord_0.05 \
  --pheno heritability/TB.phen \
  --reml \
  --thread-num 10 \
  --covar heritability/covar_cat_noTB.txt \
  --qcovar heritability/covar_cont.txt \
  --prevalence 0.0005707 \
  --out heritability/reportTB_filter_imp_maf5_concord_greml_noTB

### Results ###
# Variance Components
# V(G): genetic variance. proportion of phenotypic variance due to additive genetic factors.
# V(e): environmental variance, proportion of phenotypic variance due to environmental factors.
# Vp: total phenotypic variance, genetic + environmental variance
# 
# V(G)/Vp: heritability estimate, proportion of phenotypic variance attributable to genetic variance
# 
# Adjusted Heritability Estimate
# V(G)/Vp_L: heritability estimate transformed to the underlying scale
# 
# logL: log-likelihood of the model
# logL0: log-likelihood of the null model
# LRT (Likelihood Ratio Test): difference between the two models, testing whether including the genetic variance significantly improves the model fit.

# All covar
## 3 quantitative covariate(s)
## 4 discrete covariate(s)
## 426 cases, 594 controls (1020 total)
## 
## Source  Variance        SE
## V(G)    0.096425        0.049312
## V(e)    0.105845        0.047190
## Vp      0.202269        0.009250
## V(G)/Vp 0.476714        0.237754
## Scaled
## V(G)/Vp_L       0.158459        0.079029

# PCs only
## 2 quantitative covariate(s)
## 2 quantitative variable(s)
## 430 cases, 595 controls (1025 total)
## 
## Source  Variance        SE
## V(G)    0.030683        0.053477
## V(e)    0.210270        0.052587
## Vp      0.240953        0.010794
## V(G)/Vp 0.127338        0.220969
## Scaled
## V(G)/Vp_L       0.042274        0.073358

# No TB risk
## 3 quantitative variable(s)
## 1 discrete variable(s)
## 430 cases, 595 controls (1025 total)
## 
## Source  Variance        SE
## V(G)    0.040896        0.050955
## V(e)    0.184112        0.049871
## Vp      0.225008        0.010111
## V(G)/Vp 0.181754        0.224923
## Scaled
## V(G)/Vp_L       0.060339        0.074671
