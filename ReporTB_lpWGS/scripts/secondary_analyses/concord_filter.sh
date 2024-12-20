Rscript scripts/secondary_analyses/concord_filter.R

#Filter concordant SNPs 
/Applications/plink2 \
  --bfile result/plink/reportTB_filter_imp_maf5 \
    --make-bed --out result/plink_concord/reportTB_filter_imp_maf5_concord \
    --extract result/snp_imp_maf5_concord.txt \
    --threads 6
## 2527008 SNPs for analysis
