# Create files
Rscript scripts/secondary_analyses/magma/create_magma.R

cd /Users/kadm/Applications/magma_v1.10_mac

# Map to genes
./magma \
  --annotate window=0 \
  --snp-loc /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/plink_concord/reportTB_filter_imp_maf5_concord.bim \
  --gene-loc /Users/kadm/Applications/magma_v1.10_mac/ref/NCBI38/NCBI38.gene.loc \
  --out /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReporTB_magma

# Gene level analysis
for i in {1..22}
do
  echo $i
  
  ./magma --bfile /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/plink_concord/reportTB_filter_imp_maf5_concord \
    --pval /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReportTB_magma_pval.txt \
    ncol=NOBS \
    --gene-annot /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReporTB_magma.genes.annot \
    --gene-model multi \
    --batch $i chr \
    --out /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/batch/ReporTB_magma_gene
done

## Combine
./magma --merge /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/batch/ReporTB_magma_gene \
    --out /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReporTB_magma_gene

# Gene set analysis
./magma --gene-results /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReporTB_magma_gene.genes.raw \
  --set-annot /Users/kadm/Applications/magma_v1.10_mac/ref/MSigDB/msigdb_7.5.1_H.C2.C5BP.txt \
  col=2,1 \
  --out /Users/kadm/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/result/magma/ReporTB_magma_geneset

# FDR correction
cd ~/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS/
Rscript scripts/secondary_analyses/magma/magma_fdr.R
