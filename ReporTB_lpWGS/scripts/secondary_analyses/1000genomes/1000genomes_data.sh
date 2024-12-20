# PCA of 1000 genomes + ReporTB
# Adapted from Laura Budurlean
# https://github.com/laura-budurlean/PCA-Ethnicity-Determination-from-WGS-Data

setwd result/1000genomes/

# 1. Filter ReporTB
mkdir -p 1_raw/
## Remove AT GC SNPs
## See https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html#ancestry-estimation-1
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    ../plink_concord/reportTB_filter_imp_maf5_concord.bim  > \
    1_raw/report.ac_gt_snps
    
/Applications/plink2 --bfile ../plink_concord/reportTB_filter_imp_maf5_concord \
      --exclude 1_raw/report.ac_gt_snps \
      --make-bed \
      --out 1_raw/reportTB_filter_imp_maf5_concord_filter

# 2. Convert ReporTB to vcf
/Applications/plink2 --bfile 1_raw/reportTB_filter_imp_maf5_concord_filter \
  --recode vcf \
  --out 1_raw/reportTB_filter_imp_maf5_concord
bgzip 1_raw/reportTB_filter_imp_maf5_concord.vcf
bcftools index 1_raw/reportTB_filter_imp_maf5_concord.vcf.gz

# 3. Create "site2.bed" file which is SNP positions present in your samples, 
# that you will then use to filter out variants in the 1000 genomes samples. 
awk -F '\t' '{print $1"\t"$4-1"\t"$4}' 1_raw/reportTB_filter_imp_maf5_concord_filter.bim > site2.bed

# 4. Download 1000 genomes vcf data
prefix="ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr"
suffix=".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

for chr in {1..22}; do
  echo $chr
  curl "${prefix}""${chr}""${suffix}" --output 1_raw/chr"${chr}".vcf.gz
  curl "${prefix}""${chr}""${suffix}".tbi --output 1_raw/chr"${chr}".vcf.gz.tbi
done

# 5. Filter 1000 genome against ReporTB sites
mkdir -p 2_filter

for chr in {1..22}; do
  echo $chr
  bcftools filter -T site2.bed -Oz \
    -o 2_filter/chr${chr}.1kg.sites.vcf.gz \
    1_raw/chr${chr}.vcf.gz 
done

# 6. Combine all the chromosome files 
mkdir -p 3_merged
bcftools concat 2_filter/chr*.1kg.sites.vcf.gz > 3_merged/merged.1kg.sites.vcf
bgzip 3_merged/merged.1kg.sites.vcf
bcftools index 3_merged/merged.1kg.sites.vcf.gz

# 7. Intersect the reporTB.vcf.gz and merged.1kg.sites.vcf.gz
bcftools isec --threads 4 -p isec -n =2 1_raw/reportTB_filter_imp_maf5_concord.vcf.gz 3_merged/merged.1kg.sites.vcf.gz

bgzip --threads 4 isec/0000.vcf
bgzip --threads 4 isec/0001.vcf 
bcftools index isec/0000.vcf.gz
bcftools index isec/0001.vcf.gz
bcftools merge -m none isec/0000.vcf.gz isec/0001.vcf.gz \
  -Oz -o 3_merged/merged.1kg.reportb.vcf.gz

# 8. Set missing genotypes to 0|0
mkdir -p 4_merge_filter

bcftools +setGT 3_merged/merged.1kg.reportb.vcf.gz -- -t . -n 0p\n \
  > 4_merge_filter/merged.1kg.reportb.setGT.vcf.gz
## Filled 0 alleles

/Applications/plink2 --vcf 4_merge_filter/merged.1kg.reportb.setGT.vcf.gz \
  --set-missing-var-ids @:#[b38]\$r\$a --make-bed \
  --out 4_merge_filter/merged.1kg.reportb.setGT
  
/Applications/plink2 --bfile 4_merge_filter/merged.1kg.reportb.setGT \
  --indep-pairwise 100kb 1 .1

/Applications/plink2 --bfile 4_merge_filter/merged.1kg.reportb.setGT \
  --extract plink2.prune.in \
  --maf 0.1 --pca biallelic-var-wts --make-bed \
  --out 4_merge_filter/merged.1kg.reportb.setGT.prune
 
 wc -l 4_merge_filter/merged.1kg.reportb.setGT.prune.bim
 
# 9. PCA scoring
mkdir -p 5_pca

/Applications/plink2 --bfile 4_merge_filter/merged.1kg.reportb.setGT \
  --score 4_merge_filter/merged.1kg.reportb.setGT.prune.eigenvec.var 2 3 5 \
  --out 5_pca/merged.1kg.reportb.PC1
  
/Applications/plink2 --bfile 4_merge_filter/merged.1kg.reportb.setGT \
  --score 4_merge_filter/merged.1kg.reportb.setGT.prune.eigenvec.var 2 3 6 \
  --out 5_pca/merged.1kg.reportb.PC2

/Applications/plink2 --bfile 4_merge_filter/merged.1kg.reportb.setGT \
  --score 4_merge_filter/merged.1kg.reportb.setGT.prune.eigenvec.var 2 3 7 \
  --out 5_pca/merged.1kg.reportb.PC3
  
# 10. plot
Rscript ../../scripts/secondary_analyses/1000genomes/1000genomes_pca.R
