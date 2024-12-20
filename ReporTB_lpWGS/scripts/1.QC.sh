#### Directory setup ####
mkdir -p ~/project/result
mkdir -p ~/project/result/logs
mkdir -p ~/project/result/plink/
mkdir -p ~/project/result/gds/
mkdir -p ~/project/result/other/
mkdir -p ~/project/scripts
## Make temp directory 
mkdir -p ~/project/tmp

#### Convert to PLINK ####
~/apps/plink2 --allow-extra-chr \
        --make-bed --out ~/project/result/plink/reportTB \
        --vcf ~/project/data/7212-BA-imputed-merged.vcf.gz \
        --set-all-var-ids @:#:\$r:\$a \
        --chr 1-26, X, Y, XY --split-par hg38 \
        --snps-only --max-alleles 2 \
        --threads 50

# Save original fam file
cp ~/project/result/plink/reportTB.fam ~/project/result/plink/reportTB_orig.fam

# Add case-control info to fam file
nohup Rscript ~/project/scripts/QC/case_control_classification.R >> ~/project/result/logs/case_control_classification.log 2>&1 &

#### Filter patients ####
# Remove patients not in metadata
nohup Rscript ~/project/scripts/QC/pt_remove.R >> ~/project/result/logs/pt_remove.log 2>&1 &

~/apps/plink2 -bfile ~/project/result/plink/reportTB \
        --make-bed --out ~/project/result/plink/reportTB_filter \
        --remove ~/project/result/other/pt_remove.txt \
        --threads 50

#### Filter imputation score > 0.6 ####
nohup bash ~/project/scripts/qual_from_vcf.sh >> ~/project/result/logs/qual_vcf.log 2>&1 &

~/apps/plink2 -bfile ~/project/result/plink/reportTB_filter \
        --make-bed --out ~/project/result/plink/reportTB_filter_imp \
        --extract ~/project/result/snp_impute_pf06.txt \
        --threads 50

#### LD and MAF filter for ancestry analysis ####
~/apps/plink2 -bfile ~/project/result/plink/reportTB_filter_imp \
        --make-bed --out ~/project/result/plink/reportTB_filter_imp_LDprune \
        --maf 0.01 \
        --indep-pairwise 500kb 0.2 \
        --threads 50

#### Calc IBD and kinship ####
nohup Rscript ~/project/scripts/QC/ancestry_PCs.R >> ~/project/result/logs/ancestry_PCs.log 2>&1 &

#### HWE and MAF filters ####
# MAF 1% or 5%
# HWE P < 1E-6 in controls only
# Remove duplicates and unconfirmed twins

# 2763 samples loaded
# 952 cases, 1811 controls
# 31088846 variants loaded 

~/apps/plink2 --bfile ~/project/result/plink/reportTB_filter_imp \
    --make-bed --out ~/project/result/plink/reportTB_filter_imp_maf1 \
    --maf 0.01 --hwe 1e-6 \
    --remove ~/project/result/other/dups_remove.txt \
    --threads 50
# 192373 variants removed due to hwe
# 18503146 variants removed due to maf
# 12393327 variants remaining

~/apps/plink2 --bfile ~/project/result/plink/reportTB_filter_imp \
    --make-bed --out ~/project/result/plink/reportTB_filter_imp_maf5 \
    --maf 0.05 --hwe 1e-6 \
    --remove ~/project/result/other/dups_remove.txt \
    --threads 50
# 192373 variants removed due to hwe
# 23952020 variants removed due to maf
# 6944453 variants remaining after

# 2754 samples remaining
# 947 cases and 1807 controls remaining 

#### Filter case-control subset ####
# Needed for glmm score
~/apps/plink2 --bfile ~/project/result/plink/reportTB_filter_imp_maf1 \
    --make-bed --out ~/project/result/plink/reportTB_filter_imp_maf1_cc \
    --keep ~/project/result/other/cc_keep.txt \
    --threads 50

~/apps/plink2 --bfile ~/project/result/plink/reportTB_filter_imp_maf5 \
    --make-bed --out ~/project/result/plink/reportTB_filter_imp_maf5_cc \
    --keep ~/project/result/other/cc_keep.txt \
    --threads 50

#### Save ####
aws s3 sync ~/project/result/ s3://hawn-reporttb-results2/
