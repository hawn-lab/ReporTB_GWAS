#### Directory setup ####
mkdir -p result/plink
mkdir -p result/gds
mkdir -p result/kinship
mkdir -p result/logs
## Make temp directory 
# mkdir -p ~/project/tmp

#### Convert to PLINK ####
# Use PLINK file created in original case-control analyses

#### Filter patients ####
# Keep patients with single cell data
/Applications/plink2 -bfile ../ReporTB_lpWGS/result/plink/reportTB \
        --make-bed --out result/plink/reportTB_sc_filter \
        --keep result/sc_patients.txt \
        --threads 4

#### Filter imputation score > 0.6 ####
/Applications/plink2 -bfile result/plink/reportTB_sc_filter \
        --make-bed --out result/plink/reportTB_sc_filter_imp \
        --extract ../ReporTB_lpWGS/result/snp_impute_pf06.txt \
        --threads 4

#### LD and MAF filter for ancestry analysis ####
/Applications/plink2 -bfile result/plink/reportTB_sc_filter_imp \
        --make-bed --out result/plink/reportTB_sc_filter_imp_LDprune \
        --maf 0.01 \
        --indep-pairwise 500kb 0.2 \
        --threads 4

#### Calc IBD and kinship ####
nohup Rscript scripts/snp/ancestry_PCs.R >> result/logs/ancestry_PCs.log 2>&1 &

#### HWE and MAF filters ####
# MAF 1% or 5%
# HWE P < 1E-6 in controls only

# 58 samples loaded
# 0 cases, 33 controls
# 31088846 variants loaded 

/Applications/plink2 --bfile result/plink/reportTB_sc_filter_imp \
    --make-bed --out result/plink/reportTB_sc_filter_imp_maf1 \
    --maf 0.01 --hwe 1e-6 \
    --threads 4
# 1699 variants removed due to hwe
# 18695039 variants removed due to maf
# 12392108 variants remaining

/Applications/plink2 --bfile result/plink/reportTB_sc_filter_imp \
    --make-bed --out result/plink/reportTB_sc_filter_imp_maf5 \
    --maf 0.05 --hwe 1e-6 \
    --threads 4
# 1699 variants removed due to hwe
# 23301123 variants removed due to maf
# 7786024 variants remaining after

echo "FIN"
