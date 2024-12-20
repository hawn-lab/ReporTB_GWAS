# 1. Use SEAsnake for initial data cleaning
## https://bigslu.github.io/SEAsnake/vignette/SEAsnake_vignette.html
## Version 1.1 which includes bam indexing

cores=60
threads=60
nohup snakemake --snakefile Snakefile_step1 --cores $cores >> log/SEAsnake_step1.log 2>&1 &

## trim5p: 10
## adapter sequences: Illumina Universal adapters
## genome: 'Homo_sapiens.GRCh38'
## release: '112'

## Skip Flagstat, Picard, Subread, and combine steps
nohup snakemake --snakefile Snakefile_step2 --cores $cores \
    --allowed-rules adapterremoval fastqc_trim STAR_index STAR_load STAR_align STAR_remove align_filter \
    >> log/SEAsnake_step2.log 2>&1 &

# aws s3 sync ~/SEAsnake/result/3_bam_filter/ s3://hawn-reporttb-wgs2/3_bam_filter/ 

# 2. Call SNPs
## Index human genome
mkdir -p ~/SEAsnake/ref2/STARref
cp ~/SEAsnake/ref/release112/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa ~/SEAsnake/ref2/STARref/
cd ~/SEAsnake/ref2/STARref/
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

## Then get SNP from filtered BAM
mkdir -p ~/SEAsnake/result/4_vcf/
touch ~/SEAsnake/log/mpileup.log

for bam in ~/SEAsnake/result/3_bam_filter/*.bam;
do
  NAME=$(basename "$bam" _S1_Aligned.sortedByCoord.filter.bam)
  echo $NAME >> ~/SEAsnake/log/mpileup.log
  parallel -j 22 'bcftools mpileup --fasta-ref ~/SEAsnake/ref2/STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --min-BQ 20 --min-MQ 30 --skip-indels --threads {3} -r {1} {2} | bcftools call --output-type v --consensus-caller --threads {3} --ploidy-file ~/SEAsnake/ploidy.txt | bcftools view --output {2}_chr{1}.vcf.gz --output-type z --threads {3}' ::: {1..22} ::: $bam ::: 10
  mv ~/SEAsnake/result/3_bam_filter/*.vcf.gz ~/SEAsnake/result/4_vcf/
done

## Fix names
cd ~/SEAsnake/result/4_vcf/
rename '_S1_Aligned.sortedByCoord.filter.bam' '' *

# aws s3 sync ~/SEAsnake/result/4_vcf/ s3://hawn-reporttb-wgs2/4_vcf/ 

# 3. Index vcf
touch ~/SEAsnake/log/vcf_index.log
cd ~/SEAsnake/result/4_vcf/

for bam in ~/SEAsnake/result/3_bam_filter/*.bam;
do
  NAME=$(basename "$bam" .bam)
  echo $NAME >> ~/SEAsnake/log/vcf_index.log
  parallel -j 22 'bcftools index -t {2}_chr{1}.vcf.gz' ::: {1..22} ::: $NAME
done

# aws s3 sync ~/SEAsnake/result/4_vcf/ s3://hawn-reporttb-wgs2/4_vcf/ 

# 4. Merge files by chr
mkdir -p ~/SEAsnake/result/5_combined/
parallel -j 11 'bcftools merge ~/SEAsnake/result/4_vcf/*_chr{1}.vcf.gz -Oz -o ~/SEAsnake/result/5_combined/reportb_wgs_chr{1}.vcf.gz' ::: {1..22}  

# aws s3 sync ~/SEAsnake/result/5_combined/ s3://hawn-reporttb-wgs2/5_combined/ 

# 5. Convert to PLINK
mkdir -p ~/SEAsnake/result/6_plink/

for vcf in ~/SEAsnake/result/5_combined/*vcf.gz;
do
  NAME=$(basename "$vcf" .vcf.gz)
  echo $NAME

 plink2 --allow-extra-chr \
        --make-bed --out ~/SEAsnake/result/6_plink/"$NAME" \
        --vcf $vcf \
        --set-all-var-ids @:#:\$r:\$a \
        --chr 1-22 --split-par hg38 \
        --snps-only --max-alleles 2 \
        --threads 60
done

# aws s3 sync ~/SEAsnake/result/6_plink/ s3://hawn-reporttb-wgs2/6_plink/

# 6. Filter SNP present in low-pass data, no MAF filter
awk '{print $1, $4, $4}' ~/SEAsnake/result/reportTB_filter.bim > ~/SEAsnake/result/lpWGS_snp_pos_maf0.txt

for i in {1..22};
do
  echo $i
  plink2 --bfile ~/SEAsnake/result/6_plink/reportb_wgs_chr"$i" \
    --extract range ~/SEAsnake/result/lpWGS_snp_pos_maf0.txt \
    --out ~/SEAsnake/result/7_plink_filter/reportb_wgs_chr"$i"_lpFilter \
    --make-bed --threads 10
done

aws s3 sync ~/SEAsnake/result/7_plink_filter/ s3://hawn-reporttb-wgs2/7_plink_filter/

## Merge into 1
rm ~/SEAsnake/result/7_plink_filter/files_to_merge.txt
touch ~/SEAsnake/result/7_plink_filter/files_to_merge.txt

for file in ~/SEAsnake/result/7_plink_filter/*_lpFilter.bed;
do
  NAME=$(basename "$file" .bed)
  echo $NAME >> ~/SEAsnake/result/7_plink_filter/files_to_merge.txt
done

cd ~/SEAsnake/result/7_plink_filter/
plink --merge-list ~/SEAsnake/result/7_plink_filter/files_to_merge.txt \
    --make-bed --threads 10 \
    --out ~/SEAsnake/result/7_plink_filter/reportb_wgs_all_maf0

cd ~/SEAsnake/result/
aws s3 sync ~/SEAsnake/result/7_plink_filter/ s3://hawn-reporttb-wgs2/7_plink_filter/

# 7. Filter lpWGS data for samples in hpWGS
Rscript scripts/lpWGS_hpWGS_sample_overlap.R

plink2 -bfile ~/SEAsnake/ReporTB_lpWGS/plink/reportTB_filter \
        --make-bed --out ~/SEAsnake/result/lpWGS/reportTB_filter_hpFilter \
        --keep ~/SEAsnake/result/lpWGS/hp_samples.txt \
        --threads 10

# FIN #
