# Pull header
zcat < ~/project/data/7212-BA-imputed-merged.vcf.gz | head -n 41 > ~/project/result/vcf/7212-BA-imputed-merged.vcf.header.tsv

# Pull columns with per SNP information
# Note that QUAL (column 6) and FILTER (coloumn 7) are all ".", so we're removing them
zcat < ~/project/data/7212-BA-imputed-merged.vcf.gz | tail -n +42 | cut -f1,2,4,5,8 > ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.tsv

# Separate INFO (column 6) into 3 columns
# #CHROM	POS	REF	ALT	QUAL	INFO 
# INFO looks like this AF=0.00611801;RAF=0.000399361;INFO=0.087
awk -F '\t' 'NR==1 {
    # Print the original column names and append AF, RAF, and INFO column names
    print $1, $2, $3, $4, "AF", "RAF", "INFO"
    next
} 
{
    # Split the INFO column 6 by ";" and extract AF, RAF, INFO
    split($5, arr, ";")
    for (i in arr) {
        split(arr[i], kv, "=")
        if (kv[1] == "AF") af = kv[2]
        else if (kv[1] == "RAF") raf = kv[2]
        else if (kv[1] == "INFO") info = kv[2]
    }
    # Print the original columns followed by AF, RAF, and INFO
    print $1, $2, $3, $4, af, raf, info
}' OFS='\t' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.tsv > ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.tsv

# Remove chr from chr name
sed 's/^chr//' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.tsv > ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename.tsv
# Add snpID
awk 'BEGIN {OFS="\t"} {if (NR==1) {print "snpID", $5, $6, $7} else {print $1":"$2":"$3":"$4, $5, $6, $7}}' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename.tsv > result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv 

# Filter impute score > 0.5
awk -F"\t" '$4 > 0.5 || NR==1' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv |
  awk -F"\t" 'NR > 1 {print $1}' > ~/project/result/snp_impute_pf05.txt
# Filter impute score > 0.6
awk -F"\t" '$4 > 0.6 || NR==1' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv |
  awk -F"\t" 'NR > 1 {print $1}' > ~/project/result/snp_impute_pf06.txt
# Filter impute score > 0.7
awk -F"\t" '$4 > 0.7 || NR==1' ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv |
  awk -F"\t" 'NR > 1 {print $1}' > ~/project/result/snp_impute_pf07.txt

# Compress results files
gzip ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.tsv
gzip ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.tsv
gzip ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename.tsv
gzip ~/project/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv

# Save results to S3. Remove files from instance
aws s3 sync ~/project/result/vcf/ s3://hawn-reporttb-results/result/vcf/
