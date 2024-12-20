library(tidyverse)
library(SNPRelate)
library(SeqVarTools)
library(foreach)
library(doParallel)
registerDoParallel(cores=5)
span <- 10E3
signif <- 5E-8
# span2 <- 5000

#### All p-values ####
maf1_concord <- read_tsv("result/snp_imp_maf1_concord.txt", col_names=FALSE) %>% 
  pull(X1)

maf5_concord <- read_tsv("result/snp_imp_maf5_concord.txt", col_names=FALSE) %>% 
  pull(X1)

fdr <- read_csv("result/model_final/model_full.csv") %>% 
  filter(snpID %in% maf1_concord) %>% 
  filter(term == "value") %>% 
  select(snpID, pval)

scores <- read_tsv("result/model_final/glmm_score_full.txt") %>% 
  filter(SNP %in% maf1_concord) %>% 
  filter(!SNP %in% fdr$snpID) %>% 
  rename(snpID=SNP, pval=PVAL) %>% 
  select(snpID, pval)

dat <- bind_rows(fdr,scores) %>% 
  separate(snpID, into=c("CHR","POS"), sep=":",
           remove = FALSE, extra = "drop") %>% 
  mutate(CHR = as.numeric(CHR),
         POS = as.numeric(POS)) 

#non-signif SNP in dataset
dat_ns <- dat %>% 
  filter(pval>=signif)

#signif SNP
dat_sig <- dat %>% 
  filter(pval<signif) %>% 
  filter(snpID %in% maf5_concord) %>% 
  #Add LD region
  mutate(min.POS = POS-span,
         max.POS = POS+span)

#### SNP within span of signif SNP ####
snp.regions <- dat %>% 
  filter(CHR %in% dat_sig$CHR) 

snp.regions.signif <- snp.regions %>% 
  inner_join(dat_sig %>% 
               select(CHR, max.POS, min.POS), 
             by="CHR", relationship = "many-to-many") %>% 
  filter(POS >= min.POS & POS <= max.POS)

#### SNP data ####
# Extract and format genotypes
## MAF > 1%
GDS_file1 <- "result/gds/reportTB_filter_imp_maf1.gds"
genofile1 <- snpgdsOpen(GDS_file1)
#List SNP in this gds file
snp.id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))
snp.id1 <- snp.id1[snp.id1 %in% unique(snp.regions.signif$snpID)]

geno.dat <- snpgdsGetGeno(genofile1, snp.id=snp.id1, with.id = TRUE)

## Add row and column names
geno.num <- geno.dat$genotype
colnames(geno.num) <- geno.dat$snp.id
rownames(geno.num) <- sapply(strsplit(geno.dat$sample.id,"_"), `[`, 1)
snpgdsClose(genofile1)

all.snp <- geno.dat$snp.id

rm(fdr, scores, dat, geno.dat)
gc()

#### Calc LD ####
SNPs <- sort(dat_sig$snpID)

foreach(snp=SNPs) %dopar% {
  print(snp)
  signif.temp <- dat_sig %>% filter(snpID == snp)
  
  #### Find all SNP pairs within span kb of lead snp ####
  print("Find SNP pairs")
  #Get SNPs in region around gene
  chr.temp <- paste0("^",signif.temp$CHR,":")
  snp.temp <- data.frame(snpID = all.snp[grepl(chr.temp, all.snp)]) %>% 
    separate(snpID, into=c("CHR","POS"), sep=":", remove = FALSE,
             extra = "drop") %>% 
    mutate(POS = as.numeric(POS)) %>% 
    filter(POS >= signif.temp$min.POS & 
             POS <= signif.temp$max.POS)
  
  #### Convert genotypes to alleles ####
  geno.temp <- as.data.frame(geno.num[,snp.temp$snpID]) %>% 
    rownames_to_column() %>% 
    pivot_longer(-rowname, names_to = "snpID") %>% 
    separate(snpID, into=c("CHR","POS","REF","ALT"),
             sep=":", remove=FALSE, extra="drop") %>% 
    mutate(allele = case_when(value==0~paste(REF,REF,sep="/"),
                              value==1~paste(REF,ALT,sep="/"),
                              value==2~paste(ALT,ALT,sep="/"))) %>% 
    select(rowname, snpID, allele) %>% 
    pivot_wider(names_from = "snpID", values_from = "allele") %>% 
    column_to_rownames() %>% 
    as.matrix()
  
  #### LD just to lead SNP ####
  #Get all comparisons to lead SNP
  snps.to.LD <- data.frame(lead.snp = snp, snpID = snp.temp$snpID)

  print("Calculate LD to lead")
  LD.df <- snps.to.LD %>%
    mutate(R2="", Dprime="")

  for(i in 1:nrow(LD.df)){
    # for(i in 1:3){
    #Get SNP pair for i index
    SNP.pair <- snps.to.LD[i, ] %>% unlist(use.names = FALSE)

    #Filter allele data to SNP pairs
    if(all(SNP.pair %in% colnames(geno.temp))){
      geno.allele.sub <- geno.temp[,SNP.pair]

      #Calculate LD R^2
      #Only run if > 1 alleles in both SNPs
      no.snp1 <- unique(geno.allele.sub[,1][!is.na(geno.allele.sub[,1])]) %>%
        length()
      no.snp2 <- unique(geno.allele.sub[,2][!is.na(geno.allele.sub[,2])]) %>%
        length()

      if(no.snp1 > 1 & no.snp2 > 1){
        LD.temp <- genetics::LD(
          genetics::as.genotype(geno.allele.sub[,1]),
          genetics::as.genotype(geno.allele.sub[,2]))
        LD.df[i, "R2"] <- LD.temp[["R^2"]]
        LD.df[i, "Dprime"] <- LD.temp[["D'"]]
      } else{
        LD.df[i, "R2"] <- NA
        LD.df[i, "Dprime"] <- NA
      }
    }}

  save(LD.df, file=paste0("result/LD/", gsub(":","_",snp),".RData"))
  
  #### All pairwise comparisons####
  # All pairwise for eQTL SNPs
  if(snp%in%c("3:75816280:A:T","22:39566026:A:T","17:14155330:A:G","4:139801643:G:T")){
    print("Calculate LD to all")
    snps.to.LD2 <- as.data.frame(t(combn(snp.temp$snpID, m=2))) %>%
      filter(V1 != V2)
    LD.df.all <- snps.to.LD2 %>%
      mutate(R2="", Dprime="")

    for(i in 1:nrow(LD.df.all)){
      # for(i in 1:3){
      #Get SNP pair for i index
      SNP.pair2 <- snps.to.LD2[i, ] %>% unlist(use.names = FALSE)

      #Filter allele data to SNP pairs
      if(all(SNP.pair2 %in% colnames(geno.temp))){
        geno.allele.sub2 <- geno.temp[,SNP.pair2]

        #Calculate LD R^2
        #Only run if > 1 alleles in bpth SNPs
        no.snp3 <- unique(geno.allele.sub2[,1][!is.na(geno.allele.sub2[,1])]) %>%
          length()
        no.snp4 <- unique(geno.allele.sub2[,2][!is.na(geno.allele.sub2[,2])]) %>%
          length()

        if(no.snp3 > 1 & no.snp4 > 1){
          LD.temp2 <- genetics::LD(
            genetics::as.genotype(geno.allele.sub2[,1]),
            genetics::as.genotype(geno.allele.sub2[,2]))
          LD.df.all[i, "R2"] <- LD.temp2[["R^2"]]
          LD.df.all[i, "Dprime"] <- LD.temp2[["D'"]]
        } else{
          LD.df.all[i, "R2"] <- NA
          LD.df.all[i, "Dprime"] <- NA
        }
      }}

    #### Save ####
    save(LD.df.all, file=paste0("result/LD/", gsub(":","_",snp),".all.RData"))
  }

  print("Done")
}
beepr::beep()
