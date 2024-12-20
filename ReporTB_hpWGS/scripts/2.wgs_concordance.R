# nohup Rscript 2.wgs_concordance.R >> R.log 2>&1 &

library(tidyverse)
# library(patchwork)
library(SNPRelate)
library(foreach)
library(parallel)
library(doParallel)
# library(ggbeeswarm)

setwd("/home/ec2-user/SEAsnake/result/")

#### High pass data ####
print("Load hpWGS")
GDS_file <- "8_gds/reportb_wgs_all_maf0.gds"

test.files <- list.files(path="8_gds/", pattern=".gds")
if(! "reportb_wgs_all_maf0.gds" %in% test.files){
  SNPRelate::snpgdsBED2GDS(
    bed.fn="7_plink_filter/reportb_wgs_all_maf0.bed",
    fam.fn="7_plink_filter/reportb_wgs_all_maf0.fam",
    bim.fn="7_plink_filter/reportb_wgs_all_maf0.bim",
    out.gdsfn=GDS_file)
}

genofile <- snpgdsOpen(GDS_file)
geno.info <- snpgdsGetGeno(genofile, with.id=TRUE)
geno.num <- geno.info$genotype
## Add row and column names
colnames(geno.num) <- geno.info$snp.id
sample.id <- gsub("/home/ec2-user/SEAsnake/result/3_bam_filter/|_S1_Aligned.sortedByCoord.filter.bam",
                  "",geno.info$sample.id)
rownames(geno.num) <- sample.id
snpgdsClose(genofile)

#### Low pass data ####
print("Load lpWGS")
GDS_file2 <- "lpWGS/reportTB_filter.gds"

test.files2 <- list.files(path="lpWGS/", pattern=".gds")
if(! "reportTB_filter.gds" %in% test.files2){
  SNPRelate::snpgdsBED2GDS(
    bed.fn="lpWGS/reportTB_filter_hpFilter.bed",
    fam.fn="lpWGS/reportTB_filter_hpFilter.fam",
    bim.fn="lpWGS/reportTB_filter_hpFilter.bim",
    out.gdsfn=GDS_file2)
}

genofile2 <- snpgdsOpen(GDS_file2)
geno.info2 <- snpgdsGetGeno(genofile2, with.id=TRUE)
geno.num2 <- geno.info2$genotype
## Add row and column names
colnames(geno.num2) <- geno.info2$snp.id
sample.id2 <- sapply(strsplit(sample.id,"_"), `[`, 1)
rownames(geno.num2) <- sample.id2
snpgdsClose(genofile2)

dim(geno.num)
dim(geno.num2)

#Clean up env
rm(geno.info2)
rm(geno.info)
gc()

#### Compare genotypes ####
#1=hpWGS
#2=lpWGS

cl <- parallel::makeCluster(60)
doParallel::registerDoParallel(cl)

total_chunk <- 1000
total_snp <- ncol(geno.num2)
length_of_chunk <- round(total_snp/total_chunk)

print("Test concordance")
for(chunk in 1:total_chunk){
  print(chunk)
  start<-Sys.time()
  #Select columns to process
  if(chunk==1){
    start <- 1
  } else{
    start <- length_of_chunk*(chunk-1)+1
  }
  if(chunk==total_chunk){
    end <- total_snp
  } else{
    end <- length_of_chunk*chunk
  }
  
  #Low pass data
  geno.num2.temp <- geno.num2[,start:end]
  lowSNP_alt <- gsub("/:[A-Z]$","/:.",colnames(geno.num2.temp))
  lowSNP <- c(colnames(geno.num2.temp),lowSNP_alt)
  #high pass data
  geno.num.temp <- geno.num[,colnames(geno.num) %in% lowSNP]
  
  results <- foreach(i=1:ncol(geno.num2.temp), .combine=rbind,
                     .packages = c("dplyr", "magrittr", "stats", "stringr", 
                                   "tidyr", "foreach", "doParallel",
                                   "tibble"),
                     .noexport=c("geno.num","geno.num2")) %dopar% {
                       #Define SNP name in low pass
                       snp <- colnames(geno.num2.temp)[i]
                       pos <- paste(stringr::str_split(snp, pattern=":")[[1]][1:3], collapse=":")
                       
                       #check if exact SNP name in high pass
                       if(snp %in% colnames(geno.num.temp)){
                         temp <- as.data.frame(geno.num2.temp[,snp]) %>% 
                           rownames_to_column("sample.id") %>% 
                           rename("low"=2) %>% 
                           full_join(as.data.frame(geno.num.temp[,snp]) %>% 
                                       rownames_to_column("sample.id") %>% 
                                       rename("high"=2),
                                     by = join_by(sample.id)) %>% 
                           mutate(snpLow=snp, snpHigh=snp, exact="exact") %>% 
                           mutate(concordance = case_when(is.na(low)|is.na(high)~"NA",
                                                          low==high~"TRUE",
                                                          low!=high~"FALSE",
                                                          TRUE~"NA"),
                                  concordance = factor(concordance, 
                                                       levels=c("TRUE","FALSE","NA"))) %>% 
                           count(snpLow,snpHigh,exact,concordance, .drop=FALSE) %>% 
                           pivot_wider(names_from = concordance, values_from = n) %>% 
                           mutate(tot_calls = `TRUE`+`FALSE`,
                                  pct = `TRUE`/tot_calls*100)
                       } else {
                         # If no exact match, look for SNP in same POS with missing ALT allele
                         snp2 <- colnames(geno.num.temp)[
                           grepl(paste0("^",pos), colnames(geno.num.temp))]
                         
                         # If 1 alternate SNP name found
                         if(length(snp)==1 & length(snp2)==1){
                           temp <- as.data.frame(geno.num2.temp[,snp]) %>% 
                             rownames_to_column("sample.id") %>% 
                             rename("low"=2) %>% 
                             full_join(as.data.frame(geno.num.temp[,snp2]) %>% 
                                         rownames_to_column("sample.id") %>% 
                                         rename("high"=2),
                                       by = join_by(sample.id)) %>% 
                             mutate(snpLow=snp, snpHigh=snp2, 
                                    exact = case_when(snpLow==snpHigh~"exact", 
                                                      snpHigh==paste0(pos,":.")~"missing_alt",
                                                      TRUE~"mismatch")) %>% 
                             mutate(concordance = case_when(is.na(low)|is.na(high)~"NA",
                                                            low==high~"TRUE",
                                                            low!=high~"FALSE",
                                                            TRUE~"NA"),
                                    concordance = factor(concordance, levels=c("TRUE","FALSE","NA"))) %>% 
                             count(snpLow,snpHigh,exact,concordance, .drop=FALSE) %>% 
                             pivot_wider(names_from = concordance, values_from = n) %>% 
                             mutate(tot_calls = `TRUE`+`FALSE`,
                                    pct = `TRUE`/tot_calls*100)
                         } else if(length(snp)>1 | length(snp2)>1) {
                           # If multiple alternate SNPs found, pivot data
                           temp <- as.data.frame(geno.num2.temp[,snp]) %>% 
                             rownames_to_column("sample.id") %>% 
                             pivot_longer(-sample.id, names_to = "snpLow", values_to = "low") %>% 
                             full_join(as.data.frame(geno.num.temp[,snp2]) %>% 
                                         rownames_to_column("sample.id") %>% 
                                         pivot_longer(-sample.id, names_to = "snpHigh", 
                                                      values_to = "high"),
                                       by = join_by(sample.id)) %>% 
                             mutate(exact = case_when(snpLow==snpHigh~"exact", 
                                                      snpHigh==paste0(pos,":.")~"missing_alt",
                                                      TRUE~"mismatch")) %>% 
                             mutate(concordance = case_when(is.na(low)|is.na(high)~"NA",
                                                            low==high~"TRUE",
                                                            low!=high~"FALSE",
                                                            TRUE~"NA"),
                                    concordance = factor(concordance, 
                                                         levels=c("TRUE","FALSE","NA"))) %>% 
                             count(snpLow,snpHigh,exact,concordance, .drop=FALSE) %>% 
                             pivot_wider(names_from = concordance, values_from = n) %>% 
                             mutate(tot_calls = `TRUE`+`FALSE`,
                                    pct = `TRUE`/tot_calls*100)
                         } else{
                           #If exact match and alternatve not found in high pass
                           temp <- data.frame(snpLow=snp,
                                              snpHigh=NA,
                                              exact="missing",
                                              `TRUE`=NA,
                                              `FALSE`=NA,
                                              `NA`=NA,
                                              tot_calls=0,
                                              pct=NA)
                           colnames(temp) <- c("snpLow","snpHigh","exact","TRUE","FALSE",
                                               "NA","tot_calls","pct")
                         }
                       }
                       # write_csv(temp, file="concordance.csv", append = TRUE)
                       return(temp)
                     }
  # print(results)
  # save your foreach results
  if (chunk==1) {
    write.table(results, file= "concordance.csv", row.names=FALSE, 
                na = "", sep = ",", quote=FALSE)
  } else {
    write.table(results, file="concordance.csv", append=TRUE, 
                row.names=FALSE, col.names = FALSE, na="", sep = ",", 
                quote=FALSE)
  }
  end <- Sys.time()-start
  
  write(paste(c(paste("chunk",as.character(chunk)), "in",
                paste(as.character(end),"seconds")),
              collapse=" "), 
        file="time.log",
        append=TRUE)
}

# concord <- read_csv("concordance.csv")
