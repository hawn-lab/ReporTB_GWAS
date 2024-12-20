library(tidyverse)

#### Illumina (il) ####
il_core_exome <- read_csv("../metadata/snp_manifests/Illumina_CoreExome_v1.4_hg37.csv",
                          skip=7) %>% 
  mutate(chip="il_core_exome") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_610Quad <- read_csv("../metadata/snp_manifests/Illumina_Human610Quad_v1_hg37.1.csv",
                       skip=7) %>% 
  mutate(chip="il_610Quad") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_1M <- read_csv("../metadata/snp_manifests/Illumina_HumanHap1M_v3_hg37.1.csv",
                  skip=7) %>% 
  mutate(chip="il_1M") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_660Quad <- read_csv("../metadata/snp_manifests/Illumina_Human660WQuad_v1_hg37.1.csv",
                       skip=7) %>% 
  mutate(chip="il_660Quad") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

# il_megaex37 <- read_csv("../metadata/snp_manifests/Illumina_megaex_v1_hg37.csv",
#                         skip=7) %>% 
#   mutate(chip="il_megaex37") %>% 
#   select(Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
#   drop_na() %>% 
#   rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)
il_megaex38 <- read_csv("../metadata/snp_manifests/Illumina_megaex_v1_hg38.csv",
                        skip=7) %>% 
  mutate(chip="il_megaex38") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_omni5 <- read_tsv("../metadata/snp_manifests/Illumina_omni5_v1_hg37.txt") %>%
  mutate(chip="il_omni5", genome=37) %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo) %>% 
  select(chip, everything(), -`deCODE(cM)`)

il_omni2.5 <- read_csv("../metadata/snp_manifests/Illumina_omni2.5_v1.5_hg37.csv",
                       skip=7) %>% 
  mutate(chip="il_omni2.5") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_omni1 <- read_csv("../metadata/snp_manifests/Illumina_omni1_v1.0_hg37.csv",
                       skip=7) %>% 
  mutate(chip="il_omni1") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

il_omniExpress2 <- read_csv("../metadata/snp_manifests/Illumina_omniexpress_v1.2_hg37.csv",
                           skip=7) %>% 
  mutate(chip="il_omniExpress2") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)
il_omniExpress1 <- read_csv("../metadata/snp_manifests/Illumina_omniexpress_v1.1_hg37.csv",
                            skip=7) %>% 
  mutate(chip="il_omniExpress1") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)
il_omniExpress0 <- read_csv("../metadata/snp_manifests/Illumina_omniexpress_v1.0_hg37.csv",
                              skip=7) %>% 
  mutate(chip="il_omniExpress0") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

#Cannot get v1.1 anymore
il_ominiZhongHua <- read_csv("../metadata/snp_manifests/Illumina_omni_ZhongHua_v1.4_hg37.csv",
                        skip=7) %>% 
  mutate(chip="il_ominiZhongHua") %>% 
  select(chip, Name, Chr, MapInfo, RefStrand, GenomeBuild) %>% 
  drop_na() %>% 
  rename(snpID=Name, chr=Chr, pos=MapInfo, strand=RefStrand, genome=GenomeBuild)

#### Affymetrix (af) ####
af_500k <- read_tsv("../metadata/snp_manifests/Affymetrix_5.0_500K_hg36.1.txt",
                    skip=11) %>% 
  mutate(chip="af_500k") %>% 
  select(chip, ID, Chromosome, `Physical Position`, STRAND) %>% 
  rename(snpID=ID, chr=Chromosome, pos=`Physical Position`, strand=STRAND) %>% 
  drop_na(pos) %>% 
  filter(pos != "---") %>% 
  mutate(genome = 36.1)

af_chb1 <- read_tsv("../metadata/snp_manifests/Affymetrix_CHB1_hg37.txt",
                    skip=30) %>% 
  mutate(chip="af_chb1") %>% 
  select(chip, ID, Chromosome, `Physical Position`, Strand) %>% 
  rename(snpID=ID, chr=Chromosome, pos=`Physical Position`, strand=Strand) %>% 
  drop_na(pos) %>% 
  filter(pos != "---") %>% 
  mutate(genome = 37)

af_chb2 <- read_tsv("../metadata/snp_manifests/Affymetrix_CHB2_hg37.bed",
                    skip=11, col_names = FALSE)
colnames(af_chb2) <- c("chr","pos","end",
                       "snpID","score","strand") 
af_chb2 <- af_chb2 %>% 
  mutate(chip="af_chb2") %>% 
  select(chip, snpID, chr, pos, strand) %>% 
  mutate(genome=37)

af_6.0 <- read_tsv("../metadata/snp_manifests/Affymetrix_6.0_hg37.txt",
                   skip=10) %>% 
  mutate(chip="af_6.0") %>% 
  select(chip, ID, Chromosome, `Physical Position`, STRAND) %>% 
  rename(snpID=ID, chr=Chromosome, pos=`Physical Position`, strand=STRAND) %>% 
  drop_na(pos) %>% 
  filter(pos != "---") %>% 
  mutate(genome = 37)

af_100k <- read_tsv("../metadata/snp_manifests/Affymetrix_100K_Hind240_hg35.txt",
                    skip=26) %>% 
  bind_rows(read_tsv("../metadata/snp_manifests/Affymetrix_100K_Xba240_hg35.txt",
                     skip=26)) %>% 
  mutate(chip="af_100k") %>% 
  select(chip, ID, Chromosome, `Physical Position`, Strand) %>% 
  rename(snpID=ID, chr=Chromosome, pos=`Physical Position`, strand=Strand) %>% 
  drop_na(pos) %>% 
  filter(pos != "---") %>% 
  mutate(genome = 35)

af_limaa <- read_delim("../metadata/snp_manifests/LIMAA.all.assoc.txt.gz",
                       delim=" ", col_names = FALSE) %>%
  rename(chr=X1, pos=X2) %>%
  mutate(snpID = paste(chr,pos,sep=":")) %>%
  mutate(chip="af_limaa", genome=37) %>%
  distinct(chip, snpID, chr, pos, genome)

af_axiompel <- read_delim("../metadata/snp_manifests/AYL_AxiomPEL_11-19-16_dbSNPids.bim", col_names = FALSE) %>% 
  rename(chr=X1, pos=X4) %>% 
  mutate(snpID = paste(chr,pos,sep=":")) %>%
  mutate(chip="af_axiompel", genome=37) %>%
  distinct(chip, snpID, chr, pos, genome)

#### Combine ####
chip_all <- data.frame()

for(df in ls(pattern="^il_|^af_")){
  print(df)
  chip_all <- get(df) %>% 
    mutate(chr = as.character(chr),
           chr = gsub("^chr","",chr),
           pos = as.numeric(pos)) %>% 
    bind_rows(chip_all)
  rm(list=df)
}
gc() #free unused RAM

unique(chip_all$genome)
chip_all <- chip_all %>% 
  mutate(pos = as.numeric(pos))

#### Liftover signif SNP ####
#Signif SNP
# genomes needed: 35.0 36.1 37.0 37.1 38.0
signif38 <- read_csv("result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  distinct(snpID, CHR, POS) %>% 
  rename(chr=CHR, pos=POS) %>% 
  mutate(genome=38) %>% 
  mutate(chr=as.character(chr)) %>% 
  mutate(snpID2 = paste0("chr",chr,":",pos+1,"-",pos))

#Save bed format for liftover
signif38 %>% 
  mutate(liftCHR = paste0("chr",chr),
         liftstart=pos,
         liftend=pos) %>% 
  select(contains("lift")) %>% 
  write_tsv("result/liftover/signif_snp_hg38.bed", 
            col_names = FALSE, quote = "none")
  
key <- signif38 %>% select(snpID,snpID2)

signif37 <- read_tsv("result/liftover/hglft_hg37_genome_40735_5c9bf0.bed",
                 col_names = c("chr","pos","trash","snpID2", "trash2")) %>% 
  select(snpID2, chr, pos) %>% 
  mutate(genome=37) %>% 
  mutate(chr=as.character(chr)) %>% 
  left_join(key) %>% 
  select(snpID, chr, pos, genome)

signif37.1 <- signif37 %>% 
  mutate(genome=37.1) %>% 
  mutate(chr=as.character(chr))

#Have to go from 37 to these, NOT directly from 38
signif36.1 <- read_tsv("result/liftover/hglft_hg36_genome_150be_5ca0a0.bed",
                       col_names = c("chr","pos","trash","snpID2", "trash2")) %>% 
  select(snpID2, chr, pos) %>% 
  mutate(genome=36.1) %>% 
  mutate(chr=as.character(chr)) %>% 
  left_join(key) %>% 
  select(snpID, chr, pos, genome)

signif35 <- read_tsv("result/liftover/hglft_hg35_genome_158e7_5ca220.bed",
                     col_names = c("chr","pos","trash","snpID2", "trash2")) %>% 
  select(snpID2, chr, pos) %>% 
  mutate(genome=35) %>% 
  mutate(chr=as.character(chr)) %>% 
  left_join(key) %>% 
  select(snpID, chr, pos, genome)
 
#Combine in a key for use elsewhere
signif_key <- signif38 %>% 
  rename(chr_hg38=chr, pos_hg38=pos) %>% 
  select(-genome, -snpID2) %>% 
  full_join(signif37 %>% 
              rename(chr_hg37=chr, pos_hg37=pos) %>% 
              select(-genome)) %>% 
  full_join(signif36.1 %>% 
              rename(chr_hg36=chr, pos_hg36=pos) %>% 
              select(-genome)) %>% 
  full_join(signif35 %>% 
              rename(chr_hg35=chr, pos_hg35=pos) %>% 
              select(-genome)) %>% 
  mutate(across(contains("chr"), ~gsub("chr","",.))) %>% 
  arrange(as.numeric(chr_hg38),pos_hg38)

write_csv(signif_key, "result/liftover/ReporTB_signif_snp_liftover.csv")

#Combine for this analysis
sig_all <- data.frame()
for(df in ls(pattern = "^signif3")){
  print(df)
  sig_all <- get(df) %>% 
    mutate(chr = as.character(chr),
           chr = gsub("^chr","",chr)) %>% 
    bind_rows(sig_all)
  # rm(list=df)
}

gc() #free unused RAM

#### Compare ####
# Skipped
# Illumina HumanCNV370 - copy number variants only
# Illumina HumanHap550 - could not find
# Illumina HumanHap300 - can't download manifest https://support.illumina.com/ko-kr/downloads/humanhap300_product_file.html

temp <- sig_all %>% 
  select(-snpID2) %>% 
  left_join(chip_all %>% mutate(chip=as.factor(chip)),
            by=c("chr","pos","genome")) %>% 
  select(snpID.x, snpID.y, chip) %>% 
  arrange(snpID.x, chip) %>% 
  pivot_wider(names_from = chip, values_from = snpID.y,
              values_fn = list, names_expand = TRUE) %>% 
  select(-`NA`) %>% 
  rename(snpID=snpID.x) %>% 
  select(snpID, everything()) %>% 
  mutate(across(af_100k:il_omniExpress2, ~as.character(.))) %>% 
  mutate(across(af_100k:il_omniExpress2, ~ifelse(.=="NULL" | is.null(.),
                                           NA, .))) %>% 
  mutate(none = rowSums(!is.na(select(., -snpID))),
         none = ifelse(none==0,"Y","N")) %>% 
  select(snpID, none, everything()) %>% 
  arrange(none,snpID)

#Add concordance
write_csv(temp, file = "result/liftover/compare_to_SNPchip.csv", na="")

