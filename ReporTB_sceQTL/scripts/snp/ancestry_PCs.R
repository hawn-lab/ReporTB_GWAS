#### ANCESTRY ####
library(tidyverse)
library(SNPRelate)
library(GENESIS)
library(GWASTools)
dir.create("result/kinship/", showWarnings = FALSE)

print("Data")
#### Convert to GDS ####
# No HWE filtering!!
# LD filtered
# MAF 1% 
GDS_file <- "result/gds/reportTB_sc_filter_imp_LDprune.gds"
#Check if file exists
test.files <- list.files(path="result/gds/", 
                         pattern=".gds")

if(! "reportTB_sc_filter_imp_LDprune.gds" %in% test.files){
  SNPRelate::snpgdsBED2GDS(bed.fn="result/plink/reportTB_sc_filter_imp_LDprune.bed", 
                           fam.fn="result/plink/reportTB_sc_filter_imp_LDprune.fam", 
                           bim.fn="result/plink/reportTB_sc_filter_imp_LDprune.bim",
                           out.gdsfn=GDS_file)
}

#### Read in GDS ####
genofile <- SNPRelate::snpgdsOpen(GDS_file)

#### IBD ####
print("IBD")
# Estimate IBD coefficients
king <- SNPRelate::snpgdsIBDKING(genofile, 
                                 type="KING-robust",
                                 maf=0.01, missing.rate=0.10,
                                 num.thread = 50)
king.mat <- GENESIS::kingToMatrix(king)

king_k0k1 <- SNPRelate::snpgdsIBDKING(genofile,
                                 type="KING-homo",
                                 maf=0.01, missing.rate=0.10,
                                 num.thread = 50)

save(king, king_k0k1, file="result/kinship/reportTB_sc_KING.RData")

#### PCA ####
print("PCA")
PC.air <- GENESIS::pcair(genofile, kinobj = king.mat, divobj = king.mat,
                         num.cores=50)
## 53 unrelated

#Get variance proportions
pc.percent <- PC.air$varprop*100
pc.percent <- pc.percent[!is.na(pc.percent)]

#Make dataframe of all PCs with % > 0
PC.dat <- PC.air$vectors
rownames(PC.dat) <- PC.air$sample.id
colnames(PC.dat) <- paste0("PC", 1:ncol(PC.dat))

save(PC.air, pc.percent, PC.dat,
     file="result/kinship/reportTB_sc_PCair.RData")

SNPRelate::snpgdsClose(genofile)

#### Kinship ####
print("Kinship")
load("result/kinship/reportTB_sc_KING.RData")
load("result/kinship/reportTB_sc_PCair.RData")
# create a GenotypeData class object
geno <- GWASTools::GdsGenotypeReader(GDS_file)
genoData <- GWASTools::GenotypeData(geno)
#Iterate over SNP
geno_iter <- GWASTools::GenotypeBlockIterator(genoData)

#Calc GRM using PCs 1 and 2
BPPARAM <- BiocParallel::bpparam("MulticoreParam")
BPPARAM$workers <- 6 

PC.relate <- GENESIS::pcrelate(geno_iter, pcs = PC.dat[,1:2], 
                               training.set = PC.air$unrels,
                               BPPARAM = BPPARAM)
save(PC.relate, file="result/kinship/reportTB_sc_PCrelate.RData")
## 53 training, 12393807 SNPs

kinship <- as.data.frame(as.matrix(
  pcrelateToMatrix(PC.relate, scaleKin=2))) %>% 
  #clean names
  rownames_to_column() %>% 
  separate(rowname, into=c("rowname"), sep="_", extra="drop")
colnames(kinship) <- gsub("_[0-9]{1,10}","",colnames(kinship))

write_csv(kinship, "result/kinship/reportTB_sc_kin.csv")

#### Plot results ####
print("Plots")
#Plot k0k1
#Make dataframe
k0 <- king_k0k1$k0
colnames(k0) <- king_k0k1$sample.id
rownames(k0) <- king_k0k1$sample.id
k0[upper.tri(k0)] <- NA
k0_long <- as.data.frame(k0) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, values_to = "k0") %>% 
  filter(rowname != name) %>% 
  drop_na(k0)

k1 <- king_k0k1$k1
colnames(k1) <- king_k0k1$sample.id
rownames(k1) <- king_k0k1$sample.id
k1[upper.tri(k1)] <- NA
k1_long <- as.data.frame(k1) %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, values_to = "k1") %>% 
  filter(rowname != name) %>% 
  drop_na(k1)

k0k1 <- full_join(k0_long,k1_long) %>% 
  distinct() %>% 
  arrange(k0,k1)

p1 <- ggplot(k0k1, aes(x=k0, y=k1)) + geom_point()+
  coord_fixed() +
  xlab("k0")+ ylab("k1")+
  geom_abline(intercept = 1, slope = -1, colour='red')+
  ggtitle("KING method: Identity by Descent") +
  theme_classic()

ggsave(p1, filename="result/kinship/reportTB_sc_k0vk1.png",
       width=6, height=6)

#Plot PC contributions
p2 <- data.frame(PC = 1:length(pc.percent),
                 perc = pc.percent) %>% 
  # mutate(PC=factor(PC, levels=paste0("PC",1:length(pc.percent)))) %>% 
  ggplot(aes(x=PC, y=perc)) + 
  geom_bar(stat = "identity") +
  labs(x="PC",y="Percent variation explained") +
  theme_classic()

ggsave(p2, filename="result/kinship/reportTB_sc_PCperc.png",
       width=6, height=4)

#Plot first 5 eigenvectors
library(GGally)
#metadata
load("../metadata/ReportTB_metadata_clean.RData")
meta_sub <- as.data.frame(PC.dat) %>% 
  rownames_to_column("sample_id") %>% 
  separate(sample_id, into=c("sample_id"), sep="_", extra="drop") %>% 
  inner_join(meta_FULL) %>% 
  mutate(race = as.character(race),
         race = recode(race, "6"="other"),
         race = ifelse(is.na(race),"other",race),
         race = as.factor(race),
         race = fct_recode(race, "other"="indian",
                           "other"="asian")) %>% 
  #rename EV with % explained
  select(sample_id, PC1:PC10, race) %>% 
  pivot_longer(PC1:PC10) %>% 
  mutate(name = paste0(name, "\n", round(pc.percent[1:5],2), "%")) %>% 
  pivot_wider()

#Have to drop 1 other individual
p3 <-meta_sub %>% 
  filter(race !="other") %>%
  ggpairs(columns = c(3:7), ggplot2::aes(colour=race)) + theme_bw()
# p3

ggsave(p3, filename="result/kinship/reportTB_sc_PCAx5.png",
       width=10, height=10)

p3b <- meta_sub %>% 
  filter(race !="other") %>%
  ggpairs(columns = c(3:12), ggplot2::aes(colour=race)) + theme_bw()

ggsave(p3b, filename="result/kinship/reportTB_sc_PCAx10.png",
       width=10, height=10)

#### FIN ####
print("FIN")
