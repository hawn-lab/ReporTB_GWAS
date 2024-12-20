library(tidyverse)
library(ggrepel)
setwd("~/Documents/_Hawn/ReportTB/ReporTB_GWAS/ReporTB_lpWGS")

# Metadata
## 1000 genomes
## https://github.com/laura-budurlean/PCA-Ethnicity-Determination-from-WGS-Data
kg <- read.csv("result/1000genomes/sample_info_file.csv",
             header=TRUE, stringsAsFactors=FALSE) %>% 
  filter(!Sample %in% c("your_sample_1","your_sample_2"))

## ReporTB
load("../metadata/ReportTB_metadata_clean.RData")

tb <- meta_FULL %>% 
  drop_na(sample_id) %>% 
  rename(Sample=sample_id, Gender=sex, Race=race) %>% 
  mutate(Family.ID=Sample, Population="ReporTB",
         Gender = recode(Gender, "M"="male", "F"="female")) %>% 
  select(Sample, Family.ID, Population, Gender, Race)

meta <- bind_rows(kg, tb) %>% 
  #larger population groups
  mutate(SuperPopulation = case_when(
    Population %in% c("CHB","JPT","CHS","CDX","KHV") ~ "EAS",
    Population %in% c("CEU","TSI","FIN","GBR","IBS") ~ "EUR",
    Population %in% c("YRI","ASW","ACB","LWK","GWD","MSL","ESN") ~ "AFR",
    Population %in% c("MXL","PUR","CLM","PEL") ~ "AMR",
    Population %in% c("GIH","PJL","BEB","STU","ITU") ~ "SAS",
    TRUE ~ Population
  )) %>% 
  mutate(SuperPopulation = factor(SuperPopulation),
         SuperPopulation = fct_relevel(SuperPopulation, "ReporTB", after = Inf))

# PCA data
pc1 <- read.table("result/1000genomes/5_pca/merged.1kg.reportb.PC1.sscore",
                header=FALSE, stringsAsFactors=FALSE) %>% 
  separate(V2, into=c("Sample"), sep="_", extra = "drop") %>% 
  rename(PC1=V5) %>% 
  select(Sample,PC1)
pc2 <- read.table("result/1000genomes/5_pca/merged.1kg.reportb.PC2.sscore",
                header=FALSE, stringsAsFactors=FALSE) %>% 
  separate(V2, into=c("Sample"), sep="_", extra = "drop") %>% 
  rename(PC2=V5) %>% 
  select(Sample,PC2)

# Combine data in data frame
dat <- full_join(meta, pc1, by="Sample") %>% 
  full_join(pc2, by="Sample") %>% 
  select(Family.ID, Sample, PC1, PC2, SuperPopulation, Race) %>% 
  drop_na(PC1) %>% 
  mutate(group = ifelse(SuperPopulation == "ReporTB", "ReporTB","1000 genomes"))

# Plot 
p1 <- dat %>% 
  # filter(SuperPopulation != "ReporTB") %>% 
  
  ggplot(aes(x = PC1, y = PC2, color = SuperPopulation)) +
  # geom_point(size = 3, shape = 16, alpha = 0.2) +
  geom_point(data = dat %>% filter(SuperPopulation != "ReporTB"),
             size = 3, shape = 16, alpha = 0.2) +
  geom_point(data = dat %>% filter(SuperPopulation == "ReporTB"),
             size = 3, shape = 16, alpha = 0.2) +
  scale_color_manual(values=c("#AA4499","#DDCC77","#88CCEE",
                              "#CC6677","#117733","black")) + # #332288 alt blue
  coord_fixed() +
  # facet_wrap(~group) +
  theme_bw()

# p1
ggsave("result/1000genomes/PCA.1kg.reportb_overlay.png", p1,
       width=10, height=5) 
