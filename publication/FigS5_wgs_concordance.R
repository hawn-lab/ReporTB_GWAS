library(tidyverse)
library(patchwork)

#### PCA ####
#1000 genomes PCA
# Metadata
## 1000 genomes
## https://github.com/laura-budurlean/PCA-Ethnicity-Determination-from-WGS-Data
kg <- read.csv("../ReporTB_lpWGS/result/1000genomes/sample_info_file.csv",
               header=TRUE, stringsAsFactors=FALSE) %>% 
  filter(!Sample %in% c("your_sample_1","your_sample_2"))

## ReporTB
load("../metadata/ReportTB_metadata_clean.RData")

tb <- meta_FULL %>% 
  drop_na(sample_id) %>% 
  rename(Sample=sample_id, Gender=sex, Race=race) %>% 
  mutate(Family.ID=Sample, Population="RePORT-Brazil",
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
  )) 

# PCA data
pc1 <- read.table("../ReporTB_lpWGS/result/1000genomes/5_pca/merged.1kg.reportb.PC1.sscore",
                  header=FALSE, stringsAsFactors=FALSE) %>% 
  separate(V2, into=c("Sample"), sep="_", extra = "drop") %>% 
  rename(PC1=V5) %>% 
  select(Sample,PC1)
pc2 <- read.table("../ReporTB_lpWGS/result/1000genomes/5_pca/merged.1kg.reportb.PC2.sscore",
                  header=FALSE, stringsAsFactors=FALSE) %>% 
  separate(V2, into=c("Sample"), sep="_", extra = "drop") %>% 
  rename(PC2=V5) %>% 
  select(Sample,PC2)

#30X samples
hpWGS <- read_csv("../ReporTB_lpWGS/result/other/ReportTB_samp_30X.csv")

# Combine data in data frame
dat <- full_join(meta, pc1, by="Sample") %>% 
  full_join(pc2, by="Sample") %>% 
  select(Family.ID, Sample, PC1, PC2, SuperPopulation, Race) %>% 
  drop_na(PC1) %>% 
  mutate(group = ifelse(SuperPopulation == "RePORT-Brazil", "RePORT-Brazil","1000 Genomes")) %>% 
  mutate(SuperPopulation = recode(SuperPopulation,
                                  "AFR"="African (AFR)",
                                  "AMR"="Admixed American (AMR)",
                                  "EAS"="Asian, East (EAS)",
                                  "EUR"="European (EUR)",
                                  "SAS"="Asian, South (SAS)")) %>% 
  mutate(Race = case_when(is.na(Race)~"Unknown",
                          Race=="6"~"Unknown",
                          Race=="asian"~"Other (Asian, Indian)",
                          Race=="indian"~"Other (Asian, Indian)",
                          TRUE ~ str_to_title(Race)),
         Race = factor(Race, levels=c("Brown","Black",
                                      "Other (Asian, Indian)",
                                      "White",
                                      "Unknown"))) %>% 
  mutate(hpWGS = case_when(
    Sample %in% c(hpWGS$contact, hpWGS$culture_confirm) ~ Race,
    TRUE~"Not sequenced\nin 30X hpWGS"),
    hpWGS = factor(hpWGS, levels=c("Brown","Black",
                                   "Other (Asian, Indian)",
                                   "White",
                                   "Unknown", 
                                   "Not sequenced\nin 30X hpWGS")))

# Plot 
x_set <- c(min(dat$PC1), max(dat$PC1))
y_set <- c(min(dat$PC2), max(dat$PC2))

p0 <- dat %>% 
  filter(group=="RePORT-Brazil") %>% 
  # mutate(group = "RePORT-Brazil\nhpWGS samples") %>% 
  mutate(group_ord = case_when(hpWGS != "Not sequenced\nin 30X hpWGS"~1,
                               TRUE~2)) %>% 
  arrange(-group_ord) %>% 
  ggplot(aes(x = PC1, y = PC2, color = hpWGS)) +
  geom_point(size = 2, shape = 16, 
             aes(alpha = as.character(group_ord))) +
  scale_color_manual(values=c("#DDCC77","#AA4499",
                              "#117733",
                              "#CC6677","grey")) +
  coord_fixed() +
  theme_classic()  +
  # facet_wrap(~group) +
  labs(color = "Self-reported race") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  lims(x=x_set, y=y_set) +
  # theme(legend.position = "bottom")+
  scale_alpha_manual(values=c(1, 0.2), guide="none")
# p0

#### Global plots ####
maf0 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter.bim",
                 col_names=FALSE) %>% pull(X2)
maf0_imp <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_imp.bim",
                 col_names=FALSE) %>% pull(X2)
#MAF > 5% and imputation score > 0.6
maf5 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_imp_maf5.bim",
                 col_names=FALSE) %>% pull(X2)

#concordance
concord_maf0 <- read_csv("../ReporTB_hpWGS/result/concordance.csv.gz") %>%
  filter(snpLow %in% maf0)
concord_maf5 <- concord_maf0 %>% filter(snpLow %in% maf5)

#### presence/absence ####
nrow(concord_maf0)
p1 <- concord_maf0 %>% 
  mutate(maf_group = case_when(snpLow %in% maf5 ~ "Impute & MAF pass-filter",
                               snpLow %in% maf0_imp ~ "Impute pass-filter",
                               snpLow %in% maf0 ~ "Fail",
                               TRUE~"None"),
         maf_group = factor(maf_group, levels=c(
           "Impute & MAF pass-filter",
           "Impute pass-filter",
           "Fail"))) %>% 
  #Count participants with successful call in high-pass
  count(maf_group, tot_calls) %>% 
  #Split N=50 all called, from N<50 some called
  mutate(group = case_when(tot_calls>0~"present",
                           tot_calls==0~"missing")) %>% 
  group_by(group, maf_group) %>% 
  summarise(n=sum(n)) %>% 
  ungroup() %>% 
  #Calc percent of total SNP
  mutate(pct = n/sum(n)*100) %>% 
  #Reorder x axis
  mutate(group = factor(group, levels=rev(c("present",
                                            "missing")))) %>%
  #Create % label
  mutate(pct_lab = ifelse(pct>=1,
                          paste0(round(pct, digits=0),"%"),
                          paste0(round(pct, digits=0),"%")))%>% 
  arrange(desc(pct_lab)) %>% 
  ggplot(aes(x=group,y=n)) + 
  geom_bar(aes(fill=maf_group), stat="identity") +
  # geom_text(data=dat_p1_summ, aes(label=pct_lab), vjust=-0.2, size=3) +
  geom_text(aes(label=pct_lab), position=position_stack(vjust=0.5), 
            size=3) +
  theme_bw() +
  labs(y="lpWGS SNPs\n(N = 49,493,544)",
       x="Presence in\n30X hpWGS",
       fill="lpWGS filters") +
  lims(y=c(0,3.5E7)) +
  scale_fill_manual(values=c("#117733","#88CCEE","#CC6677")) +
  theme(legend.position = "none")
# p1

#### Concordance plot ####
# Keep only SNP called in all high-pass participants
concord2 <- concord_maf0 %>% 
  filter(tot_calls==50) %>% 
  mutate(maf_group = case_when(snpLow %in% maf5 ~ "Impute & MAF pass-filter",
                               snpLow %in% maf0_imp ~ "Impute pass-filter",
                               snpLow %in% maf0 ~ "Fail",
                               TRUE~"None"),
         maf_group = factor(maf_group, levels=c(
           "Impute & MAF pass-filter",
           "Impute pass-filter",
           "Fail")))
nrow(concord2)

## Total % of successful calls concord
mean_concord <- signif(mean(concord2$pct), digits=2)
sd_concord <- signif(sd(concord2$pct), digits=2)

#histogram
p2 <- concord2 %>% 
  # filter(snpLow %in% maf5) %>% 
  ggplot(aes(pct)) +
  geom_histogram(bins=50, aes(fill=maf_group), alpha=0.4, position="identity") +
  theme_bw() +
  labs(x="Percent concordant participants",
       y="Overlapping SNPs\n(N = 15,472,133)",
       fill="lpWGS filters") +
  #label mean
  geom_vline(xintercept = mean_concord, color="red") +
  annotate("text", label=paste0("mean = ",mean_concord,"%"),
           x=mean_concord, y=1.75e6, hjust=1,
           color="red", size=3) +
  #relabel y-axis
  scale_y_continuous(labels = function(x) formatC(x, format = "G", digits=1,
                                                 drop0trailing=TRUE)) +
  scale_fill_manual(values=c("#117733","#88CCEE","#CC6677")) +
  theme(legend.position = "none")
# p1+p2

#### Signif SNP plots ####
signif_snp <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value" & pval < 5E-8)

#all sig snp
snpSig <- unique(signif_snp$snpID)

#### presence/absence ####
p1b <- concord_maf5 %>% 
  filter(snpLow %in% snpSig) %>% 
  mutate(maf_group = case_when(snpLow %in% maf5 ~ "Impute & MAF pass-filter",
                               snpLow %in% maf0_imp ~ "Impute pass-filter",
                               snpLow %in% maf0 ~ "Fail",
                               TRUE~"None"),
         maf_group = factor(maf_group, levels=c(
           "Impute & MAF pass-filter",
           "Impute pass-filter",
           "Fail"))) %>% 
  #Count participants with successful call in high-pass
  count(maf_group, tot_calls) %>% 
  #Split N=50 all called, from N<50 some called
  mutate(group = case_when(tot_calls>0~"present",
                           tot_calls==0~"missing")) %>% 
  group_by(group, maf_group) %>% 
  summarise(n=sum(n)) %>% 
  ungroup() %>% 
  #Calc percent of total SNP
  mutate(pct = n/sum(n)*100) %>% 
  #Reorder x axis
  mutate(group = factor(group, levels=rev(c("present",
                                            "missing")))) %>%
  #Create % label
  mutate(pct_lab = ifelse(pct>=1,
                          paste0(round(pct, digits=0),"%"),
                          paste0(round(pct, digits=0),"%"))) %>% 
  arrange(desc(pct_lab)) %>% 
  ggplot(aes(x=group,y=n)) + 
  geom_bar(aes(fill=maf_group), stat="identity") +
  # geom_text(data=dat_p1_summ, aes(label=pct_lab), vjust=-0.2, size=3) +
  geom_text(aes(label=pct_lab), position=position_stack(vjust=0.5), 
            size=3) +
  theme_bw() +
  labs(y="lpWGS significant\nSNPs (N = 68)",
       x="Presence in\n30X hpWGS",
       fill="lpWGS filters") +
  # lims(y=c(0,70)) +
  scale_fill_manual(values=c("#117733","#88CCEE","#CC6677")) +
  theme(legend.position = "none")
# p1b

#### Concordance plot ####
# Keep only SNP called in all high-pass participants
concord2_maf5 <- concord_maf5 %>% 
  filter(snpLow %in% snpSig & tot_calls==50) %>% 
  mutate(maf_group = case_when(snpLow %in% maf5 ~ "Impute & MAF pass-filter",
                               snpLow %in% maf0_imp ~ "Impute pass-filter",
                               snpLow %in% maf0 ~ "Fail",
                               TRUE~"None"),
         maf_group = factor(maf_group, levels=c(
           "Impute & MAF pass-filter",
           "Impute pass-filter",
           "Fail")))
nrow(concord2_maf5)

## Total % of successful calls concord
mean_concord2 <- signif(mean(concord2_maf5$pct), digits=2)
sd_concord2 <- signif(sd(concord2_maf5$pct), digits=2)

#histogram
p2b <- concord2_maf5 %>% 
  # filter(snpLow %in% maf5) %>% 
  ggplot(aes(pct)) +
  geom_histogram(bins=50, aes(fill=maf_group), alpha=0.4, position="identity",
                 show.legend=TRUE) +
  theme_bw() +
  labs(x="Percent concordant participants",
       y="Overlapping significant\nSNPs (N = 34)",
       fill="lpWGS filters") +
  #label mean
  geom_vline(xintercept = mean_concord2, color="red") +
  annotate("text", label=paste0("mean = ",mean_concord2,"%"),
           x=mean_concord2, y=3, hjust=1,
           color="red", size=3) +
  #relabel y-axis
  scale_y_continuous(labels = function(x) formatC(x, format = "G", digits=1,
                                                  drop0trailing=TRUE)) +
  scale_fill_manual(values=c("#117733","#88CCEE","#CC6677"), drop=FALSE) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))
# p1b+p2b

#### Per signif snp plot ####
concord_sig <- concord_maf5 %>% 
  #filter significant SNPs
  filter(snpLow %in% snpSig) %>% 
  #long format to color by concordance
  pivot_longer(`TRUE`:`NA`) %>% 
  mutate(name = recode(name, "TRUE"="concordant",
                       "FALSE"="discordant",
                       "NA"="missing"),
         name = factor(name, levels=c("missing","discordant",
                                      "concordant"))) %>% 
  replace_na(list(value=0)) %>% 
  #Reorder SNP by concordance with missing on right (if included)
  mutate(snpLow = factor(snpLow),
         snpLow = fct_reorder2(snpLow, name, value)) %>% 
  #fix NA count for SNP with 0 calls
  mutate(value = case_when(exact=="missing" & name=="missing"~50,
                           TRUE~value))

#bar plot of overlapping annotated signif SNP concordance
concord_sig %>% 
  filter(snpLow %in% snpSig & name=="concordant") %>% 
  mutate(concord = ifelse(pct>=60,"high","low")) %>% 
  distinct(snpLow, exact, concord) %>% 
  count(exact,concord)

p3 <- concord_sig %>% 
  filter(snpLow %in% snpSig) %>% 
  filter(exact!="missing" & value>0) %>%
  # filter(col.group != "no annotation") %>% 
  droplevels() %>%
  
  ggplot(aes(x=snpLow,y=value,fill=name)) +
  geom_bar(stat = "identity") +
  # geom_col() +
  theme_bw() +
  labs(fill="",
       x="Overlapping significant SNPs (N = 34)",y="Total participants") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  scale_fill_manual(values=c("#AA4499","#DDCC77","#117733")) +
  geom_hline(yintercept = 0.6*50, lty="dashed")
# p3

#### Save ####
# lo <- "
# ABB
# CDD
# EEE
# "
lo <- "
AABD
CCEE
FFFF
"
p <- p0+p1+free(p2)+p1b+free(p2b)+free(p3)+
  plot_annotation(tag_levels = "A")+
  plot_layout(design=lo, heights = 1)
# p

ggsave("FigS5_wgs_concordance.png", p,
       width=9, height = 7)

