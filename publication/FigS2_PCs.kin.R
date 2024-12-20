library(tidyverse)
library(patchwork)

#### Ancestry PCA ####
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
    TRUE~"Not sequenced"),
    hpWGS = factor(hpWGS, levels=c("Brown","Black",
                                   "Other (Asian, Indian)",
                                   "White",
                                   "Unknown", "Not sequenced")))

# Plot 
x_set <- c(min(dat$PC1), max(dat$PC1))
y_set <- c(min(dat$PC2), max(dat$PC2))

p0a <- dat %>% 
  filter(group=="1000 Genomes") %>% 
  ggplot(aes(x = PC1, y = PC2, color = SuperPopulation)) +
  geom_point(size = 3, shape = 16, alpha = 0.2) +
  scale_color_manual(values=c("#DDCC77","#AA4499",
                              "#88CCEE","#117733",
                              "#CC6677","black")) + # #332288 alt blue
  coord_fixed() +
  theme_classic()  +
  facet_wrap(~group) +
  labs(color = "1000 Genomes\nSuper Population") +
  guides(colour = guide_legend(override.aes = list(alpha = 1), ncol=2)) +
  lims(x=x_set, y=y_set)+
  theme(legend.position = "bottom")

p0b <- dat %>% 
  filter(group=="RePORT-Brazil") %>% 
  mutate(group_ord = case_when(
    hpWGS != "Not sequenced"~1,
    Race %in% c("Other (Asian, Indian)","Unknown")~2,
    TRUE~3)) %>% 
  arrange(-group_ord) %>% 
  ggplot(aes(x = PC1, y = PC2, color = Race)) +
  geom_point(size = 3, shape = 16, alpha = 0.2) +
  scale_color_manual(values=c("#DDCC77","#AA4499",
                              #"#88CCEE",
                              "#117733",
                              "#CC6677","black")) + # #332288 alt blue
  coord_fixed() +
  theme_classic()  +
  facet_wrap(~group) +
  labs(color = "Self-reported\nrace") +
  guides(colour = guide_legend(override.aes = list(alpha = 1), ncol=3)) +
  lims(x=x_set, y=y_set) +
  theme(legend.position = "bottom")

# p0a+p0b

#### PC percent ####
load("../ReporTB_lpWGS/result/kinship/reportTB_PCair.RData")

#Plot PC contributions
p1 <- data.frame(PC = 1:length(pc.percent),
                 perc = pc.percent) %>% 
  # mutate(PC=factor(PC, levels=paste0("PC",1:length(pc.percent)))) %>% 
  ggplot(aes(x=PC, y=perc)) + 
  geom_bar(stat = "identity") +
  labs(x="PC",y="Percent variation explained") +
  theme_classic()
# p1

#### PC model fit ####
output <- read_csv("../ReporTB_lpWGS/result/model_fitting/model_pcs.csv")
output_PC_delta <- output %>% 
  filter(term == "value" & model %in% c("base","PC1","PC2","PC3","PC4","PC5")) %>% 
  mutate(model_fct = factor(model, levels=c("base","PC1","PC2","PC3","PC4","PC5"))) %>% 
  arrange(snpID,model_fct) %>% 
  group_by(snpID) %>% 
  mutate(delta = AIC-lag(AIC)) %>% 
  drop_na(delta)

p2 <- ggplot(output_PC_delta, aes(x=model,y=delta)) +
  geom_violin() +
  labs(x="Total PCs included in GLM",
       y="Change in AIC\nBetter fit\nwith additional PC <--   --> without additional PC") +
  geom_hline(yintercept = 0, color="red", lty="dashed") +
  theme_bw()
# p2

#### Kinship ####
kin_raw <- read_csv("../ReporTB_lpWGS/result/kinship/reportTB_kin.csv") %>% 
  column_to_rownames()

kin <- kin_raw
kin[lower.tri(kin)] <- NA

kin_long <- kin %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  drop_na() %>% 
  filter(rowname != name) %>% 
  mutate(col_group = case_when(value < 0.125 ~ "None",
                               value >= 0.125 & value < 0.25 ~ "Third-degree",
                               value >= 0.25 & value < 0.5 ~ "Second-degree",
                               value >= 0.5 ~ "First-degree")) %>% 
  mutate(col_group = factor(col_group, levels = c("First-degree",
                                                  "Second-degree",
                                                  "Third-degree",
                                                  "None")))

#Filter donors in final dataset
fam <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_imp_maf5.fam", 
                col_names = FALSE) %>% 
  select(X2) %>% 
  separate(X2, into = "sample_id", sep="_", remove = FALSE) %>% 
  rename(libID = X2)

kin_long_filter <- kin_long %>% 
  filter(rowname %in% fam$sample_id & name %in% fam$sample_id)

p3 <- kin_long_filter %>% 
  ggplot(aes(x=1, y=value)) +
  geom_violin() +
  geom_jitter(data = filter(kin_long_filter, value > 0.125), 
              aes(color = col_group)) +
  geom_hline(yintercept = c(0.125, 0.25, 0.5), 
             lty = "dashed", color = "red") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Kinship", color = "Relatedness") +
  scale_color_manual(values=c("black","grey35","grey70"))
# p3

# length(unique(c(kin_long_filter$rowname, kin_long_filter$name)))
table(kin_long_filter$col_group)

#### Kinship fit ####
output2 <- read_csv("../ReporTB_lpWGS/result/model_fitting/model_kin.csv") %>% 
  bind_rows(output)
output_kin_delta <- output2 %>% 
  filter(term == "value" & model %in% c("PC2","PC2_kinship")) %>% 
  mutate(model_fct = factor(model, levels=c("PC2","PC2_kinship"))) %>% 
  arrange(snpID,model_fct) %>% 
  group_by(snpID) %>% 
  mutate(delta = AIC-lag(AIC)) %>% 
  drop_na(delta)

p4 <- ggplot(output_kin_delta, aes(x=model,y=delta)) +
  geom_violin() +
  labs(x="GLMMkin vs GLM",
       y="Change in AIC\nBetter fit\nwith kinship <--   --> without kinship") +
  geom_hline(yintercept = 0, color="red", lty="dashed") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# p4

#### Covariate fit ####
output3 <- read_csv("../ReporTB_lpWGS/result/model_fitting/model_cov.csv") %>% 
  bind_rows(read_csv("../ReporTB_lpWGS/result/model_fitting/model_cov_minus1.csv")) %>%
  bind_rows(output2)

#calculate delta AIC for kinship
output_delta3 <- output3 %>% 
  filter(term == "value") %>% 
  filter(model %in% c("PC2_kinship", "age","sex","smoke","smoke2","hiv","diabetes","all") |
           grepl("minus", model)) %>% 
  select(model, term, snpID, AIC) %>% 
  pivot_wider(names_from = model, values_from = AIC) %>% 
  pivot_longer(-c(snpID,term,PC2_kinship,all)) %>% 
  mutate(delta = case_when(!grepl("minus", name)~value-PC2_kinship,
                           grepl("minus", name)~value-all)) %>% 
  drop_na(delta) %>% 
  mutate(group = ifelse(grepl("minus", name), "Leave-one-out", "Univariate"),
         group = factor(group, levels=c("Univariate","Leave-one-out")),
         name = gsub("_","\n",name))

p5 <- output_delta3 %>% 
  filter(!grepl("minus",name)) %>% 
  ggplot(aes(x=name,y=delta)) +
  geom_violin() +
  labs(x="GLMMkin + 2 PCs + 1 covariate", 
       y="Change in AIC\n Better fit\nwith covariate <--   --> without covariate") +
  geom_hline(yintercept = 0, color="red", lty="dashed") +
  theme_classic() +
  facet_wrap(~group, scales = "free")

p6 <- output_delta3 %>%
  filter(grepl("minus",name)) %>% 
  ggplot(aes(x=name,y=delta)) +
  geom_violin() +
  labs(x="GLMMkin + 2 PCs + sex + age + smoking + HIV + diabetes", 
       y="Change in AIC\n Better fit\nwith covariate removed <--   --> with covariate retained") +
  geom_hline(yintercept = 0, color="red", lty="dashed") +
  theme_classic() +
  facet_wrap(~group, scales = "free")

#### Save #### 
lo <- "
AAABBB
AAABBB
CCDDEF
CCDDEF
GGGHHH
GGGHHH
"
plot_all <- p0a+p0b+p1+p2+p4+p3+p5+p6 + plot_annotation(tag_levels = "A") +
  plot_layout(design=lo)
# plot_all

ggsave("FigS2_PCs.kin.png", plot_all,
       width=11, height=12.5)
