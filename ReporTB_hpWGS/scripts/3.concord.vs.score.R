# Compare hpWGS concordance vs lpWGS imputation score
library(tidyverse)
library(viridis)
library(ggExtra)

#### Data ####
#Concordance
concord <- read_csv("result/concordance.csv.gz")

#Imputation score
imp <- read_tsv("../ReporTB_lpWGS/result/vcf/7212-BA-imputed-merged.vcf.snp.split.rename2.tsv.gz")

#Add minor allele frequency info
maf5 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_maf5.bim", 
                 col_names=FALSE) %>% pull(X1)
maf1 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_maf1.bim", 
                 col_names=FALSE) %>% pull(X1)

#Combine data
identical(concord$snpLow, imp$snpID)
dat <- imp %>% 
  select(snpID, INFO) %>% 
  bind_cols(concord %>% select(snpLow, snpHigh, tot_calls, pct)) %>% 
  mutate(maf = case_when(snpID %in% maf5~"maf5",
                         snpID %in% maf1~"maf1",
                         TRUE~"maf0"))
identical(dat$snpID, dat$snpLow)

rm(concord, imp, maf1, maf5)
gc()

#Select overlap
nrow(dat)
dat %>% 
  count(tot_calls) %>% 
  arrange(-tot_calls) %>% 
  View()

dat_overlap <- dat %>% 
  filter(tot_calls>=40) %>% 
  mutate(prop = pct/100)
nrow(dat_overlap)

#### Linear model ####
# Linear model
mod <- lm(INFO~prop, data=dat_overlap)
broom::tidy(mod)

#### Determine cutoffs ####
# Spread of imputation scores
summary(dat$INFO)
dat %>% 
  ggplot(aes(INFO)) +
  geom_histogram(bins=100) +
  theme_bw() +
  geom_vline(xintercept = 0.7, color="red")

# Spread of concordance
summary(dat_overlap$prop)
dat_overlap %>% 
  ggplot(aes(prop)) +
  geom_histogram() +
  theme_bw() +
  geom_vline(xintercept = 0.7, color="red")

#### Plot all data ####
#Dot plot
p1 <- ggplot(dat_overlap, aes(x=prop, y=INFO)) +
  geom_point(alpha=0.2) +
  geom_smooth(method="lm", se=FALSE, color="red") +
  theme_bw() +
  coord_fixed() +
  labs(x="Concordance", y="Imputation score")
# p1

ggsave("figs/concord.vs.score.png", p1, height=5, width=5)

# density plot
p2 <- dat_overlap %>% 
  # slice_head(n=1000) %>%
  ggplot(aes(x=prop, y=INFO)) +
  geom_point(alpha=0) +
  geom_bin2d(bins = 50) +
  coord_fixed() +
  theme_bw() +
  labs(x="Concordance", y="Imputation score") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_fill_viridis_c(name = "Count", trans = "log",
                       breaks=c(10,100,1E3,1E4,1E5),
                       labels=formatC(c(10,100,1E3,1E4,1E5))) +
  theme(legend.text = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept=0.7) +
  geom_hline(yintercept=0.7)

p2 <- ggMarginal(p2, type="density")
# p2

ggsave("figs/concord.vs.score_density.png", p2, height=6, width=5)
beepr::beep()
