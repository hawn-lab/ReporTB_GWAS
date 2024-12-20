library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggvenn)
signif <- 5E-8

#### Data ####
maf5_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf5_concord.txt", 
                 col_names = FALSE) %>% 
  pull(X1)

#Full models
m1 <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value" & snpID %in% maf5_concord) %>% 
  select(model, snpID, pval)

m2 <- read_csv("../ReporTB_lpWGS/result/model_final/model_full_comp.csv") %>% 
  filter(term=="value" & snpID %in% maf5_concord) %>% 
  select(model, snpID, pval)

#scores
s1 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_full.txt") %>% 
  filter(SNP %in% maf5_concord & !SNP %in% m1$snpID) %>% 
  select(SNP,PVAL) %>% 
  rename(snpID=SNP, pval=PVAL) %>% 
  mutate(model="full")

s2 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_full_comp.txt") %>% 
  filter(SNP %in% maf5_concord & !SNP %in% m2$snpID) %>% 
  select(SNP,PVAL) %>% 
  rename(snpID=SNP, pval=PVAL) %>% 
  mutate(model="full_comp")

dat <- bind_rows(m1, m2, s1, s2) %>% 
  pivot_wider(names_from = model, values_from = pval) %>% 
  mutate(Significance = case_when(full<signif & full_comp<signif ~ "Both",
                                  full<signif ~ "All covariates only",
                                  full_comp<signif ~ "No TB risk covariates only",
                                  TRUE ~ "None"),
         Significance = factor(Significance,
                               levels=c("Both", "All covariates only",
                                        "No TB risk covariates only", "None")))

table(dat$Significance)

#### Compare P all genes ####
#all
p1 <- dat %>% 
  # sample_n(size = 100000) %>% 
  ggplot(aes(x=-log10(full), y=-log10(full_comp))) +
  geom_hline(yintercept = -log10(signif), lty="dashed", color="grey") +
  geom_vline(xintercept = -log10(signif), lty="dashed", color="grey") +
  stat_cor(method = "pearson", label.x = 0, label.y = 11) +
  geom_point(aes(color = Significance, alpha=Significance)) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  theme_classic() +
  coord_fixed() +
  labs(x="-log10( P )\nAll covariates",
       y="-log10( P )\nNo TB risk covariates")  +
  scale_color_manual(values = c('#88CCEE','#AA4499','grey')) +
  scale_alpha_manual(values = c(1,1,0.5))

#### Venn all vs no TB risk cov ####
venn.ls <- list()
venn.ls[["All covariates"]] <- dat %>% 
  filter(full < signif) %>% 
  pull(snpID)
venn.ls[["No TB risk covariates"]] <- dat %>% 
  filter(full_comp < signif) %>% 
  pull(snpID)

p3 <- ggvenn(venn.ls, show_percentage = FALSE,
             fill_color = c("white","white"),
             stroke_size = 0.5, set_name_size = 3, text_size = 3) +
  labs(x = paste("Significant P <",signif)) +
  theme(axis.title.x=element_text(size = 9)) +
  lims(y=c(-1,1.3))
# p3

#### Venn cov impacts ####
venn.ls2 <- list()

#Main models
venn.ls2[["All covariates"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value" & pval < signif) %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

venn.ls2[["No TB risk covariates"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_full_comp.csv") %>% 
  filter(term=="value" & pval < signif) %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

#Leave one out models
venn.ls2[["+ smoke + DM"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_loo.csv") %>% 
  filter(term=="value" & pval < signif & model == "minus_hiv") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

venn.ls2[["+ smoke + HIV"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_loo.csv") %>% 
  filter(term=="value" & pval < signif & model == "minus_diabetes") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

venn.ls2[["+ DM + HIV"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_loo.csv") %>% 
  filter(term=="value" & pval < signif & model == "minus_smoke") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

#Univariate models
venn.ls2[["+ HIV"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_univar.csv") %>% 
  filter(term=="value" & pval < signif & model == "plus_hiv") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

venn.ls2[["+ DM"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_univar.csv") %>% 
  filter(term=="value" & pval < signif & model == "plus_diabetes") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

venn.ls2[["+ smoke"]] <- read_csv("../ReporTB_lpWGS/result/model_final/model_univar.csv") %>% 
  filter(term=="value" & pval < signif & model == "plus_smoke") %>% 
  filter(snpID %in% maf5_concord) %>% 
  pull(snpID)

# v1 <- ggvenn(venn.ls2[c(2,6:8)], show_percentage = FALSE,
#              set_name_size=4) +
#   labs(x = paste("Signif P <", signif))+
#   theme(axis.title.x=element_text()) +
#   lims(y=c(-2,1.5),x=c(-3,2)) +
#   labs(title="Adding cov to\n~ 2PCs + kin + age + sex")

p4 <- ggvenn(venn.ls2[c(1,3:5)], show_percentage = FALSE,
             fill_color = c("white","white","white","white"),
             stroke_size = 0.5, set_name_size = 3, text_size = 3) +
  labs(x = paste("Significant P <", signif))+
  theme(axis.title.x=element_text(size = 9)) +
  lims(y=c(-2,1.5),x=c(-3,2.2))
p4b <- ggvenn(venn.ls2[c(3:5)], show_percentage = FALSE,
             fill_color = c("white","white","white","white"),
             stroke_size = 0.5, set_name_size = 3, text_size = 3) +
  labs(x = paste("Significant P <", signif))+
  theme(axis.title.x=element_text(size = 9)) +
  lims(y=c(-2,1.8),x=c(-3,2.2))
p4+p4b

#### Save ####
fig4 <- p1+p3 + plot_annotation(tag_levels = "A") 
# fig4
# ggsave("temp/Fig3_simple_model.pdf", fig4, width = 7, height = 5)
ggsave("temp/Fig3_simple_model.png", fig4, width = 7, height = 5)
