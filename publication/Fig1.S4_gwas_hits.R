library(tidyverse)
# library(data.table)
library(ggtext)
library(ggvenn)
# library(ggrepel)
library(patchwork)
# library(viridis)
signif <- 5E-8
sugg <- 5E-5 

#### GWAS DATA ####
imp_maf5 <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB_filter_imp_maf5.bim",
                 col_names = FALSE) %>% pull(X2)
imp_maf5_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf5_concord.txt",
                 col_names=FALSE) %>% pull(X1)

#Full cohort
full <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value") %>% 
  select(snpID, pval)
full2 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_full.txt") %>%
  filter(! SNP %in% full$snpID) %>% 
  rename(snpID=SNP,pval=PVAL) %>% 
  select(snpID, pval)

full_all <- bind_rows(full, full2) %>% 
  #MAF > 5%
  filter(snpID %in% imp_maf5) %>% 
  #positions
  separate(snpID, into = c("CHR","POS","REF","ALT"), sep=":", remove=FALSE)

# Case-control
cc <- read_csv("../ReporTB_lpWGS/result/model_final/model_cc.csv") %>% 
  filter(term=="value") %>% 
  select(snpID, pval)
cc2 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_cc.txt") %>% 
  filter(! SNP %in% cc$snpID) %>% 
  rename(snpID=SNP,pval=PVAL) %>% 
  select(snpID, pval)

cc_all <- bind_rows(cc, cc2) %>% 
  #MAF > 5%
  filter(snpID %in% imp_maf5) %>% 
  #positions
  separate(snpID, into = c("CHR","POS","REF","ALT"), sep=":", remove=FALSE)

rm(full, full2, cc, cc2)
gc()

#### QQ PLOT ####
##### full cohort ####
full.signif <- full_all %>% filter(snpID %in% imp_maf5_concord & pval < signif) 
full.sugg <- full_all %>% filter(snpID %in% imp_maf5_concord & pval < sugg) 
full.obs <- full_all %>% filter(snpID %in% imp_maf5_concord) %>% pull(pval)
full.exp <- (rank(full.obs, ties.method="first")+0.5)/(length(full.obs)+1)
full.snp <- full_all %>% filter(snpID %in% imp_maf5_concord) %>% pull(snpID)

p0 <- data.frame(snpID = full.snp,
                 exp = -log10(full.exp),
                 obs = -log10(full.obs)) %>% 
  mutate(col.group = case_when(snpID %in% full.signif$snpID ~ "GLMMkin P < 5E-8",
                               snpID %in% full.sugg$snpID ~ "GLMMkin P < 5E-5",
                               TRUE ~ "None")) %>% 
  mutate(col.group = factor(col.group, 
                            levels=c("GLMMkin P < 5E-8",
                                     "GLMMkin P < 5E-5",
                                     "None"))) %>% 
  arrange(desc(col.group)) %>%
  
  ggplot() +
  aes(x=exp, y=obs) +
  geom_point(aes(color = col.group)) + 
  theme_classic(base_size = 10) +
  labs(x = "Expected -log10( pval )",
       y = "Observed -log10( pval )",
       color = "Significance", title="Full cohort") +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = c("black","grey50","grey90"), drop=FALSE) +
  coord_fixed() 

##### cc cohort ####
cc.signif <- cc_all %>% filter(snpID %in% imp_maf5_concord & pval < signif)
cc.sugg <- cc_all %>% filter(snpID %in% imp_maf5_concord & pval < sugg)
cc.obs <- cc_all %>% filter(snpID %in% imp_maf5_concord) %>% pull(pval)
cc.exp <- (rank(cc.obs, ties.method="first")+0.5)/(length(cc.obs)+1)
cc.snp <- cc_all %>% filter(snpID %in% imp_maf5_concord) %>% pull(snpID)

p1 <- data.frame(snpID = cc.snp,
                 exp = -log10(cc.exp),
                 obs = -log10(cc.obs)) %>% 
  mutate(col.group = case_when(snpID %in% cc.signif$snpID ~ "GLMMkin P < 5E-8",
                               snpID %in% cc.sugg$snpID ~ "GLMMkin P < 5E-5",
                               TRUE ~ "None")) %>% 
  mutate(col.group = factor(col.group, 
                            levels=c("GLMMkin P < 5E-8",
                                     "GLMMkin P < 5E-5",
                                     "None"))) %>% 
  arrange(desc(col.group)) %>% 
  ggplot() +
  aes(x=exp, y=obs) +
  geom_point(aes(color = col.group)) + 
  theme_classic(base_size = 10) +
  labs(x = "Expected -log10( pval )",
       y = "Observed -log10( pval )",
       color = "Significance", title="Case-control cohort") +
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = c("black","grey50","grey90"), drop=FALSE) +
  coord_fixed()+
  theme(legend.position = "none")

##### Save QQ ####
ggsave(p1+p0, filename="FigS4_QQplot.png", height=3, width=6)
beepr::beep()

##### QQ lambda ####
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

inflation(full.obs) # 0.9853587
inflation(full.exp) # 0.9999996
inflation(cc.obs)   # 0.9940643
inflation(cc.exp)   # 0.9999996

#### VENN ####
##### Overlap full and cc models ####
#Before concordance cutoff
venn.ls0 <- list()
venn.ls0[["Full cohort\nP < 5e-08"]] <- full_all %>% filter(snpID %in% imp_maf5 & pval < signif) %>% pull(snpID)
venn.ls0[["Full cohort\nP < 5e-05"]] <- full_all %>% filter(snpID %in% imp_maf5 & pval < sugg) %>% pull(snpID)
#Zero don't include
# venn.ls0[["Case-control cohort\nP < 5e-08"]] <- cc_all %>% filter(snpID %in% imp_maf5 & pval < signif) %>% pull(snpID)
venn.ls0[["Case-control cohort\nP < 5e-05"]] <- cc_all %>% filter(snpID %in% imp_maf5 & pval < sugg) %>% pull(snpID)

p2a <- ggvenn(venn.ls0, show_percentage = FALSE,
             fill_color = c("pink","lightblue","lightblue"),
             stroke_size = 0.5, set_name_size = 3, text_size = 3) +
  labs(x = paste("Significant P <",signif), title="imp_maf5") +
  theme(axis.title.x=element_text(size = 9)) +
  lims(y=c(-2,2))

# After concordance cutoff
venn.ls <- list()
venn.ls[["Full cohort\nP < 5e-08"]] <- full_all %>% filter(snpID %in% imp_maf5_concord & pval < signif) %>% pull(snpID)
venn.ls[["Full cohort\nP < 5e-05"]] <- full_all %>% filter(snpID %in% imp_maf5_concord & pval < sugg) %>% pull(snpID)
#Zero don't include
# venn.ls[["Case-control cohort\nP < 5e-08"]] <- cc_all %>% filter(snpID %in% imp_maf5_concord & pval < signif) %>% pull(snpID)
venn.ls[["Case-control cohort\nP < 5e-05"]] <- cc_all %>% filter(snpID %in% imp_maf5_concord & pval < sugg) %>% pull(snpID)

# only_cc <- venn.ls[["Case-control cohort\nP < 5e-05"]]
# only_cc <- only_cc[!only_cc %in% venn.ls[["Full cohort\nP < 5e-05"]]]

p2b <- ggvenn(venn.ls, show_percentage = FALSE,
             fill_color = c("pink","lightblue","lightblue"),
             stroke_size = 0.5, set_name_size = 3, text_size = 3) +
  labs(x = paste("Significant P <",signif), title="imp_maf5_concord") +
  theme(axis.title.x=element_text(size = 9)) +
  lims(y=c(-2,2))
p2a+p2b

library(eulerr)
fit <- euler(c('A' = 0, 'B' = 238, 'A&B' = 19,
               'C' = 64, 'A&C' = 0, 'B&C' = 0,
               'A&B&C' = 0))
plot(fit)

#Compare pval full to cc
# full_all %>% 
#   filter(pval < signif) %>% 
#   rename(full_p=pval) %>% 
#   left_join(cc_all) %>% 
#   rename(cc_p=pval) %>% View

#### ANNOTATION #####
anno <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  filter(snpID %in% imp_maf5_concord)

any_anno <- anno %>% 
  select(snpID, contains("biotype")) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  pull(snpID) %>% unique()

no_anno <- anno %>% 
  filter(!snpID %in% any_anno) %>% 
  pull(snpID) %>% unique()

p3 <- anno %>% 
  select(snpID, contains("biotype")) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  separate(name, into=c("name"), sep="_", extra="drop") %>% 
  distinct() %>% 
  count(name, value) %>% 
  
  mutate(bio_group = case_when(
    grepl("pseudogene",value) ~"pseudogene",
    is.na(value)~"none",
    value=="protein_coding"~"protein-coding",
    grepl("RNA",value)~"RNA",
    grepl("IG_",value)~"IG",
    value=="TEC"~"pseudogene",
    TRUE ~ value)) %>%
  add_row(name="no annotation", value="none",n=length(no_anno), bio_group="none") %>% 
  mutate(bio_group = factor(bio_group, 
                            levels=rev(c("protein-coding","pseudogene",
                                         "RNA","IG","other","none")))) %>% 
  mutate(cis = fct_recode(name,
                          "cis up to 50kb"="cis"),
         cis = factor(cis, levels=c("intragenic",
                                    "cis up to 50kb",
                                    "no annotation"))) %>% 
  
  ggplot(aes(x=cis, y=n, fill=bio_group)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 10) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "left") +
  labs(y = "Total annotations", fill="Annotation") +
  scale_fill_manual(values = rev(c('#CC6677','#332288','#44AA99',
                                   'grey')))
# p3

#### MANHATTEN ####
dat_all <- full_all %>% 
  mutate(cohort="full") %>% 
  bind_rows(cc_all %>% mutate(cohort="cc")) %>% 
  filter(snpID %in% imp_maf5_concord)

dat <- dat_all %>% 
  mutate(POS = as.numeric(POS),
         CHR = as.numeric(CHR)) %>% 
  #Determine signif in one or both models
  pivot_wider(names_from = cohort, values_from = pval) %>% 
  mutate(group = case_when(
    full < signif & cc < signif ~ "Both significant",
    full < signif & cc < sugg ~ "Full cohort significant,\nCase-control suggestive",
    full < signif & cc > sugg ~ "Full cohort significant",
    full > sugg & cc < signif ~ "Case-control significant",
    
    full < sugg & cc < signif ~ "Full cohort suggestive,\nCase-control significant",
    full < sugg & cc < sugg ~ "Full cohort and\nCase-control suggestive",
    full < sugg & cc > sugg ~ "Full cohort suggestive",
    full > sugg & cc < sugg ~ "Case-control suggestive",
    # full > sugg & cc < sugg ~ "Case-control suggestive only",
    )) %>% 
  pivot_longer(full:cc, values_to = "pval") %>% 
  drop_na(pval) %>% 
  #chr colors
  mutate(chr.col = case_when(CHR %in% seq(1,22, 2) ~ "odd",
                             TRUE ~ "even")) %>% 
  mutate(col.final = case_when(!is.na(group) ~ group,
                               TRUE ~ chr.col),
         col.final = factor(
           col.final, 
           levels=rev(c("Both significant",
                        "Full cohort significant,\nCase-control suggestive",
                        "Full cohort significant",
                        "Case-control significant",
                        "Full cohort suggestive,\nCase-control significant",
                        "Full cohort and\nCase-control suggestive",
                        "Full cohort suggestive",
                        "Case-control suggestive",
                        "even","odd"))))%>% 
  arrange(col.final)

##### Manhattan plot - full ####
# X-axis cumulative
dat_cum <- dat %>% 
  #Full cohort only
  filter(name=="full") %>% 
  group_by(CHR) %>% 
  summarise(max_pos = max(POS)) %>% 
  mutate(pos_add = lag(cumsum(max_pos), default = 0)) %>% 
  select(CHR, pos_add)

dat_ID <- dat %>%
  #Full cohort only
  filter(name=="full") %>% 
  left_join(dat_cum) %>% 
  # slice_tail(n=100) %>%
  mutate(POS_cum = POS + pos_add)

# X-axis labels
axis_set <- dat_ID %>% 
  group_by(CHR) %>% 
  summarize(center = mean(POS_cum))

p4 <- dat_ID  %>% 
  arrange(col.final) %>% 
  ggplot(aes(x = POS_cum, y = -log10(pval))) +
  # geom_hline(yintercept = -log10(sugg), 
  #            color = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = -log10(signif), 
             color = "black", linetype = "dashed") +
  geom_point(data=dat_ID %>% 
               filter(col.final %in% c("odd","even")),
             alpha = 0.7, aes(color = col.final), 
             show.legend = FALSE) +
  geom_point(data= dat_ID %>% 
               filter(!col.final %in% c("odd","even")),
             aes(fill = col.final),pch=21, color="black",
             size=2) +
  #Label signif color points
  # geom_point(data = dat_ID2 %>% filter(!is.na(group)), 
  #            aes(color=group)) +
  # geom_text_repel(data = dat_ID2 %>% filter(!is.na(lab)),
  #                 aes(label = lab, color = group), #color = "black",
  #                 size=3, min.segment.length = unit(0.1, 'lines'),
  #                 show.legend = FALSE, max.overlaps = 100,
  #                 nudge_y = 0.3, direction = "y", hjust = "center") +
  #Beautify
  scale_x_continuous(label = axis_set$CHR,
                     breaks = axis_set$center) +
  scale_fill_manual(values = c("#88CCEE","#DDCC77","#1DB14F"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("grey80","grey40"))+ 
  
  labs(x = "Chromosome", 
       y = "-log<sub>10</sub>( p )",
       # fill = "Significance",
       title="Full cohort") + 
  theme_bw(base_size = 10) +
  theme( 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 10, vjust = 0.5),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.title = element_blank(), plot.title = element_text(size=10),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10)
  )
# p4

##### Manhattan plot - CC ####
# X-axis cumulative
dat_cum2 <- dat %>% 
  #cc cohort only
  filter(name=="cc") %>% 
  group_by(CHR) %>% 
  summarise(max_pos = max(POS)) %>% 
  mutate(pos_add = lag(cumsum(max_pos), default = 0)) %>% 
  select(CHR, pos_add)

dat_ID2 <- dat %>%
  #CC cohort only
  filter(name=="cc") %>% 
  left_join(dat_cum) %>% 
  # slice_tail(n=100) %>%
  mutate(POS_cum = POS + pos_add)

# X-axis labels
axis_set2 <- dat_ID2 %>% 
  group_by(CHR) %>% 
  summarize(center = mean(POS_cum))

p5 <- dat_ID2  %>% 
  arrange(col.final) %>% 
  ggplot(aes(x = POS_cum, y = -log10(pval))) +
  # geom_hline(yintercept = -log10(sugg), 
  #            color = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = -log10(signif), 
             color = "black", linetype = "dashed") +
  geom_point(data= dat_ID2 %>% 
               filter(col.final %in% c("odd","even")),
             alpha = 0.7, aes(color = col.final),
             show.legend = FALSE) +
  geom_point(data= dat_ID2 %>% 
               filter(!col.final %in% c("odd","even")),
             aes(fill = col.final), shape=21, color="black",
             size=2) +
  #Label signif color points
  # geom_point(data = dat_ID2 %>% filter(!is.na(group)), 
  #            aes(color=group)) +
  # geom_text_repel(data = dat_ID2 %>% filter(!is.na(lab)),
  #                 aes(label = lab, color = group), #color = "black",
  #                 size=3, min.segment.length = unit(0.1, 'lines'),
  #                 show.legend = FALSE, max.overlaps = 100,
  #                 nudge_y = 0.3, direction = "y", hjust = "center") +
  #Beautify
  scale_x_continuous(label = axis_set2$CHR,
                     breaks = axis_set2$center) +
  scale_fill_manual(values = c("#88CCEE","#DDCC77","#1DB14F"),
                    guide = guide_legend(reverse = TRUE)) +
  scale_color_manual(values = c("grey80","grey40"))+ 
  labs(x = "Chromosome", 
       y = "-log<sub>10</sub>( p )",
       # fill = "Significance",
       title="Case-control cohort") + 
  theme_bw(base_size = 10) +
  theme( 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    legend.position = "none", 
    axis.text.x = element_text(angle = 60, size = 10, vjust = 0.5), 
    legend.title=element_blank(),plot.title = element_text(size=10)
  )
# p5

#### Save ####

lo <- "
AAAABB
CCCCCC
DDDDDD
"

# plot_all <- p4/p5
plot_all <- plot_spacer()+p3+p4+p5+
  plot_layout(design=lo, heights=c(1,1,1)) + 
  plot_annotation(tag_levels = list(c("C","D","E")))
# plot_all

# ggsave("temp/Fig1_gwas_hits.pdf", plot_all, width=8, height=8)
ggsave("temp/Fig1_gwas_hits.png", plot_all, width=8, height=6)
beepr::beep()
