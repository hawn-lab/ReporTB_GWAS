library(tidyverse)
library(patchwork)

#### PC percent ####
load("../ReporTB_sceQTL/result/kinship/reportTB_sc_PCair.RData")

#Plot PC contributions
p1 <- data.frame(PC = 1:length(pc.percent),
                 perc = pc.percent) %>% 
  # mutate(PC=factor(PC, levels=paste0("PC",1:length(pc.percent)))) %>% 
  ggplot(aes(x=PC, y=perc)) + 
  geom_bar(stat = "identity") +
  labs(x="PC",y="Percent variation explained") +
  theme_classic()
# p1

#### Kinship ####
kin_raw <- read_csv("../ReporTB_sceQTL/result/kinship/reportTB_sc_kin.csv") %>% 
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

p3 <- kin_long %>% 
  ggplot(aes(x=1, y=value)) +
  geom_violin() +
  geom_jitter(data = filter(kin_long, value > 0.125), 
              aes(color = col_group)) +
  geom_hline(yintercept = c(0.125, 0.25, 0.5), 
             lty = "dashed", color = "red") +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Kinship", color = "Relatedness") +
  scale_color_manual(values=c("black","grey35","grey70")) +
  theme(legend.position = "left")
# p3

# length(unique(c(kin_long_filter$rowname, kin_long_filter$name)))
table(kin_long$col_group)

#### Covariate fit ####
output3 <- read_csv("../ReporTB_sceQTL/result/sceQTL/eqtl_model_fitting.csv") %>% 
  distinct()

#calculate delta AIC for kinship
output_delta3 <- output3 %>% 
  select(model, cell, gene, genotype, AIC) %>% 
  pivot_wider(names_from = model, values_from = AIC) %>% 
  pivot_longer(-c(genotype,cell, gene,base,full)) %>% 
  mutate(delta = case_when(!grepl("minus", name)~value-base,
                           grepl("minus", name)~value-full)) %>% 
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
  labs(x="GLMMkin + 2 PCs + sex + age + smoking + diabetes", 
       y="Change in AIC\n Better fit\nwith covariate removed <--   --> with covariate retained") +
  geom_hline(yintercept = 0, color="red", lty="dashed") +
  theme_classic() +
  facet_wrap(~group, scales = "free")

#### Save ####
lo <- "
AAAAB
CCCDD
"
plot_all <- p1+p3+p5+p6 + plot_annotation(tag_levels = "A") +
  plot_layout(design=lo)
plot_all

ggsave("FigS3_sceQTL.PCs.kin.png", plot_all,
       width=10, height=8)
