#Compare kinship values between controls that were calculated in full data set
#and the scRNAseq subset
library(tidyverse)
library(ggpubr)

#### Data ####
kin_sub <- read_csv("result/kinship/reportTB_sc_kin.csv")
kin_orig <- read_csv("../ReporTB_lpWGS/result/kinship/reportTB_kin.csv")

overlap <- intersect(kin_sub$rowname, kin_orig$rowname)

# Filter for overlap
kin_sub_long <- kin_sub %>% 
  filter(rowname %in% overlap) %>% 
  pivot_longer(-rowname, values_to = "kin_sub") %>% 
  filter(name %in% overlap)

kin_orig_long <- kin_orig %>% 
  filter(rowname %in% overlap) %>% 
  pivot_longer(-rowname, values_to = "kin_orig") %>% 
  filter(name %in% overlap)

dat <- full_join(kin_sub_long, kin_orig_long)
p1 <- dat %>% 
  # filter(rowname != name) %>%
  ggplot(aes(x=kin_sub, y=kin_orig)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept = 0, lty="dashed") +
  geom_smooth(method="lm", color = "grey") +
  stat_cor(digits=3) +
  lims(x=c(-0.04,1.2), y=c(-0.04,1.2)) +
  coord_fixed()
p1

ggsave("result/kinship/new_v_old_kinship.png", p1, width=3, height=3)
