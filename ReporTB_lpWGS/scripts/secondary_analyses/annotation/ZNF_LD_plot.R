library(tidyverse)
library(viridis)
library(patchwork)

load("result/LD/3_75816280_A_T.all.RData")

signif <- read_csv("result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  pull(snpID)
# "3:75816211:A:G"   "3:75816219:T:G"   "3:75816280:A:T"
# "3:75816310:C:A"   "3:75816342:A:C"

LD.df.all %>% 
  select(V1, V2, R2) %>% 
  arrange(V1,V2) %>% 
  pivot_wider(names_from = V2, values_from = R2) %>% View

#Entire region
p1 <- LD.df.all %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  mutate(V1 = ifelse(V1=="3:75816280:A:T", paste0(V1,"*"), V1),
         V2 = ifelse(V2=="3:75816280:A:T", paste0(V2,"*"), V2)) %>% 
ggplot(aes(V1, V2, fill= R2)) + 
  geom_tile() +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed()

#Signif SNPs
p2 <- LD.df.all %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  filter(V1 %in% signif & V2 %in% signif) %>% 
ggplot(aes(V1, V2, fill= R2)) + 
  geom_tile() +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed()

p3 <- LD.df.all %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  filter(V1 %in% signif & V2 %in% signif) %>% 
  
  ggplot(aes(V1, V2, fill= R2)) + 
  geom_tile() +
  scale_fill_viridis() +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed()

LD.df.all %>% 
  mutate(R2=as.numeric(R2),
         Dprime=as.numeric(Dprime)) %>% 
  arrange(V1,V2) %>% 
  filter(V1 %in% signif & V2 %in% signif) %>% 
  summary()

p_all <- p1+(p2/p3) + plot_layout(widths = c(4,1))

ggsave(p_all, filename="result/LD/ZNF717_LD.png", 
       width=11, height=6)
