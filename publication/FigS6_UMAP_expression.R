library(tidyverse)
library(Seurat)
library(ggrepel)
library(patchwork)

#### Single-cell Data ####
#Seurat
load("../ReporTB_sceQTL/data_clean/LTBI_scRNAseq_seurat_anno.RData")
sc_meta <- as.data.frame(sc@meta.data)

#expression
load("../ReporTB_sceQTL/data_clean/ReporTB_sc_voom.RData")

geneOI <- read_csv("../ReporTB_sceQTL/result/sceQTL/sceQTL_map.csv") %>% 
  pull(hgnc_symbol) %>%  unique() %>% sort()

#labels
cell_key <- read_csv("../ReporTB_sceQTL/result/sc_key.csv") %>% 
  mutate(anno_long = gsub("cluster ","",anno_long))

#### UMAP cell annotations ####
umap_anno <- sc@meta.data %>% 
  group_by(annotation) %>% 
  summarise(UMAP1=mean(UMAP1),
            UMAP2=mean(UMAP2)) %>% 
  ungroup() %>% 
  #nudge labels
  mutate(UMAP1 = case_when(annotation=="Platelet"~4,
                           # annotation=="Plasmablast"~2.5,
                           # annotation=="CD4 mem−2"~8.5,
                           TRUE~UMAP1)) %>%
  mutate(UMAP2 = case_when(annotation=="Platelet"~3.5,
                           # annotation=="CD4 mem−2"~8,
                           # annotation=="CD4 mem−3"~7,
                           TRUE~UMAP2)) %>%
  #add nice labels
  left_join(cell_key) %>% 
  mutate(ord = case_when(grepl("Mono",anno_short)~1,
                         grepl("DC",anno_short)~2,
                         grepl("HSPC",anno_short)~3,
                         anno_short%in%c("Platelet","Plasmablast")~4,
                         anno_short%in%c("B")~5,
                         grepl("NK",anno_short)~6,
                         grepl("CD8 CM|CD8 EM",anno_short)~7,
                         grepl("gdT|MAIT",anno_short)~8,
                         grepl("CD4 mem|CD4 CM|CD4 EM",anno_short)~9,
                         grepl("HSP|Treg|ISG|act",anno_short)~10,
                         grepl("naive",anno_short)~11,
                         TRUE~14),
         anno_long = fct_reorder(factor(anno_long),-ord),
         anno_short = fct_reorder(factor(anno_short),-ord))
#plot
p1 <- sc@meta.data %>% 
  filter(condition=="Media") %>% 
  left_join(umap_anno %>% select(-UMAP1,-UMAP2)) %>% 
  arrange(rev(annotation)) %>% 
  # slice_head(n=1000) %>% #For testing
  
  ggplot() +
  aes(x = UMAP1, y = UMAP2) +
  geom_point(aes(color = anno_short), size = 0.01) +
  labs(title = "", color = "") +
  theme_classic(base_size = 8) +
  # coord_fixed() +
  geom_label_repel(data=umap_anno, 
                   aes(label = anno_short, color=anno_short), 
                   direction = "both",
                   show.legend = FALSE, max.overlaps = 100,
                   size=3, fill="white", box.padding = 0.4) +
  # scale_color_manual(values=c(scales::hue_pal()(21),"grey")) +
  scale_color_manual(values=c(
    "#332288","#117733","#44AA99","#88CCEE","#DDCC77",
    "#CC6677","#AA4499","#882255","black",
    "#332288","#117733","#44AA99","#88CCEE","#DDCC77",
    "#CC6677","#AA4499","#882255","black",
    "#332288","#117733","#44AA99","#88CCEE","#DDCC77",
    "#CC6677","#AA4499","#882255","black",
    "#332288","#117733","#44AA99","#88CCEE","#DDCC77",
    "#CC6677","#AA4499","#882255")) +
  theme(legend.position = "none")#+
  # guides(color=guide_legend(nrow=18,override.aes = list(size=4)))
# p1

#### Expression ####
#missing gene
geneAll <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  select(snpID,contains("_gene_"), contains("_biotype_")) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  separate(name, into=c("anno","name","ID"), sep="_") %>% 
  pivot_wider() %>% 
  filter(biotype=="protein_coding") %>% 
  pull(gene) %>% unique()
geneAll[!geneAll %in% geneOI]

dat <- data.frame()
for(v in ls(pattern="voom_")){
  print(v)
  temp <- get(v)
  
  dat <- as.data.frame(temp$E) %>% 
    rownames_to_column("hgnc_symbol") %>% 
    filter(hgnc_symbol %in% geneAll) %>% 
    pivot_longer(-hgnc_symbol, names_to = "libID") %>% 
    left_join(temp$targets %>% select(libID, bid, condition, annotation, sample_id),
              by="libID") %>% 
    bind_rows(dat)
}

dat_summ <- dat %>%
  group_by(hgnc_symbol, condition, annotation) %>% 
  summarise(mean=mean(value,na.rm = TRUE),
            sd=sd(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  # mutate(annotation = gsub("-","_",annotation)) %>% 
  left_join(cell_key) %>% 
  mutate(col.group = case_when(mean<2.5~"< 2.5",
                               mean<5~"< 5.0",
                               mean<7.5~"< 7.5",
                               mean<10~"< 10"),
         col.group = factor(col.group, levels=c("< 2.5","< 5.0","< 7.5","< 10"))) %>% 
  mutate(ord = case_when(grepl("Mono",anno_short)~1,
                         grepl("DC",anno_short)~2,
                         grepl("HSPC",anno_short)~3,
                         anno_short%in%c("Platelet","Plasmablast")~4,
                         anno_short%in%c("B")~5,
                         grepl("NK",anno_short)~6,
                         grepl("CD8 CM|CD8 EM",anno_short)~7,
                         grepl("gdT|MAIT",anno_short)~8,
                         grepl("CD4 mem|CD4 CM|CD4 EM",anno_short)~9,
                         grepl("HSP|Treg|ISG|act",anno_short)~10,
                         grepl("naive",anno_short)~11,
                         TRUE~14),
         fct = case_when(grepl("Mono",anno_short)~"Monocyte",
                         grepl("DC",anno_short)~"DC",
                         grepl("HSPC",anno_short)~"Other",
                         anno_short%in%c("Platelet","Plasmablast")~"Other",
                         anno_short%in%c("B")~"Other",
                         grepl("NK",anno_short)~"NK",
                         grepl("CD8 CM|CD8 EM",anno_short)~"CD8+\nmemory",
                         grepl("gdT|MAIT",anno_short)~"T-cell\nother",
                         grepl("CD4 mem|CD4 CM|CD4 EM",anno_short)~"CD4+ memory",
                         grepl("HSP|Treg|ISG|act",anno_short)~"CD4+ other",
                         grepl("naive",anno_short)~"CD4+ or CD8+\nnaive",
                         TRUE~"Unclass"),
         anno_long = fct_reorder(factor(anno_long),-ord),
         fct = fct_reorder(factor(fct),ord)) %>% 
  mutate(hgnc_symbol=factor(hgnc_symbol))

p2a<- dat_summ %>%
  filter(ord <=7) %>%
  ggplot(aes(x=hgnc_symbol, y=anno_long)) +
  geom_point(aes(size=mean,color=mean)) +
  theme_bw(base_size=9) +
  facet_grid(fct~., scales="free_y", space="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None",
        strip.background = element_rect("white"),
        plot.margin = unit(c(5.5,-2,5.5,-2), "pt") #trbl
        ) +
  labs(x="",y="", size="Mean log2 CPM", color="Mean log2 CPM") +
  viridis::scale_color_viridis(option="magma",limits=c(0.4,10), 
                               breaks=c(2.5,5,7.5,10)) +
  scale_size_continuous(limits=c(0.1,10), range = c(0.1,4)) +
  scale_y_discrete(limits=rev) +
  guides(size = guide_legend(reverse=TRUE)) +
  scale_x_discrete(drop = FALSE)
# p2a

p2b<- dat_summ %>% 
  filter(ord > 7) %>%
  ggplot(aes(x=hgnc_symbol, y=anno_long)) +
  geom_point(aes(size=mean,color=mean)) +
  theme_bw(base_size=9) +
  facet_grid(fct~., scales="free_y", space="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box = "vertical",
        legend.position = "bottom",
        strip.background = element_rect("white"),
        plot.margin = unit(c(5.5,-2,5.5,-2), "pt") #trbl
        ) +
  labs(x="",y="", size="Mean log2 CPM", color="Mean log2 CPM") +
  viridis::scale_color_viridis(option="magma",limits=c(0.4,10), 
                               breaks=c(2.5,5,7.5,10)) +
  scale_size_continuous(limits=c(0.1,10), range = c(0.1,4)) +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(drop = FALSE)
  # guides(size = guide_legend(reverse=TRUE))
# p2a+p2b

#### Save ####
p_all <- free(p1)/free(p2a+p2b) +
  plot_layout(heights=c(1,1.3))+
  plot_annotation(tag_levels = list(c("A","B","")))
# p_all
ggsave(p_all,
       filename="temp/FigS6_UMAP_expression.pdf",
       width=7.5, height=9.5)

