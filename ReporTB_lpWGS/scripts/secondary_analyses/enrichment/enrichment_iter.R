library(tidyverse)
library(SEARchways)
library(viridis)
library(patchwork)
library(openxlsx)
#

#### Protein-coding ####
signif <- read_csv("result/model_final_anno/ReporTB_signif_snp_anno.csv")

gene_pc <- signif %>% 
  select(snpID, contains("gene"), contains("biotype")) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  separate(name, into=c("loc","name","ID"), sep="_") %>% 
  pivot_wider() %>% 
  filter(biotype == "protein_coding") %>% 
  select(snpID, gene)
 
h_pc <- iterEnrich(anno_df = gene_pc, 
                   anno_featCol = "snpID",
                   anno_annotationCol = "gene",
                   niter = 100,
                   category = "H", protein_coding = TRUE,
                   ncores=3)

c2_pc <- iterEnrich(anno_df = gene_pc, 
                    anno_featCol = "snpID",
                    anno_annotationCol = "gene",
                    niter = 100,
                    category = "C2", subcategory = "CP",
                    protein_coding = TRUE,
                    ncores=3)

c5_pc <- iterEnrich(anno_df = gene_pc, 
                    anno_featCol = "snpID",
                    anno_annotationCol = "gene",
                    niter = 100,
                    category = "C5", subcategory = "GO:BP",
                    protein_coding = TRUE,
                    ncores=3)

enrich_pc <- h_pc$summary %>%
  mutate(gs_cat="H") %>% 
  bind_rows(c2_pc$summary %>% mutate(gs_cat="C2", gs_subcat="CP")) %>% 
  bind_rows(c5_pc$summary %>% mutate(gs_cat="C5", gs_subcat="GO:BP"))

#### Save ####
save(enrich_pc, 
     file="result/enrich/enrich_iter_signif_concord_all.RData")
save(h_pc, c2_pc, c5_pc, 
     file="result/enrich/enrich_iter_signif_concord.RData")

#### plots ####

p1 <- enrich_pc %>% 
  filter(FDR < 0.3 & k_median > 1) %>% 
  mutate(fdr.group = case_when(FDR<0.2~"FDR<0.2",
                               TRUE~"FDR<0.3"),
         fdr.group = factor(fdr.group, levels=c("FDR<0.3","FDR<0.2"))) %>% 
  ggplot(aes(x=reorder(pathway, -log10(FDR)), 
             y="Significant protein-coding genes", 
             size=fdr.group, color=`k/K_median`*100)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  facet_grid(gs_cat~., scales="free", space="free") +
  scale_color_viridis() +
  labs(y="", color="Percent enriched", x="Pathway",
       title="CP and GO\nFDR < 0.3\noverlap > 1", size="")
# p1

ggsave(p1, file="result/enrich/enrich_iter_gene_pc_fdr0.3_overlap2.png", 
       width=7, height=3)

