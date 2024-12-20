library(tidyverse)
library(limma)
# library(Seurat)
library(ggrepel)
library(ggpp)
library(ggtext)
# library(SNPRelate)
library(patchwork)

#### Data ####
# Model results
attach("../ReporTB_sceQTL/result/sceQTL/sceQTL_df.RData")
model_noCov <- model_noCov %>% 
  filter(variable_group=="SNP" & condition=="media")
signif <- 0.05

# Nice cell names
cell_key <- read_csv("../ReporTB_sceQTL/result/sc_key.csv") %>% 
  mutate(anno_long = gsub("cluster ","",anno_long))

#Format results
model_snp <- model_noCov %>% 
  mutate(anno_short = gsub("_"," ",cell)) %>% 
  left_join(cell_key, by="anno_short") %>% 
  mutate(col.group = ifelse(FDR < signif, anno_short, "NS"),
         col.group = fct_relevel(factor(col.group), "NS", 
                                 after = Inf)) %>% 
  mutate(lab = ifelse(FDR < signif,"Y","N"))

#Total signif combos
model_snp %>% 
  filter(FDR < signif) %>% 
  nrow()
model_snp %>% 
  filter(FDR < signif) %>%
  distinct(genotype, gene) %>% 
  nrow()

model_snp %>% 
  filter(FDR < signif) %>% 
  pull(genotype) %>% unique() %>% length()
model_snp %>% 
  filter(FDR < signif) %>% 
  pull(gene) %>% unique() %>% length()

model_snp2 <- model_snp %>% 
  #collapse multi annos
  group_by(cell, gene, estimate, FDR, col.group, lab) %>% 
  mutate(reps = n(),
         gene_snps = ifelse(
           reps>1, 
           paste0(unique(gene)," [",reps,"]"), 
           unique(gene))) %>% 
  ungroup() %>% 
  distinct(cell, gene, gene_snps, estimate, FDR, col.group, lab)

#### Volcano plot ####
pal <- c("#332288","#117733","#88CCEE","#DDCC77","#CC6677")
col.vec <- c(pal, pal, pal[1:4], "grey80")
shp.vec <- c(rep(24,5), rep(22,5), rep(23,4), 21)
title <- paste0("Cell FDR < ",signif)

p1 <- model_snp2  %>% 
  droplevels() %>% 
  ggplot(aes(x=estimate,y=-log10(FDR))) +
  geom_hline(yintercept = -log10(signif), 
             lty="dashed", color="grey") +
  # geom_vline(xintercept=0, lty="dashed") +
  geom_jitter(data = model_snp2 %>% filter(col.group=="NS"),
              shape=20, color="grey80")+
  geom_jitter(data = model_snp2 %>% filter(col.group!="NS"),
              aes(shape=col.group,
                  fill=col.group), 
              size=2, width=0, height=0, color="black", stroke=0.2)+
  theme_bw(base_size = 8) +
  labs(color=title, shape=title, fill=title,
       x="Log2 estimate",
       y="-log<sub>10</sub>( FDR )") +
  geom_text_repel(data=model_snp2 %>% 
                    filter(lab=="Y" & estimate<0),
                  aes(label = gene_snps), 
                  position = position_nudge_to(x = -3),
                  min.segment.length = 0, size=3,
                  hjust="right", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf,
                  segment.size=0.2) +
  geom_text_repel(data=model_snp2 %>%
                    filter(lab=="Y" & estimate>0),
                  aes(label = gene_snps), 
                  position = position_nudge_to(x = 2.5),
                  min.segment.length = 0, size=3,
                  hjust="left", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf,
                  segment.size=0.2) +
  scale_fill_manual(values=col.vec) +
  scale_color_manual(values=col.vec) +
  scale_shape_manual(values=shp.vec) +
  theme(strip.background = element_rect(fill="white"),
        legend.position = "right",
        legend.spacing.y = unit(-0.3, 'lines'),
        legend.key.spacing.y = unit(-0.3, 'lines'),
        axis.title.y = element_markdown(),
        legend.title = element_blank(),
        legend.text=element_text(size=8))+
  guides(fill=guide_legend(ncol=1, bycol = TRUE),
         shape=guide_legend(ncol=1, bycol = TRUE))

# p1

#### sceQTL ####
#SNP-gene map
map <- read_csv("../ReporTB_sceQTL/result/sceQTL/sceQTL_map.csv")
#SNP data
geno.num <- read_csv("../ReporTB_sceQTL/data_clean/geno.num.subset.csv")
#Single cell RNAseq data
attach("../ReporTB_sceQTL/data_clean/ReporTB_sc_voom.RData")

snpOI <- model_snp %>% 
  filter(FDR<signif) %>% 
  select(gene,cell, genotype)

to_plot <- model_noCov %>% 
  filter(variable_group=="SNP" & condition=="media") %>% 
  inner_join(snpOI)

#### Boxplots ####
to_plot_lvl <- to_plot %>% 
  distinct(cell, genotype, gene) %>% 
  arrange(gene, cell, genotype) %>% 
  # mutate(first="Y")
  # #label right in gene row
  full_join(data.frame(
    cell = c("CD14_Mono-2","CD14_Mono-2"),
    genotype = c("snp.3_75816280_A_T",
                "snp.4_139801643_G_T"),
    gene = c("ZNF717","MAML3"),
    first="Y"))

plot.ls <- list()

for(i in 1:nrow(to_plot_lvl)){
  temp <- to_plot_lvl[i,]
  print(paste(temp$cell,temp$genotype,temp$gene))
  snpSignif <- temp$genotype
  snpNice <- gsub("snp.","",snpSignif)
  snpNice <- gsub("_",":",snpNice)
  cellSignif <- temp$cell
  geneSignif <- temp$gene
  
  geno.temp <- geno.num %>% 
    select(sample_id, all_of(snpNice)) %>% 
    pivot_longer(-sample_id, names_to = "snpID", values_to = "geno")
  
  voom.temp <- get(paste0("voom_", cellSignif))
  e.temp <- as.data.frame(voom.temp$E) %>% 
    rownames_to_column("gene") %>% 
    filter(gene == geneSignif) %>%
    pivot_longer(-gene, names_to = "libID", values_to = "log2CPM") %>% 
    left_join(voom.temp$targets, by="libID") %>% 
    select(sample_id,gene:log2CPM, condition, anno_short) 
  
  fdr.temp <- model_noCov %>% 
    filter(variable_group=="SNP") %>% 
    filter(cell==cellSignif & genotype==snpSignif & gene ==geneSignif) %>% 
    mutate(FDR = signif(FDR, digits=2),
           pval = signif(pval, digits=2)) %>% 
    mutate(fc= ifelse(estimate<0,"down","up")) %>% 
    distinct(genotype, cell, gene, pval, FDR, condition, fc) %>% 
    mutate(snpID = gsub("^snp.","",genotype),
           snpID = gsub("_",":",snpID))  %>% 
    select(snpID,gene,cell, condition, FDR, fc) %>% 
    mutate(anno_short = gsub("_"," ",cell))
  
  plot_dat <- geno.temp %>%
    full_join(e.temp, by="sample_id") %>% 
    left_join(fdr.temp, by=c("snpID","gene","anno_short")) %>% 
    mutate(facet_lab = paste0(anno_short, 
                              "\nFDR = ",FDR))%>% 
    mutate(geno_fct = as.factor(geno)) %>% 
    drop_na(log2CPM)
  
  #merge geno 1 and 2 if dominant model
  #Does not occur for any signif models
  if(any(geno.temp$geno==2)){
    print("CHECK 3 LEVELS")
  #   plot_dat <- plot_dat %>%
  #     mutate(geno_fct = fct_recode(geno_fct, "1 or 2"="1", "1 or 2"="2"))
  }
  
  cellNice <- unique(fdr.temp$anno_short)
  
  line.col <- ifelse(unique(plot_dat$fc)=="down","blue","red")
  
  base_plot <- plot_dat %>% 
    ggplot(aes(x=geno_fct, y=log2CPM, color=geno_fct)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    geom_smooth(method="lm", aes(x=as.numeric(geno_fct)), 
                formula = 'y ~ x', se=FALSE, color=line.col) +
    theme_bw(base_size = 8) +
    scale_color_manual(values=c("grey70","black")) +
    facet_wrap(~facet_lab)
  
  if(is.na(temp$first)){
    plot.ls[[i]] <- base_plot +
      labs(x=snpNice, y="") +
      theme(strip.background=element_rect(fill="white"), 
            legend.position = "none",
            text=element_text(size=8),
            axis.title.x = element_text(size = 8),
            strip.text.x = element_text(size = 8),
            axis.title.y=element_blank())
  } else{
    plot.ls[[i]] <- base_plot+
      theme(strip.background=element_rect(fill="white"), 
            legend.position = "none",
            text=element_text(size=8),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            strip.text.x = element_text(size = 8)) +
      labs(x=snpNice, y="Normalized log2 CPM")
  }
}

#### Plot all ####
p2 <- 
  #ZNF (6-13) 
  plot.ls[[8]]+plot.ls[[12]]+plot.ls[[9]]+plot.ls[[13]]+plot_spacer()+
  #MAML3
  plot.ls[[3]]+plot.ls[[4]]+plot.ls[[5]]+
  #COX10
  plot.ls[[2]]+
  #ATF4
  plot.ls[[1]] +
  plot_layout(ncol=5) 
# p2

#### Not in main plot ####
#Not lead ZNF SNP
#plot.ls[[6]]+plot.ls[[7]]+
#plot.ls[[10]]+plot.ls[[11]]+

#### Save ######
p <- free(p1)/p2+ 
  plot_layout(heights = c(1,1.2))+
  plot_annotation(tag_levels = list(c("A","B ZNF717","","","",
                                      "C MAML3","","",
                                      "D COX10", "E ATF4"
                                      ))) & 
  theme(plot.tag.position = c(0.3, 1.075),
        plot.margin = margin(10, 5, 5, 5))
# p

ggsave(p, filename="Fig4_sceQTL_media.png",
       width=5.1, height=6)

#### Summary numbers ####
model_snp %>% distinct(genotype,gene) %>% nrow()
model_snp %>% distinct(genotype) %>% nrow()
model_snp %>% distinct(gene) %>% nrow()
# model_snp %>% distinct(condition) %>% nrow()
model_snp %>% distinct(anno_long) %>% nrow()

model_snp %>% filter(FDR<signif)
model_snp %>% filter(FDR<signif) %>% nrow()
model_snp %>% filter(FDR<signif) %>% count(genotype)
model_snp %>% filter(FDR<signif) %>% count(gene)
# model_snp %>% filter(FDR<signif) %>% count(condition)
model_snp %>% filter(FDR<signif) %>% count(anno_long)

