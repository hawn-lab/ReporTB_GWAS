library(tidyverse)
library(viridis)
library(patchwork)

signif <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv") %>% 
  pull(snpID)

#### LD ZNF ####
load("../ReporTB_lpWGS/result/LD/3_75816280_A_T.all.RData")

LD.df.znf <- LD.df.all %>% 
  # select(-Dprime) %>% 
  separate(V1, into=c("chr1","pos1"), sep=":", remove = FALSE, extra = "drop") %>% 
  separate(V2, into=c("chr2","pos2"), sep=":", remove = FALSE, extra = "drop") 

#Entire region
p1 <- LD.df.znf %>% 
  mutate(across(c(pos1,pos2), ~ifelse(.%in%c("75816280","75816219",
                                             "75816211","75816310"),
                                      paste0(.,"**"),.))) %>% 
  mutate(R2=as.numeric(R2)) %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  ggplot(aes(pos1, pos2, fill= R2)) + 
  geom_tile() +
  geom_hline(yintercept = c(29.5,34.5), color="white",size=0.25) +
  geom_vline(xintercept = c(30.5,35.5), color="white",size=0.25) +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="", title="ZNF717") +
  theme_void(base_size = 8) +
  coord_fixed() +
  theme(legend.position = "top", legend.title.position = "top",
        legend.title.align = 0.5) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5))
# p1

#Signif SNPs
p2 <- LD.df.znf %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  filter(V1 %in% signif & V2 %in% signif) %>% 
  mutate(across(c(V1,V2), ~gsub("3:75816","",.))) %>% 
  
  ggplot(aes(V1, V2, fill= R2)) + 
  geom_tile() +
  scale_fill_viridis() +
  labs(x="",y="") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        legend.margin = margin(-5, 0, -5, 0), # trbl
        legend.box.margin = margin(-5, 0, -5, 0) ) +
  coord_fixed() +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 4))
# p2
nudge<-0.08
p_ld <- p1+ 
  annotation_custom(grid::linesGrob(
    x = unit(c(0.578, 0.68), "npc"),
    y = unit(c(0.575, 0.4), "npc"),
    gp = grid::gpar(col = "black",lwd = 1.5)))+
  annotation_custom(grid::linesGrob(
    x = unit(c(0.655, 0.93), "npc"),
    y = unit(c(0.652, 0.43), "npc"),
    gp = grid::gpar(col = "black",lwd = 1.5)))+
  inset_element(p2, 0.38,-0.6,1.03,1)  # trbl
#
#### LD other genes ####
load("../ReporTB_lpWGS/result/LD/22_39566026_A_T.all.RData")

LD.df.atf <- LD.df.all %>% 
  separate(V1, into=c("chr1","pos1"), sep=":", remove = FALSE, extra = "drop") %>% 
  separate(V2, into=c("chr2","pos2"), sep=":", remove = FALSE, extra = "drop") %>% 
  mutate(across(c(pos1,pos2), ~ifelse(.=="39566026",
                                      paste0(.,"**"),.))) 

#Entire region
p2a <- LD.df.atf %>% 
  mutate(R2=as.numeric(R2)) %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  ggplot(aes(pos1, pos2, fill= R2)) + 
  geom_tile() +
  geom_hline(yintercept = c(18.5,19.5), color="white",size=0.25) +
  geom_vline(xintercept = c(19.5,20.5), color="white",size=0.25) +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="", title="ATF4") +
  theme_void(base_size = 8) +
  coord_fixed() +
  theme(legend.position = "bottom", legend.title.position = "top",
        legend.title.align = 0.5) +
  guides(fill = guide_colorbar(barwidth = 8, barheight = 0.5))
p2a

#####
load("../ReporTB_lpWGS/result/LD/17_14155330_A_G.all.RData")

LD.df.cox <- LD.df.all %>% 
  separate(V1, into=c("chr1","pos1"), sep=":", remove = FALSE, extra = "drop") %>% 
  separate(V2, into=c("chr2","pos2"), sep=":", remove = FALSE, extra = "drop") %>% 
  mutate(across(c(pos1,pos2), ~ifelse(.=="14155330",
                                      paste0(.,"**"),.))) 

#Entire region
p2b <- LD.df.cox %>% 
  mutate(R2=as.numeric(R2)) %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  ggplot(aes(pos1, pos2, fill= R2)) + 
  geom_tile() +
  geom_hline(yintercept = c(27.5,28.5), color="white",size=0.25) +
  geom_vline(xintercept = c(28.5,29.5), color="white",size=0.25) +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="", title="COX10") +
  theme_void(base_size = 8) +
  coord_fixed() +
  theme(legend.position = "none")
p2b

#####
load("../ReporTB_lpWGS/result/LD/4_139801643_G_T.all.RData")

LD.df.mam <- LD.df.all %>% 
  separate(V1, into=c("chr1","pos1"), sep=":", remove = FALSE, extra = "drop") %>% 
  separate(V2, into=c("chr2","pos2"), sep=":", remove = FALSE, extra = "drop") %>% 
  mutate(across(c(pos1,pos2), ~ifelse(.=="139801643",
                                      paste0(.,"**"),.)))

#Entire region
p2c <- LD.df.mam  %>% 
  mutate(R2=as.numeric(R2)) %>% 
  arrange(V1,V2) %>% 
  ggplot(aes(pos1, pos2, fill= R2)) + 
  geom_tile() +
  geom_hline(yintercept = c(14.5,15.5), color="white",size=0.25) +
  geom_vline(xintercept = c(15.5,16.5), color="white",size=0.25) +
  scale_fill_viridis(limits=c(0,1)) +
  labs(x="",y="", title="MAML3") +
  theme_void(base_size = 8) +
  coord_fixed() +
  theme(legend.position = "none")
p2c
p2a+p2b+p2c
#
#### Co-localization ####
# GWAS results
maf1_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf1_concord.txt",
                         col_names=FALSE) %>% pull(X1)

#Full cohort GWAS
full <- read_csv("../ReporTB_lpWGS/result/model_final/model_full.csv") %>% 
  filter(term=="value") %>% 
  select(snpID, pval)
full2 <- read_tsv("../ReporTB_lpWGS/result/model_final/glmm_score_full.txt") %>%
  filter(! SNP %in% full$snpID) %>% 
  rename(snpID=SNP,pval=PVAL) %>% 
  select(snpID, pval)

full_all <- bind_rows(full, full2) %>% 
  #MAF > 1%
  filter(snpID %in% maf1_concord) %>% 
  #positions
  separate(snpID, into = c("CHR","POS","REF","ALT"), sep=":", remove=FALSE) %>% 
  mutate(POS = as.numeric(POS))

#eQTL results
load("../ReporTB_sceQTL/result/sceQTL/colocal.RData")

colocal_noCov <- colocal_noCov %>% 
  filter(variable_group=="SNP")

dat1 <- full_all %>% 
  select(CHR, POS, pval) %>% 
  rename(GWAS=pval)

dat2 <- colocal_noCov %>% 
  # filter(gene == "ZNF717") %>% 
  mutate(annotation=gsub("_"," ",cell)) %>% 
  separate(variable, into=c("trash","CHR","POS"), sep="[.]|_", 
           remove=FALSE, extra="drop") %>% 
  mutate(POS = as.numeric(POS)) %>% 
  select(annotation, gene, CHR, POS, pval) %>% 
  pivot_wider(names_from = annotation, values_from = pval)

dat <- inner_join(dat1,dat2) %>% 
  pivot_longer(-c(gene,CHR,POS)) %>% 
  mutate(group = ifelse(name=="GWAS","GWAS Full cohort",
                        "sceQTL")) %>% 
  mutate(cutoff = ifelse(name=="GWAS",5E-8,0.05)) %>% 
  drop_na(value) %>% 
  mutate(name = fct_relevel(factor(name), "GWAS", after=0))

plot.ls2 <- list()
for(g in sort(unique(dat$gene))){
  dat.temp <- dat %>% filter(gene==g)
  xlab <- paste("Chromosome",unique(dat.temp$CHR))
  xbreak <- dat.temp %>% 
    slice_min(value) %>% pull(POS)
  
  plot.ls2[[g]] <- dat.temp %>% 
    
    ggplot(aes(x=POS, y=-log10(value))) +
    geom_point(aes(color=name), alpha=0.5) +
    geom_smooth(aes(color=name), se=FALSE, method="loess", 
                formula = 'y ~ x')+
    labs(y="-log10( P )", color="", x=xlab) +
    theme_bw() +
    facet_wrap(~group, scales="free", ncol=4)+ 
    geom_hline(aes(yintercept = -log10(cutoff)), lty="dashed") +
    geom_hline(aes(yintercept = -log10(cutoff)), lty="dashed") +
    scale_color_manual(values=c("GWAS"="grey30",
                                "CD14 Mono-2"="#E66100",
                                "CD4 CM"="#56B4E9",
                                "CD4 naive-3"="#A8016B",
                                "CD4 Treg"="#DDCC77",
                                "cDC"="#0072B2",
                                "MAIT"="#CC79A7",
                                "NK-1"="#E69F00",
                                "NK-2"="#009E73"),
                       breaks=c("CD14 Mono-2","CD4 CM",
                                "CD4 naive-3","CD4 Treg",
                                "cDC","MAIT","NK-1","NK-2")) +
    theme(legend.position = "bottom", legend.direction="horizontal",
          legend.margin = margin(-5, 0, -5, 0), # trbl
          legend.box.margin = margin(-5, 0, -5, 0) ) +
    guides(color = guide_legend(nrow = 1)) +
    scale_x_continuous(
      breaks = c(xbreak-7500,xbreak,xbreak+7500)
    )
}

# wrap_plots(plot.ls2)

### save ####
lo <-"
AEE
BFF
CGG
DHH
"
p_all <- p_ld+p2c+p2b+p2a + 
  free(plot.ls2[[4]])+free(plot.ls2[[3]])+
  free(plot.ls2[[2]])+free(plot.ls2[[1]])+
  plot_annotation(tag_levels = list(c("A","","B","C","D"))) +
  plot_layout(design = lo, heights = 1)
# p_all

ggsave(p_all, filename="FigS8_LD.png", 
       width=7.5, height=10)

