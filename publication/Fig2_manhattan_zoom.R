library(tidyverse)
library(ggrepel)
library(viridis)
library(ggtext)
library(patchwork)

wrapper <- function(x, ...) { 
  paste(strwrap(x, ...), collapse = "\n")
}

#### Zoomed Manhattan with LD ####
maf1_concord <- read_tsv("../ReporTB_lpWGS/result/snp_imp_maf1_concord.txt",
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
  #MAF > 1%
  filter(snpID %in% maf1_concord) %>% 
  #positions
  separate(snpID, into = c("CHR","POS","REF","ALT"), sep=":", remove=FALSE)

# Signif protein-coding snps
signif <- read_csv("../ReporTB_lpWGS/result/model_final_anno/ReporTB_signif_snp_anno.csv")

biotypes <- colnames(signif)[grepl("biotype",colnames(signif))]

signif_snp <- signif %>% 
  filter(if_any(biotypes, ~.=="protein_coding")) %>% 
  select(snpID, rsID, full_pval, cc_pval,
         CHR, POS, contains("gene"), contains("biotype")) %>% 
  pivot_longer(intragenic_gene_1:cis_biotype_30) %>% 
  separate(name, into=c("cis","name","ID"), sep="_") %>% 
  pivot_wider() %>% 
  filter(biotype == "protein_coding") %>% 
  #collapse annotations
  group_by(snpID, rsID, full_pval, cc_pval, CHR,POS,cis) %>% 
  arrange(rsID, gene) %>% 
  summarise(genes = paste(gene, collapse=", ")) %>% 
  ungroup()

# Plots
plot.ls1 <- list()
plot.ls2 <- list()
span <- 10E3

to_plot_int <- signif_snp %>% 
  filter(cis=="intragenic") %>% pull(snpID) %>% sort()
to_plot_cis <- signif_snp %>% 
  filter(cis=="cis") %>% pull(snpID) %>% sort()
#remove non lead ZNF snps
to_plot_cis <- to_plot_cis[!to_plot_cis %in% c("3:75816211:A:G","3:75816219:T:G",
                                               "3:75816310:C:A","3:75816342:A:C")]

for(snp in to_plot_int){
  print(snp)
  #### Read in LD data ####
  load(paste0("../ReporTB_lpWGS/result/LD/",gsub(":","_",snp),".RData"))
  
  #### format data ####
  signif.temp <- signif_snp %>% filter(snpID==snp)
  
  LD.temp <- LD.df %>% 
    mutate(across(R2:Dprime, as.numeric)) %>% 
    filter(lead.snp == snp)
  
  dat.temp <- full_all %>% 
    filter(snpID %in% unique(c(LD.temp$lead.snp, LD.temp$snpID))) %>%
    select(snpID, pval) %>% 
    separate(snpID, into=c("CHR","POS","REF","ALT"), sep=":", remove=FALSE) %>% 
    mutate(POS = as.numeric(POS),
           CHR = as.numeric(CHR)) %>% 
    left_join(LD.temp, by="snpID")
  
  #### plot ####
  gene_lab <- unique(signif.temp$genes)
  #SARADH is intragenic
  if(snp=="9:133722315:T:A") { gene_lab <- gene_lab[gene_lab != "VAV2"] }
  x_lab <- wrapper(paste(sort(gene_lab),collapse=", "),24)
  x_lab <- paste0("Chromosome ", unique(signif.temp$CHR),
                  "\n(", x_lab, ")")
  
  LD1 <- dat.temp %>%
    arrange(desc(R2)) %>% 
    
    ggplot(aes(x=POS, y=-log10(pval))) +
    geom_point(data = dat.temp %>% filter(!snpID%in%signif.temp$snpID),
               aes(color=R2), shape=16) +
    #lead snp
    geom_point(data = dat.temp %>% filter(snpID%in%signif.temp$snpID),
               shape=18, size=2, color="red") +
    geom_text_repel(data = dat.temp %>% filter(snpID%in%signif.temp$snpID),
                    aes(label = snpID), 
                    segment.color = NA, size=3,
                    show.legend = FALSE, max.overlaps = 100,
                    nudge_y = 0.7, direction = "y", hjust = "center") +
    #signif line
    geom_hline(yintercept = -log10(5E-8), lty="dashed") +
    #beautify
    theme_classic(base_size = 9) +
    labs(x = x_lab,
         color="LD to lead SNP", y="-log<sub>10</sub>( p )") +
    scale_color_viridis(option = "viridis", limits=c(0,1)) +
    # scale_x_continuous(breaks=c(minC,maxC), 
    #                    labels = function(x) format(x, scientific = TRUE)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_markdown(),
          plot.title = element_text(size=10),
          # legend.direction = "horizontal", 
          # legend.position = "bottom"
          ) +
    # guides(color = guide_colorbar(title.position = "bottom",  
    #                             title.hjust = 0.5)) +
    lims(y=c(0,11.5),x=c(unique(signif.temp$POS)-10E3,unique(signif.temp$POS)+10E3))
  # LD1
  
  plot.ls1[[snp]] <- LD1
}

for(snp in to_plot_cis){
  print(snp)
  #### Read in LD data ####
  load(paste0("../ReporTB_lpWGS/result/LD/",gsub(":","_",snp),".RData"))
  
  #### format data ####
  signif.temp <- signif_snp %>% filter(snpID==snp)
  
  LD.temp <- LD.df %>% 
    mutate(across(R2:Dprime, as.numeric)) %>% 
    filter(lead.snp == snp)
  
  dat.temp <- full_all %>% 
    filter(snpID %in% unique(c(LD.temp$lead.snp, LD.temp$snpID))) %>%
    select(snpID, pval) %>%
    separate(snpID, into=c("CHR","POS","REF","ALT"), sep=":", remove=FALSE) %>% 
    mutate(POS = as.numeric(POS),
           CHR = as.numeric(CHR)) %>% 
    left_join(LD.temp, by="snpID")
  
  #### plot ####
  gene_lab <- unique(signif.temp$genes)
  #SARDH is intragenic
  if(snp=="9:133722315:T:A") { gene_lab <- gene_lab[gene_lab != "SARDH"] }
  x_lab <- wrapper(paste(sort(gene_lab),collapse=", "),24)
  x_lab <- paste0("Chromosome ", unique(signif.temp$CHR),
                  "\n(", x_lab, ")")
  
  LD1 <- dat.temp %>%
    arrange(desc(R2)) %>% 
    
    ggplot(aes(x=POS, y=-log10(pval))) +
    geom_point(data = dat.temp %>% filter(!snpID%in%signif.temp$snpID),
               aes(color=R2), shape=16) +
    #lead snp
    geom_point(data = dat.temp %>% filter(snpID%in%signif.temp$snpID),
               shape=18, size=2, color="red") +
    geom_text_repel(data = dat.temp %>% filter(snpID%in%signif.temp$snpID),
                    aes(label = snpID), 
                    segment.color = NA, size=3,
                    show.legend = FALSE, max.overlaps = 100,
                    nudge_y = 0.7, direction = "y", hjust = "center") +
    #signif line
    geom_hline(yintercept = -log10(5E-8), lty="dashed") +
    #beautify
    theme_classic(base_size = 9) +
    labs(x = x_lab,
         color="LD to lead SNP", y="-log<sub>10</sub>( p )") +
    scale_color_viridis(option = "viridis", limits=c(0,1)) +
    # scale_x_continuous(breaks=c(minC,maxC),
    #                    labels = function(x) format(x, scientific = TRUE)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_markdown(),
          plot.title = element_text(size=10),
          legend.position = "none") +
    lims(y=c(0,11.5), x=c(unique(signif.temp$POS)-10E3,unique(signif.temp$POS)+10E3))
  # LD1
  
  plot.ls2[[snp]] <- LD1
}

#### Save ####
p1 <- plot.ls1[[5]] + plot.ls1[[1]] + plot.ls1[[2]] + plot.ls1[[3]] + plot.ls1[[4]] + 
  plot.ls2[[3]] + plot.ls2[[1]] + plot.ls2[[2]] + 
  plot_layout(guides = "collect", nrow=2) +
  plot_annotation(tag_levels = list(c("A","","","","","B","","")))
# p1

ggsave(plot = p1, 
       filename="Fig2_manhattan_zoom.png",
       width=8.5, height=3.5)

