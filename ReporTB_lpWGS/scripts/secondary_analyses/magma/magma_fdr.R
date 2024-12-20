library(tidyverse)
library(data.table)

#### Gene analysis ####
gene <- fread("result/magma/ReporTB_magma_gene.genes.out",sep=" ")

#Annotate to HGNC symbol
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene.key <- biomaRt::getBM(attributes = c('entrezgene_id', 'ensembl_gene_id', 
                                 'hgnc_symbol', 'gene_biotype'), 
                  filters = 'entrezgene_id', 
                  values = unique(gene$GENE), 
                  mart = ensembl) %>% 
  group_by(entrezgene_id) %>% 
  summarise(ensembl_gene_id = list(unique(ensembl_gene_id)),
            hgnc_symbol = list(unique(hgnc_symbol)),
            gene_biotype = list(unique(gene_biotype))) %>% 
  ungroup() %>% 
  mutate(across(ensembl_gene_id:gene_biotype, ~as.character(.)))

gene.format <- gene %>% 
  left_join(gene.key, by=c("GENE"="entrezgene_id")) %>% 
  rename(entrezgene_id=GENE)  %>% 
  #FDR correct
  mutate(FDR_B = p.adjust(P_MULTI, method="bonferroni"),
         FDR_BH = p.adjust(P_MULTI, method="fdr")) %>% 
  select(entrezgene_id, ensembl_gene_id, hgnc_symbol, CHR:P_MULTI, 
         FDR_B, FDR_BH,
         everything()) %>% 
  arrange(P_MULTI)

write_csv(gene.format, file="result/magma/ReporTB_magma_gene_fdr.csv", na = "")  

#### Gene set analysis ####
geneset <- fread("result/magma/ReporTB_magma_geneset.gsa.out",sep=" ", skip=4) %>% 
  separate(FULL_NAME, into=c("db","pathway"), sep="_", extra="merge", 
           remove=FALSE) %>% 
  mutate(pathway = gsub("_"," ",pathway)) %>% 
  filter(NGENES >= 2)

#Get genes in terms
library(msigdbr)

msig <- msigdbr("human", "H") %>% 
  bind_rows(msigdbr("human", "C2")) %>% 
  bind_rows(msigdbr("human", "C5")) %>% 
  select(gs_cat:entrez_gene) %>% 
  # genes in GWAS dataset
  filter(entrez_gene %in% gene.format$entrezgene_id) %>% 
  # gene sets in MAGMA results
  filter(gs_name %in% geneset$FULL_NAME) %>% 
  # summarise pathways
  arrange(gs_cat, gs_subcat, gene_symbol) %>% 
  group_by(gs_cat, gs_subcat, gs_name) %>% 
  summarise(GENES = list(unique(gene_symbol))) %>% 
  ungroup() 

geneset.format <- geneset %>% 
  left_join(msig, by=c("FULL_NAME"="gs_name")) %>% 
  group_by(gs_cat, gs_subcat) %>% 
  mutate(FDR_B = p.adjust(P, method="bonferroni"),
         FDR_BH = p.adjust(P, method="fdr")) %>% 
  select(FULL_NAME, gs_cat, gs_subcat, pathway, BETA:P, 
         FDR_B, FDR_BH, NGENES,GENES) %>% 
  arrange(P)
  
geneset.format %>% 
  mutate(GENES = as.character(GENES)) %>% 
  write_csv(file="result/magma/ReporTB_magma_geneset_fdr.csv")  

#### Signif #### 
gene.format %>% 
  filter(FDR_BH < 0.2)

geneset.format %>% 
  filter(FDR_BH < 0.2 & NGENES > 1) 
