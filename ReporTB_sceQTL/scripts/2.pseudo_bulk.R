#Single cell eQTL analysis
#Create pseudo bulk normalized counts
library(tidyverse)
library(readxl)
library(Seurat)
library(edgeR)
library(limma)
library(RNAetc)
# library(biomaRt)
library(BIGpicture)
library(patchwork)

#### Meta data ####
#Study metadata
load("../metadata/ReportTB_metadata_clean.RData")
#Add progressors
prog <- read_excel("data_raw/metadata/progressor metadata Kim May6 2024 MCF.xlsx") %>% 
  rename(subjid=pid, tb=TB, sex=Sex, bl_age=Age, bl_hiv=HIV,
         diabetes="Diabetes - self reported") %>% 
  mutate(sex = recode(sex, "Male"="M", "Female"="F"),
         bl_age = round(bl_age),
         bl_hiv = recode(bl_hiv, "Negative"="N"),
         diabetes = recode(diabetes, "No"="N", "Yes"="Y"),
         smokhx = fct_recode(factor(smokhx), "never"="Never",
                             "past"="Used to smoke",
                             "current"="Currently smoke"))

meta_prog <- bind_rows(meta_FULL, prog) %>% 
  mutate(smokhx2 = case_when(smokhx %in% c("current","past")~"Y",
                             smokhx == "never"~"N",
                             TRUE~smokhx),
         smokhx2 = factor(smokhx2, levels=c("N","Y")))

#Fill in more prog data
meta_addtl <- read_excel("data_raw/metadata/Add data Kim May21 RePORTBrasilFase1Coo-OffStudyPIDsDevelope_DATA_LABELS_2024-05-21_0812.xlsx") %>% 
  janitor::clean_names() %>% 
  group_by(participante_id_coorte_b) %>% 
  fill(raca:tb_activation, .direction="downup") %>% 
  mutate(across(baseline_visit:tb_activation, ~as.Date(.))) %>% 
  filter(event_name=="Visita Inicial") %>% 
  #recode to English
  mutate(race = recode(raca, "Branco"="white", "Negro"="black", "Pardo"="brown")) %>% 
  #calc other variables
  mutate(bmi = weight_b0_kg/((height_b0_cm/100)^2),
         days_followup = as.numeric(left_study-baseline_visit)) %>% 
  rename(bid=participante_id_coorte_b) %>% 
  select(bid, race, bmi, days_followup)

meta_prog <- meta_prog %>% 
  left_join(meta_addtl, by=c("subjid"="bid")) %>% 
  mutate(race = coalesce(race.x, race.y), 
         bmi = coalesce(bmi.x, bmi.y),
         days_followup = coalesce(days_followup.x, days_followup.y), 
         .before = resregion) %>% 
  select(!contains(".x") & !contains(".y"))

#### Single cell RNAseq data ####
sc <- readRDS("data_raw/scRNAseq/LTBI_Seurat_Object_12DEC2024.rds")

# Check counts
mem.maxVSize(1000000)
# library(Matrix)
# temp <- sc@assays[["RNA"]]@layers[["counts"]]
# temp <- as.matrix(temp)
# temp[1:3,1:3]
# rm(temp)

#Add UMAP coordinates
umap <- read_csv("data_raw/scRNAseq/LTBI_supplemental_metadata.csv") %>% 
  select(UMAP1, UMAP2)
sc@meta.data <- as.data.frame(sc@meta.data) %>% 
  bind_cols(umap)

sc_meta <- as.data.frame(sc@meta.data)

save(sc, file="data_clean/LTBI_scRNAseq_seurat_anno.RData")

# Save cell annotation key
cell_key <- sc_meta %>% 
  distinct(annotation) %>% 
  mutate(anno_long = gsub("_|-"," ",annotation)) %>% 
  #Recode reactive monocytes
  mutate(anno_long = recode(anno_long,
                            "Mono GRV reactive"="CD14 Mono 2",
                            "Mono MTB300 reactive"="CD14 Mono 1",
                            "CD14 Mono"="CD14 Mono 3")) %>% 
  #make long names
  mutate(anno_long = recode(
    anno_long,
    "act CD4"="CD4 act: T-cell CD4+ activated",
    "B"="B: B cell",
    "CD14 Mono 3"="CD14 Mono-3: Monocyte CD14+, cluster 3",
    "CD16 Mono"="CD16 Mono: Monocyte CD16+",
    "CD4 Naive"="CD4 naive-4: T-cell CD4+ naive, cluster 4",
    "CD4 TCM"="CD4 CM: T-cell CD4+ central memory",
    "CD4 TEM"="CD4 EM: T-cell CD4+ effector memory",
    "CD4 Treg"="CD4 Treg: T-cell CD4+ regulatory",
    "CD4M 1"="CD4 mem-1: T-cell CD4+ memory, cluster 1",
    "CD4M 2"="CD4 mem-2: T-cell CD4+ memory, cluster 2",
    "CD4M 3"="CD4 mem-3: T-cell CD4+ memory, cluster 3",
    "CD4N 1"="CD4 naive-1: T-cell CD4+ naive, cluster 1",
    "CD4N 2"="CD4 naive-2: T-cell CD4+ naive, cluster 2",
    "CD4N 3"="CD4 naive-3: T-cell CD4+ naive, cluster 3",
    "CD8 Naive"="CD8 naive: T-cell CD8+ naive",
    "CD8 TCM"="CD8 CM: T-cell CD8+ central memory",
    "CD8 TEM"="CD8 EM: T-cell CD8+ effector memory",
    "cDC"="cDC: Dendritic cell, conventional",
    "gdT"="gdT: T-cell gamma delta",
    "HSP CD4"="CD4 HSP: T-cell CD4+ heat shock protein+",
    "HSPC"="HSPC: Hematopoietic stem & progenitor cell",
    "Int Mono"="CD14 Mono int: Monocyte CD14+, intermediate",
    "ISG CD4"="CD4 ISG: T-cell CD4+ IFN stimulated genes+",
    "MAIT"="MAIT: T-cell mucosal-associated invariant",
    "CD14 Mono 2"="CD14 Mono-2: Monocyte CD14+, cluster 2",
    "CD14 Mono 1"="CD14 Mono-1: Monocyte CD14+, cluster 1",
    "NK"="NK-3: Natural killer cell, cluster 3",   
    "NK 1"="NK-1: Natural killer cell, cluster 1",
    "NK 2"="NK-2: Natural killer cell, cluster 2",
    "pDC"="pDC: Dendritic cell, plasmacytoid",
    "Plasmablast"="Plasmablast",
    "Platelet"="Platelet"
   ),
   #Make short names
    anno_short = gsub(":.*","",anno_long))

write_csv(cell_key, "result/sc_key.csv")

#### Explore sc metadata ####
##Total donors per group
sc_meta %>% 
  distinct(bid, progr_contr) %>% 
  count(progr_contr)

sc_meta %>% 
  distinct(bid, progr_contr) %>% 
  left_join(meta_prog %>% select(subjid,tb), by=c("bid"="subjid")) %>% 
  count(tb)

#Unique samples
sc_meta %>% 
  rownames_to_column() %>% 
  distinct(bid,  condition) %>%
  count(bid, condition)

## Cell types
sort(unique(sc_meta$annotation))
length(unique(sc_meta$annotation))

#### Pseudo bulk ####
# pseudo_sc <- AggregateExpression(sc, assays = "RNA",
#                                  return.seurat = TRUE,
#                                  group.by = c("bid",
#                                               "condition",
#                                               "annotation"))
# save(pseudo_sc, sc_meta,
#      file="data_clean/LTBI_scRNAseq_seurat_pseudobulk.RData")
# beepr::beep()

#Keep environment size down
rm(sc) 
gc()
load("data_clean/LTBI_scRNAseq_seurat_pseudobulk.RData")

#Extract counts
count <- as.data.frame(pseudo_sc[["RNA"]]$counts)

#### Filter protein coding ####
#Protein coding
ensembl <- biomaRt::useEnsembl(biomart="ensembl", 
                               dataset="hsapiens_gene_ensembl")

#Format gene key
key <- biomaRt::getBM(attributes=c("ensembl_gene_id",
                                   "hgnc_symbol", 
                                   "gene_biotype"), 
                      mart=ensembl) %>% 
  #Filter protein coding genes
  filter(gene_biotype == "protein_coding")

key.filter <- key %>% 
  #Filter protein coding genes in count table
  filter(hgnc_symbol %in% rownames(count)) %>% 
  #collapse multiannotations.
  group_by(hgnc_symbol, gene_biotype) %>% 
  summarise(ensembl_gene_id = list(unique(ensembl_gene_id)), .groups = "drop")

count_pc <- count %>% 
  rownames_to_column() %>% 
  filter(rowname %in% key.filter$hgnc_symbol) %>% 
  column_to_rownames() %>% 
  t() %>% as.data.frame()

#### Create metadata ####
##Sample ID key
key_pt <- read_csv("data_raw/metadata/ReportTB_sample_sheet_clean.csv") %>% 
  distinct(sample_id, sample_description) %>% 
  rename(bid=sample_description) %>% 
  #pt in data
  filter(bid %in% sc_meta$bid) 

bulk_meta <- data.frame(libID = rownames(count_pc)) %>% 
  #Info in library names
  separate(libID, into = c("bid","condition","annotation"), sep="_",
           remove = FALSE) %>% 
  #add single cell metadata
  full_join(sc_meta %>% distinct(bid, condition, progr_contr),
            by = c("bid","condition")) %>% 
  mutate(TB = factor(progr_contr, levels = c("controller","progressor")),
         condition = factor(condition, levels=c("Media","GRV","MTB300"))) %>% 
  #Add snp metadata
  left_join(key_pt) %>% 
  left_join(meta_prog)

#### Add PCs to metadata ####
attach("result/kinship/reportTB_sc_PCair.RData")

pc <- as.data.frame(PC.dat) %>%
  rownames_to_column("sample_id") %>% 
  separate(sample_id, into=c("sample_id"), sep="_", extra="drop") %>% 
  select(sample_id:PC10)

bulk_meta_pc <- bulk_meta %>% 
  left_join(pc, by="sample_id")

#### DGE object ####
##Remove libraries missing in one dataset
bulk_meta_ord <- bulk_meta_pc %>% 
  arrange(libID) %>% 
  filter(libID %in% rownames(count_pc)) %>% 
  filter(condition=="Media")

write_csv(file = "data_clean/ReporTB_metadata_sc_clean.csv",
          bulk_meta_ord)

count_pc_ord <- count_pc %>% 
  rownames_to_column() %>% 
  arrange(rowname) %>% 
  filter(rowname %in% bulk_meta_ord$libID) %>% 
  column_to_rownames() %>% 
  t()

key.filter.ord <- key.filter %>% 
  arrange(hgnc_symbol)

dat <- DGEList(
  #count table
  counts=count_pc_ord,
  #metadata
  samples=bulk_meta_ord,
  #gene info
  genes=key.filter.ord)

#Save
save(dat, file="data_clean/ReporTB_sc_dat_unnorm.RData")

#### PER CELL TYPE ####
types <- sort(unique(bulk_meta_ord$annotation))
removed_genes <- c()
retained_genes <- c()

for(cell in types){
  meta.temp <- bulk_meta_ord %>%
    filter(annotation == cell) %>% 
    left_join(cell_key, by = join_by(annotation))
  
  counts.temp <- as.data.frame(count_pc_ord) %>%
    select(all_of(meta.temp$libID))
  
  cell_name <- unique(meta.temp$anno_short)
  cell_name <- gsub(" ","_",cell_name)
  print(cell_name)
  
  #DGE object
  dat.temp <- DGEList(
    #count table
    counts=counts.temp,
    #metadata
    samples=meta.temp,
    #gene info
    genes=key.filter.ord)
  
  assign(paste0("dat_",cell_name), dat.temp, envir = .GlobalEnv)
  
  #Filter 0 genes
  low_count <- min(dat.temp$counts[dat.temp$counts>0])-0.0001
  dat.abund <- RNAetc::filter_rare(dat.temp, min_CPM = low_count, min.sample = 3,
                                   gene.var="hgnc_symbol")
  removed_genes <- c(nrow(dat.temp$counts)-nrow(dat.abund$counts), removed_genes)
  names(removed_genes)[1] <- cell
  retained_genes <- c(nrow(dat.abund$counts), retained_genes)
  names(retained_genes)[1] <- cell
  
  #### Normalize ####
  #TMM norm
  print("TMM")
  dat.abund.norm <- calcNormFactors(dat.abund, method = "TMM")
  
  #Log2 CPM + weights
  print("voom")
  dat.abund.norm.voom <- voomWithQualityWeights(
    dat.abund.norm,
    design=model.matrix(~ tb,
                        data=dat.abund.norm$samples),
    plot=TRUE)
  
  assign(paste0("voom_", cell_name),
         dat.abund.norm.voom, envir = .GlobalEnv)
}

#Total genes
min(retained_genes)
max(retained_genes)

#Save
obs1 <- ls(pattern="dat_")
obs2 <- ls(pattern="voom_")

save(list=obs1, file="data_clean/ReporTB_sc_dat.RData")
save(list=obs2, file="data_clean/ReporTB_sc_voom.RData")
beepr::beep()

#### Sample list for SNP data ####
#Actual SNP sample names
fam <- read_tsv("../ReporTB_lpWGS/result/plink/reportTB.fam", col_names = FALSE) %>% 
  select(X2) %>% 
  separate(X2, into=c("sample_id"), sep="_", remove=FALSE)

#Save file for PLINK filtering
bulk_meta %>% 
  distinct(sample_id) %>% 
  left_join(fam) %>% 
  
  mutate(X1=0) %>% 
  select(X1,X2) %>% 
  write_tsv(col_names = FALSE, 
            file = "result/sc_patients.txt")
