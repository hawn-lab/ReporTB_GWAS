library(tidyverse)
# library(biomaRt)
# select <- dplyr::select
score_P <- 1E-4
signif <- 5E-8
sugg <- 5E-5

#MAF 5% and concordant
maf5_concord <- read_tsv("result/snp_imp_maf5_concord.txt", col_names=FALSE) %>% 
  pull(X1)

#Full models
## Full cohort
full <- read_csv("result/model_final/model_full.csv") 
## Case-control
cc <- read_csv("result/model_final/model_cc.csv")

#### Full models ####
#combine
dat <- bind_rows(full, cc) %>% 
  filter(term == "value" & pval < sugg & snpID %in% maf5_concord) %>% 
  select(model, snpID, pval, OR) %>% 
  #Get positions and alleles
  separate(snpID, into=c("CHR","POS","REF","ALT"), sep=":", remove = FALSE) %>% 
  mutate(POS = as.numeric(POS),
         CHR = as.numeric(CHR))

#### rsID annotation ####
## get CHR ID
chr.key <- read_tsv("../metadata/seq_region.txt", col_names = FALSE) %>% 
  rename(CHR_ID=X1, CHR=X2) %>% 
  select(-X3) %>% 
  filter(CHR %in% dat$CHR) %>% 
  arrange(CHR_ID)
summary(chr.key$CHR_ID)

#### Bash ####
## Keep columns of interest
# awk -F '\t' '{ print $2,$3,$7,$9,$13 }' ../metadata/snp_key/variation_feature.txt > ../metadata/snp_key/variation_feature_select.txt 

## Split by CHR
# for i in {131537..131560};
# do
#   echo $i
#   awk -F '\t' -v i="$i" 'index($1, i) == 1' ../metadata/snp_key/variation_feature_select.txt > ../metsdata/snp_key/$i.txt
# done

#### R ####
# library(R.utils)
# gzip('../metadata/snp_key/variation_feature_select.txt',
#      destname='../metadata/snp_key/variation_feature_select.txt.gz')
# 
# # Read in and filter each CHR by positions of SNP in data
# files <- list.files(path = "../metadata/snp_key/", pattern = ".txt",
#                     full.names = TRUE)
# 
# key.ls <- list()
# for(f in files){
#   print(f)
#   #Load key for CHR of interest
#   temp <- data.table::fread(f,
#                             col.names = c("CHR_ID","POS",
#                                           "allele","rsID","anno"),
#                             fill=TRUE)
# 
#   #Get actual CHR name
#   ID <- unique(temp$CHR_ID)
#   ID.CHR <- chr.key %>% filter(CHR_ID == ID) %>% pull(CHR)
# 
#   if(length(ID.CHR) > 0){
#     #Filter data to CHR of interest
#     dat.temp <- dat %>%
#       filter(CHR == ID.CHR) %>%
#       select(CHR,POS) %>%
#       distinct()
# 
#     #Filter key to positions in data
#     key.ls[[as.character(ID.CHR)]] <- temp %>%
#       mutate(CHR = as.numeric(ID.CHR)) %>%
#       inner_join(dat.temp, by=c("CHR","POS"))
#   }}
# 
# save(key.ls, file="../metadata/snp_key/snp_key.RData")
# 
# ## convert to df
# key.df <- as.data.frame(do.call(bind_rows, key.ls))
# 
# ## Match to our data
# snp_db_clean1 <- key.df %>%
#   #Keep only SNP
#   filter(grepl("^[A-Z]{1}\\/[A-Z]{1}$", allele)) %>%
#   separate(allele, into = c("REF","ALT"), sep="/", remove = FALSE) %>%
#   distinct() %>%
#   mutate(snp_ref = "orig")
# 
# #Reverse ref and alt alleles
# snp_db_clean2 <- snp_db_clean1 %>%
#   mutate(REF2 = REF, ALT2 = ALT,
#          ALT = REF, REF = ALT2,
#          allele = paste(REF, ALT, sep="/")) %>%
#   select(-REF2, -ALT2) %>%
#   distinct() %>%
#   mutate(snp_ref = "reverse")
# 
# snp_db_clean <- bind_rows(snp_db_clean1, snp_db_clean2) %>%
#   distinct() %>%
#   group_by(rsID, allele, REF, ALT, CHR, POS, anno) %>%
#   summarise(snp_ref = paste(unique(snp_ref), collapse=",")) %>%
#   ungroup() %>%
#   mutate(snpID = paste(CHR,POS,REF,ALT,sep=":"))
# save(key.ls, snp_db_clean, file="../metadata/snp_key/snp_key.RData")

load("../metadata/snp_key/snp_key.RData")

#### Gene annotation ####
#Suggestive SNP
sugg.snp <- dat %>% 
  distinct(snpID, CHR, POS, REF, ALT) %>% 
  #add 50kb cis
  mutate(POS1 = POS-50E3, 
         POS2 = POS+50E3,
         POS1 = ifelse(POS1 < 0, 0, POS1),
         POS2 = ifelse(POS2 < 0, 0, POS2)) %>% 
  mutate(chrID =paste(CHR,POS1,POS2, sep = ":"))

#### Get reference genome ####
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene.key <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                                          'gene_biotype',
                                          'chromosome_name', 'start_position', 
                                          'end_position'), 
                           filters = 'chromosomal_region', 
                           values = unique(sugg.snp$chrID), 
                           mart = ensembl)

gene.key.format <- full_join(sugg.snp, gene.key, by=c("CHR"="chromosome_name"),
                             relationship = "many-to-many") %>% 
  filter(POS1 <= end_position & POS2 >= start_position) %>% 
  select(-c(POS1,POS2,chrID)) %>% 
  #mark cis vs internal to a gene
  mutate(cis = case_when(POS <= end_position & 
                           POS >= start_position ~ "internal",
                         TRUE ~ "cis")) %>% 
  select(snpID, ensembl_gene_id, hgnc_symbol,gene_biotype, cis,
         POS, start_position, end_position)

#proten coding only
pc_biotypes <- c("protein_coding", "IG_C_gene", "IG_V_gene", 
                 "TR_C_gene", "TR_V_gene")
gene.key.pc <- gene.key %>% 
  filter(gene_biotype %in% pc_biotypes)

gene.key.format.pc <- full_join(sugg.snp, gene.key.pc, 
                                by=c("CHR"="chromosome_name"),
                                relationship = "many-to-many") %>% 
  filter(POS1 <= end_position & POS2 >= start_position) %>% 
  select(-c(POS1,POS2,chrID)) %>% 
  #mark cis vs internal to a gene
  mutate(cis = case_when(POS <= end_position & POS >= start_position ~ "internal",
                         TRUE ~ "cis")) %>% 
  select(snpID, ensembl_gene_id, hgnc_symbol,gene_biotype, cis,
         POS, start_position, end_position)

#### Map snp to genes ####
#map to internal within pc gene
gene.key.int.pc <- gene.key.format.pc %>% 
  filter(cis == "internal") %>% 
  #replace gene with ensembl if no symbol available
  mutate(intragenic_gene = ifelse(!is.na(hgnc_symbol) & 
                                    !hgnc_symbol %in% c(""," "),
                                  hgnc_symbol, ensembl_gene_id)) %>% 
  rename(intragenic_biotype=gene_biotype,
         intragenic_ensembl=ensembl_gene_id) %>% 
  select(snpID, intragenic_gene, intragenic_ensembl, intragenic_biotype)

#map to internal within non-pc gene
gene.key.int.nonpc <- gene.key.format %>% 
  filter(cis == "internal" & ! gene_biotype %in% pc_biotypes) %>% 
  #replace gene with ensembl if no symbol available
  mutate(intragenic_gene = ifelse(!is.na(hgnc_symbol) & 
                                    !hgnc_symbol %in% c(""," "),
                                  hgnc_symbol, ensembl_gene_id)) %>% 
  rename(intragenic_biotype=gene_biotype,
         intragenic_ensembl=ensembl_gene_id) %>% 
  select(snpID, intragenic_gene, intragenic_ensembl, intragenic_biotype)

#map cis to pc genes
gene.key.cis.pc <- gene.key.format.pc %>% 
  #remove intragenic pc
  # filter(!snpID %in% gene.key.int.pc$snpID) %>% 
  filter(cis == "cis") %>% 
  #Calculate distance to gene
  rowwise() %>% 
  mutate(distance = min(c(abs(start_position-POS), abs(POS-end_position)))) %>%
  #determine closest
  group_by(snpID) %>% 
  mutate(cis = case_when(cis != "internal" & distance == min(distance) ~ "cis_closest",
                         TRUE ~ cis)) %>% 
  ungroup() %>% 
  #replace gene with ensembl if no symbol available
  mutate(cis_gene = ifelse(!is.na(hgnc_symbol) & 
                             !hgnc_symbol %in% c(""," "),
                           hgnc_symbol, ensembl_gene_id)) %>% 
  #label
  rename(cis_biotype=gene_biotype, distance_to_cis=distance,
         cis_ensembl=ensembl_gene_id) %>% 
  select(snpID, cis_gene, cis_ensembl, cis_biotype, distance_to_cis)

#map cis to non-pc gene
gene.key.cis.nonpc <- gene.key.format %>% 
  #remove intragenic pc
  # filter(!snpID %in% gene.key.int.pc$snpID) %>% 
  filter(cis == "cis") %>% 
  #Calculate distance to gene
  rowwise() %>% 
  mutate(distance = min(c(abs(start_position-POS), abs(POS-end_position)))) %>%
  #determine closest
  group_by(snpID) %>% 
  mutate(cis = case_when(cis != "internal" & distance == min(distance) ~ "cis_closest",
                         TRUE ~ cis)) %>% 
  ungroup() %>% 
  #replace gene with ensembl if no symbol available
  mutate(cis_gene = ifelse(!is.na(hgnc_symbol) & 
                             !hgnc_symbol %in% c(""," "),
                           hgnc_symbol, ensembl_gene_id)) %>% 
  #label
  rename(cis_biotype=gene_biotype, distance_to_cis=distance,
         cis_ensembl=ensembl_gene_id) %>% 
  select(snpID, cis_gene, cis_ensembl, cis_biotype, distance_to_cis)

#COMBINE
gene.key.cis.wide <- full_join(gene.key.cis.pc,gene.key.cis.nonpc) %>% 
  arrange(snpID, distance_to_cis) %>% 
  mutate(distance_to_cis = as.character(distance_to_cis)) %>% 
  rename(cis_distance = distance_to_cis)  %>% 
  pivot_longer(-c(snpID)) %>% 
  group_by(snpID, name) %>% 
  mutate(id = row_number()) %>% 
  mutate(name = paste(name, id, sep="_")) %>% 
  select(-id) %>% 
  ungroup() %>% 
  pivot_wider() %>% 
  mutate(across(contains("distance"), as.numeric))

gene.key.int.wide <- full_join(gene.key.int.pc, gene.key.int.nonpc) %>% 
  pivot_longer(-c(snpID)) %>% 
  group_by(snpID, name) %>% 
  mutate(id = row_number()) %>% 
  mutate(name = paste(name, id, sep="_")) %>% 
  select(-id) %>% 
  ungroup() %>% 
  pivot_wider()


gene.key.all <- full_join(gene.key.int.wide, gene.key.cis.wide) %>% 
  distinct() #%>% 
  # #Create main anno group
  # mutate(anno_group = case_when(
  #   #in pc gene
  #   intragenic_biotype %in% pc_biotypes ~ "intra_protein",
  #   #cis to pc gene
  #   cis_biotype %in% pc_biotypes ~ "cis_protein",
  #   #in nonpc gene
  #   !intragenic_biotype %in% pc_biotypes & !is.na(intragenic_biotype) ~
  #     "intra_nonprotein",
  #   #cis nonpc gene
  #   ! cis_biotype %in% pc_biotypes & !is.na(cis_biotype) ~
  #     "cis_nonprotein",
  #   #not in a gene
  #   is.na(intragenic_gene) & is.na(cis_gene)~"no_anno"))  %>% 
  # select(snpID, anno_group, everything())

#### Format SNP results ####
# Add key anno group
gene.snp.key.all <- snp_db_clean %>% 
  select(snpID, rsID, anno, snp_ref) %>% 
  # full_join(gene.key.all %>% select(snpID, anno_group)) %>%
  filter(snpID %in% sugg.snp$snpID) %>% 
  #Get positions
  separate(snpID, into=c("CHR","POS"), sep=":", remove = FALSE) %>% 
  select(snpID,rsID,CHR,POS) 

save(key.ls, snp_db_clean, gene.snp.key.all,
     file="../metadata/snp_key/snp_key.RData")

#### Summary tables ####
fdr <- bind_rows(full, cc) %>% 
  select(model, term, snpID, pval, OR) %>% 
  filter(snpID %in% dat$snpID & term=="value") %>% 
  select(-term) %>% 
  pivot_longer(-c(model, snpID)) %>% 
  mutate(name = paste0(model, "_", name)) %>% 
  select(-model) %>% 
  pivot_wider() %>% 
  mutate(full_group = case_when(full_pval < signif ~"signif",
                                full_pval < sugg ~ "suggestive",
                                TRUE ~ "NS"),
         cc_group = case_when(cc_pval < signif ~"signif",
                              cc_pval < sugg ~ "suggestive",
                              TRUE ~ "NS"))  %>% 
  full_join(gene.snp.key.all %>% select(snpID, rsID)) %>% 
  filter(full_pval < sugg | cc_pval < sugg) %>% 
  mutate(full_group = ifelse(is.na(full_group),"not_tested",full_group)) %>% 
  #Get position columns and alleles
  separate(snpID, into=c("CHR","POS","REF","ALT"), sep=":", remove = FALSE) %>% 
  mutate(POS = as.numeric(POS),
         CHR = as.numeric(CHR)) %>% 
  select(snpID, rsID, full_pval, cc_pval, full_group, cc_group,
         CHR, POS, REF, ALT, everything()) %>% 
  arrange(full_pval)

#### MAF within group ####
library(SNPRelate)
GDS_file <- "result/gds/reportTB_filter_imp_maf5.gds"
genofile <- snpgdsOpen(GDS_file)
geno.dat <- snpgdsGetGeno(genofile, 
                          snp.id = unique(fdr$snpID),
                          with.id=TRUE)
geno.num <- geno.dat[["genotype"]]
snpgdsClose(genofile)

## Add row and columnnames
colnames(geno.num) <- geno.dat[["snp.id"]]
rownames(geno.num) <- sapply(strsplit(geno.dat[["sample.id"]],"_"), `[`, 1)

attach("../metadata/ReportTB_metadata_clean.RData")

geno.num.select <- as.data.frame(geno.num) %>% 
  rownames_to_column("sample_id") %>% 
  inner_join(meta_FULL %>% 
               select(sample_id,tb) %>% 
               drop_na()) %>% 
  select(sample_id,tb, everything())

maf_summ <- geno.num.select %>% 
  #total participants per group
  group_by(tb) %>% 
  mutate(n = n()) %>%
  pivot_longer(-c(sample_id,tb,n), names_to = "snpID") %>% 
  #minor allele per group
  group_by(tb,snpID,n) %>% 
  summarise(minor = sum(value)) %>% 
  ungroup() %>% 
  #entire cohort
  group_by(snpID) %>% 
  mutate(minor2 = sum(minor)) %>% 
  ungroup() %>% 
  #calc maf
  rowwise() %>% 
  mutate(maf_group=minor/(2*n),
         MAF_overall=minor2/(2*nrow(geno.num))) %>% 
  select(-n,-minor,-minor2) %>% 
  pivot_wider(names_from = tb, values_from = maf_group) %>% 
  rename(MAF_control=contact,
         MAF_case=culture_confirm)

# count alleles per group
ale_summ <- geno.num.select %>%
  pivot_longer(-c(sample_id,tb), names_to = "snpID") %>% 
  group_by(tb, snpID, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(tb = recode(tb, "culture_confirm"="case",
                     "contact"="control"),
         name = paste(tb, value, sep="_")) %>% 
  select(-tb,-value) %>% 
  pivot_wider(values_from = n)
ale_summ[is.na(ale_summ)] <- 0

fdr_maf <- fdr %>% 
  left_join(maf_summ) %>% 
  left_join(ale_summ)

#### Add cov notes ####
fdr_maf_notes <- fdr_maf %>% 
  left_join(read_csv("result/model_final/ReporTB_signif_snp_notes.csv")) %>% 
  left_join(read_csv("result/model_final/ReporTB_signif_cov_notes.csv"))

#### Add gene annotation ####
fdr_maf_notes_anno <- fdr_maf_notes %>% 
  left_join(gene.key.all)

#### Save ####
dir.create("result/model_final_anno/", showWarnings = FALSE)
write_csv(fdr_maf_notes_anno,
          file = "result/model_final_anno/ReporTB_sugg_snp_anno.csv", 
          na = "")

fdr_maf_notes_anno %>% 
  filter(full_group=="signif"|cc_group=="signif") %>% 
  write_csv(file = "result/model_final_anno/ReporTB_signif_snp_anno.csv", 
          na = "")

#List all genes, pc
temp <- gene.key.all %>% 
  filter(snpID %in% fdr_maf_notes_anno$snpID) %>% 
  select(snpID, contains("gene"), contains("biotype")) %>% 
  pivot_longer(-snpID) %>% 
  drop_na(value) %>% 
  separate(name, into=c("anno","name","id"), sep="_") %>% 
  pivot_wider()

temp %>% 
  filter(biotype %in% pc_biotypes) %>% 
  distinct(gene) %>% 
  arrange(gene) %>% 
  write_csv("result/model_final_anno/ReporTB_signif_snp_pc_genes.csv", 
            col_names = FALSE)

#List all genes, any biotype
temp %>% 
  distinct(gene) %>% 
  arrange(gene) %>% 
  write_csv("result/model_final_anno/ReporTB_signif_snp_all_genes.csv")

