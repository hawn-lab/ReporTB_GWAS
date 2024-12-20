# EXAMPLE
# test <- reportb_eqtl(dat = ls(pattern = "voom_"),
#              models = c("~condition*genotype + PC1 + PC2 + (1|sample_id)",
#                         "~condition*genotype + PC1 + PC2 + sex + (1|sample_id)"),
#              group = "univariate",
#              processors = 6)

reportb_eqtl <- function(dat = NULL, models = NULL, group = NULL, processors = 6){
  require(parallel)
  require(foreach)
  require(tidyverse)
  require(kimma)
  
  print("Running models")
  model_ls <- NULL
  cl <- parallel::makeCluster(processors)
  doParallel::registerDoParallel(cl)
  
  model_ls <- foreach::foreach(voom = dat,
                               .export = ls(globalenv()),
                               .packages = c("dplyr", "magrittr", "stats", "broom", 
                                            "lme4", "car", "tibble", "lme4qtl", "utils", "emmeans", 
                                            "foreach", "doParallel","kimma")) %dopar% {
    print(voom)
    voom.temp <- get(voom)
    
    #Media only samples
    ##meta
    voom.temp$targets <- voom.temp$targets %>%
      filter(condition=="Media")
    ##weights
    colnames(voom.temp$weights) <- colnames(voom.temp$E)
    rownames(voom.temp$weights) <- rownames(voom.temp$E)
    voom.temp$weights <- as.data.frame(voom.temp$weights) %>%
      select(all_of(voom.temp$targets$libID)) %>%
      as.matrix()
    ##expression
    voom.temp$E <- as.data.frame(voom.temp$E) %>%
      select(all_of(voom.temp$targets$libID)) %>%
      as.matrix()
    
    #Rm genes not in expression. Differs for indiv cell types
    map.rename.temp <- map.rename %>% 
      filter(hgnc_symbol %in% rownames(voom.temp$E))
    
    #Rm genotypes with only 1 level and without matching gene
    geno.num.rename.temp <- geno.num.rename %>%
      filter(sample_id %in% voom.temp$targets$sample_id) %>% 
      select(sample_id, all_of(unique(map.rename.temp$snpID)))%>%
      select_if(function(col) length(unique(col))>1) 
    
    #Run models
    for(m in models){
      #Short model name
      if(group=="univariate"){
        model_name <- strsplit(m, split=" [+] ")[[1]][4]
        if(model_name == "(1|sample_id)"){ model_name <- "base" }
        
      } else if(group=="loo"){
        all <- "~genotype + PC1 + PC2 + sex + bl_age + smokhx2 + diabetes + (1|sample_id)"
        
        if(m==all){
          model_name <- "full"
        } else{
          all_split <- strsplit(all, split=" [+] ")[[1]]
          
          model_name <- strsplit(m, split=" [+] ")[[1]]
          model_name <- all_split[!all_split %in% model_name]
          model_name <- paste0("minus.",model_name)
        }
      } else{
        stop("group parameter not recognized.")
      }
      
      model_ls[[paste(voom,model_name,sep="_")]] <- kmFit_eQTL(
        dat_snp=geno.num.rename.temp, 
        dat=voom.temp, 
        dat_map=map.rename.temp,
        kin=kin.matrix,
        geneID="hgnc_symbol", genotypeID="snpID",
        patientID="sample_id",
        model=m,
        use_weights=TRUE, run_lmerel=TRUE, run_contrast = FALSE,
        metrics=TRUE, processors=1)
    }
    return(model_ls)
  }
  
  parallel::stopCluster(cl)
  
  print("Saving models")
  # Combine results into df
  # Only keep fit results since need to rerun model for contrasts anyway
  model_uni_fit <- data.frame()
  model_names <- map_depth(model_ls, 1, names) %>% unlist(use.names = FALSE)
  
  #One level per processor
  for(i in 1:length(model_ls)){
    model.ls.temp <- model_ls[[i]]
    names.temp <- names(model.ls.temp)
    
    for(m2 in names.temp){
      model.temp <- model.ls.temp[[m2]]
      temp <- strsplit(m2, split="voom_")[[1]][2]
      temp <- strsplit(temp, split="_(?=[^_]*$)", perl = TRUE)[[1]]
      model_name <- temp[2]
      cell_name <- temp[1]
      
      model_uni_fit <- model.temp$lmerel.fit %>% 
        mutate(cell=cell_name, .before="genotype") %>% 
        mutate(model=model_name, .before="gene") %>% 
        bind_rows(model_uni_fit)
    }
  }
  return(model_uni_fit)                                          
}
