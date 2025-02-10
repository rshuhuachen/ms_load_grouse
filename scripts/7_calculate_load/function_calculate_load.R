calculate_load_snpeff <- function(vcf, output_vcf, loadtype){
  ## metadata on filenames and ids
  filenames <- fread("data/genomic/raw/metadata/idnames.fam")
  ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
  
  #merge
  idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
  
  names(vcf) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id)#rename columns

  # only select 29 largest autosomal scaffolds
  scaf <- fread("data/genomic/refgenome/30_largest.scafs.tsv")
  scaf$scaf <- gsub(":", ";", scaf$scaf)
  scaf$scaf <- gsub("\\.", "=", scaf$scaf)
  scaf <- subset(scaf, scaf_no != 4)
  
  vcf <- subset(vcf, CHROM %in% scaf$scaf)
  
  # exclude warning messages
  
  vcf <- subset(vcf, !grepl("WARNING", INFO))
  vcf <- as.data.frame(vcf)
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(10:ncol(vcf))
  select_n3 <- function(x){x = substr(x,1,3)}
  vcf[gt] <- lapply(vcf[gt], select_n3)
  
  # calculate load
  load <- list()
  # loop over ids
  for( id in 10:(ncol(vcf))){
    # subset per id
    subset_id <- vcf[,c(1:9, id)]
    
    # filter for snps in het and hom
    het_data <- subset(subset_id, subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    hom_data <- subset(subset_id, subset_id[[10]] == "1/1")
    
    # count amount of snps in het and hom
    het_load_sum <- nrow(het_data)
    hom_load_sum <- nrow(hom_data)
    
    # count no of snps successfully genotyped
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[10]] == "./."))
    n_total <- nrow(subset_id)
    
    # collect data in df
    df <- data.frame(id = colnames(vcf[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     het_load = het_load_sum / n_genotyped,
                     hom_load = hom_load_sum / n_genotyped,
                     total_load = (het_load_sum*0.5 + hom_load_sum) / n_genotyped,
                     loadtype = loadtype)
    load[[id]] <- df
  }
  # convert list to df
  load <- do.call(rbind.data.frame, load)
  
  if(output_vcf == TRUE){
    out <- list(load = load, vcf = vcf)
    return(out)
  }
  
  if(output_vcf==FALSE){
    return(load)}
}

calculate_load_gerp <- function(vcf, output_vcf, loadtype){
  
  ## metadata on filenames and ids
  filenames <- fread("data/genomic/raw/metadata/idnames.fam")
  ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
  
  #merge
  idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
  
  names(vcf) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns
  
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(11:ncol(vcf))
  select_n3 <- function(x){x = substr(x,1,3)}
  vcf[gt] <- lapply(vcf[gt], select_n3)
  
  # calculate load
  load <- list()
  # loop over ids
  for( id in 11:(ncol(vcf))){
    # subset per id
    subset_id <- vcf[,c(1:10, id)]
    
    # filter for snps in het and hom
    het_data <- subset(subset_id, subset_id[[11]] == "1/0" | subset_id[[11]] == "0/1")
    hom_data <- subset(subset_id, subset_id[[11]] == "1/1")
    
    # count amount of snps in het and hom
    het_load_sum <- nrow(het_data)
    hom_load_sum <- nrow(hom_data)
    
    # count no of snps successfully genotyped
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[11]] == "./."))
    n_total <- nrow(subset_id)
    
    # collect data in df
    df <- data.frame(id = colnames(vcf[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     het_load = het_load_sum / n_genotyped,
                     hom_load = hom_load_sum / n_genotyped,
                     total_load = (het_load_sum*0.5 + hom_load_sum) / n_genotyped,
                     loadtype = loadtype)
    load[[id]] <- df
  }
  # convert list to df
  load <- do.call(rbind.data.frame, load)
  
  if(output_vcf == TRUE){
    out <- list(load = load, vcf = vcf)
    return(out)
  }
  
  if(output_vcf==FALSE){
  return(load)}
}


model_load <- function(load, ndraw){
  #join with phenotypes
  load("data/phenotypes/phenotypes_lifetime.RData") #LMS
  pheno <- pheno_wide %>% mutate(core = as.factor(case_when(is.na(LMS) ~ "no core", !is.na(LMS) ~ "core")))
  
  data_pheno <- left_join(load, pheno[,c("id", "LMS", "LMS_min", "core", "site", "born")], by = "id")
  
  model <- glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), family = "poisson", ziformula = ~1, data = data_pheno)
  summary <- summary(model)
  
  beta <- summary$coefficients$cond["scale(total_load)","Estimate"]
  se <- summary$coefficients$cond["scale(total_load)","Std. Error"]
  zval <- summary$coefficients$cond["scale(total_load)","z value"]
  pval <- summary$coefficients$cond["scale(total_load)","Pr(>|z|)"]
  
  out <- data.frame(ndraw = ndraw, 
                    beta = beta,
                    se = se,
                    zval = zval,
                    pval = pval)
  
  return(out)
}
