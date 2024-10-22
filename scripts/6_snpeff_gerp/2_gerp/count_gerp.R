#### packages #####
library(dplyr); library(data.table); library(readr)
args <- commandArgs(trailingOnly = TRUE)
gerp_vcf <- args[[1]]
scafno <- args[[2]]
gerp_out <- args[[3]]

calculate_gerp_load <- function(gerp_vcf, scafno){
  ## metadata on filenames and ids
  filenames <- fread("data/genomic/raw/metadata/idnames.fam")
  ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
  
  #merge
  idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
  
  file <- read_tsv(gerp_vcf, col_names = c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) )#rename columns
  
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(11:ncol(file))
  select_n3 <- function(x){x = substr(x,1,3)}
  file[gt] <- lapply(file[gt], select_n3)
  
  # replace genotype with RS value but separate per zygosity, do per ID
  gerp_load <- list()
  for( id in 11:ncol(file)){
    subset_id <- file[,c(1:10, id)]
    subset_id <- subset_id %>% mutate(gerp_cat = as.factor(case_when(
        rs_score < 0 ~ "< 0", #changed
        rs_score >= 0 & rs_score < 1 ~ "0-1",
        rs_score >= 1 & rs_score < 2 ~ "1-2",
        rs_score >= 2 & rs_score < 3 ~ "2-3",
        rs_score >= 3 & rs_score < 4 ~ "3-4",
        rs_score >= 4 ~ "4-5"
    )))
    gerp_load_id <- list()
    for (i in c("< 0", "0-1", "1-2", "2-3", "3-4", "4-5")){#changed
        cat_subset <- subset(subset_id, gerp_cat == i)
        het_data <- subset(cat_subset, cat_subset[[11]] == "1/0" | cat_subset[[11]] == "0/1")
        hom_data <- subset(cat_subset, cat_subset[[11]] == "1/1")
        n_genotyped <- nrow(cat_subset) - nrow(subset(cat_subset, cat_subset[[11]] == "./."))
        n_total <- nrow(cat_subset)
        df <- data.frame(id = colnames(file[id]),
                         gerp_cat = i,
                         scafno = scafno,
                         n_total = n_total,
                         n_genotyped = n_genotyped,
                         het_data = nrow(het_data),
                         hom_data = nrow(hom_data))
        
        gerp_load_id[[i]] <- df
        
        }
        gerp_load_id <- do.call(rbind.data.frame, gerp_load_id)
        rownames(gerp_load_id) <- NULL
    
    gerp_load[[id]] <- gerp_load_id
    }
  gerp_load <- do.call(rbind.data.frame, gerp_load)

    return(gerp_load)}


gerp_load <- calculate_gerp_load(gerp_vcf, scafno)
write_tsv(gerp_load, file = gerp_out, col_names=TRUE)