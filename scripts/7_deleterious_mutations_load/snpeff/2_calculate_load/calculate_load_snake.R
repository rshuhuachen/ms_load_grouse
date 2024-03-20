#### packages #####
library(dplyr); library(data.table); library(readr)
args <- commandArgs(trailingOnly = TRUE)
snpeff_vcf <- args[[1]]
load_out <- args[[2]]

calculate_load <- function(snpeff_vcf){
  ## metadata on filenames and ids
  filenames <- fread("/vol/cluster-data/rchen/git/genetic_load_ltet/data/metadata/ltet_snps_noindel_minDP20_maxmis0.7_maxdp60_minQC30_pruned_hwe.fam")
  ids <- fread("/vol/cluster-data/rchen/git/genetic_load_ltet/data/metadata/file_list_all_bgi_clean.csv")
  
  #merge
  idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
  
  file <- read.table(snpeff_vcf)
  names(file) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns
  
  #if filtering for only the 30 biggest scaffolds:
  scafs <- fread("data/genomic/refgenome/30_largest.scafs_reduced.tsv")
  scafs$V1 <- gsub(":", ";", scafs$V1)
  scafs$V1 <- gsub("\\.", "=", scafs$V1)

  file <- subset(file, CHROM %in% scafs$V1)
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(10:ncol(file))
  select_n3 <- function(x){x = substr(x,1,3)}
  file[gt] <- lapply(file[gt], select_n3)
  
    # replace genotype with RS value but separate per zygosity, do per ID
  load <- list()

  for (id in 10:ncol(file)){
    subset_id <- file[,c(1:9, id)]
    subset_id <- subset(subset_id, "CHROM" != "ScEsiA3_16870;HRSCAF=19426") #remove scaf 4
    potential_data <- subset(subset_id, subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    realized_data <- subset(subset_id, subset_id[[10]] == "1/1")
    total_data <- subset(subset_id, subset_id[[10]] == "1/1" | subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    potential_load_sum <- nrow(potential_data)
    realized_load_sum <- nrow(realized_data)
    total_load_sum <- nrow(total_data)
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[10]] == "./."))
    n_total <- nrow(subset_id)
    df <- data.frame(id = colnames(file[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     load_p = potential_load_sum / n_genotyped,
                     load_r = realized_load_sum / n_genotyped,
                     load_t_codom = total_load_sum / n_genotyped)
    load[[id]] <- df
  }
  load <- do.call(rbind.data.frame, load)
  return(load)
  
}
    
snpeff_load <- calculate_load(snpeff_vcf)
write_tsv(snpeff_load, file = load_out, col_names=TRUE)