#install.packages(glmmTMB)
library(dplyr); library(data.table); library(readr); library(glmmTMB)
args <- commandArgs(trailingOnly = TRUE)
allgerp <- args[[1]]
window_file <- args[[2]]
scaf_nr <- args[[3]]
out <- args[[4]]

setwd("/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse")

print(allgerp)
print(window_file)
print(scaf_nr)

### read in the bed ####
tryCatch(window <- fread(file = window_file), error=function(e) NULL)

loads_per_scaf <- list()

for (i in 1:nrow(window)){
  # read in the window and name the chr, start and end pos
  chr = window$V1[i]
  start = window$V2[i]
  end = window$V3[i]
   
  # make directory for each scaf
  system(paste0("if [ ! -d output/gerp/windows/scaf_", scaf_nr, " ]; then mkdir output/gerp/windows/scaf_", scaf_nr, "; fi"))
  # extract gerp mutations within the window
  tmp_gerps_window = paste0("output/gerp/windows/scaf_", scaf_nr, "/gerps_100k_scaf_", scaf_nr, "_window_", i, "_tmp.vcf")
  system(paste0("bcftools view ", allgerp, " --regions '", chr, "':", start, "-", end, " -H > ", tmp_gerps_window))
  
  # load in window
  window_sub_n <- system(paste0("cat ", tmp_gerps_window, " | wc -l"), intern = TRUE)
  window_sub_n <- as.numeric(window_sub_n)
  
  # 
  if(window_sub_n > 0){
    tryCatch(window_sub <- read.table(file = tmp_gerps_window), error=function(e) NULL)
    # load function to calc load
    source("scripts/7_calculate_load/function_calculate_load.R")
    # calculate load, although gerp this is snpeff format
    load_window <- calculate_load_snpeff(vcf = window_sub, output_vcf = FALSE, loadtype = paste0("scaf_", scaf_nr, "_window_", i))
    
    
  }
  
  if(window_sub_n == 0){
    load_window <- data.frame(id = NA,
                              n_total = 0,
                              n_genotyped = 0,
                              het_load = NA,
                              hom_load = NA,
                              total_load = NA,
                              loadtype = paste0("scaf_", scaf_nr, "_window_", i))
    # summary_total <- data.frame(scaf_nr = scaf_nr,
    #                             start = start,
    #                             end = end,
    #                             window_nr = i,
    #                             n_snps_window = nrow(window),
    #                             est = NA,
    #                             se = NA,
    #                             zval = NA,
    #                             pval = NA,
    #                             load= "Total")
    }
  
  # output
  #summary_scaf <- rbind(summary_scaf, summary_total)
  loads_per_scaf[[i]] <- load_window
  
  
}

save(loads_per_scaf, file = out)
