library(dplyr); library(data.table); library(readr)
args <- commandArgs(trailingOnly = FALSE)
scaf_nr <- args[[1]]
allgerp <- args[[3]]
window_file <- args[[3]]
out <- args[[4]]

### read in the bed ####
window <- fread(window_file)

for (i in 1:nrow(window)){
  # read in the window and name the chr, start and end pos
  chr = window$V1[i]
  start = window$V2[i]
  end = window$V3[i]
   
  # # extract the gerp mutations from the vcf (in vcf format) within each window
  # tmp_gerps = paste0("output/gerp/windows/gerps_10k_", scaf_nr, "window_", i, "_tmp.vcf")
  # system(paste0("vcftools --vcf ", vcf, " --bed ", allgerp, " --out output/load/gerp/gerps_over4_correct.vcf --recode --keep-INFO-all"))
   
  # now further restrict within each window
  tmp_gerps_window = paste0("output/gerp/windows/gerps_10k_scaf_", scaf_nr, "_window_", i, "_tmp.vcf")
  system(paste0("bcftools view ", allgerp, " --regions '", chr, "':", start, "-", end, " > ", tmp_gerps_window))
  
}

system(paste0("tabix -p vcf ", gsub(".tsv.gz", ".vcf", tsv)))
system(paste0("bedtools intersect -a ", tmp_onewindow, " -b ", tsv, " > ", tmp_gerps))

# export just the locations of the gerp mutations on this scaffold
system(paste0("zcat ", tsv, " | cut > ", gsub(".tsv.gz", "_locations.tsv", tsv)))
# extract the gerp locations from the full vcf
system(paste0("zcat ", tsv, " > ", gsub(".tsv.gz", ".vcf", tsv)))
