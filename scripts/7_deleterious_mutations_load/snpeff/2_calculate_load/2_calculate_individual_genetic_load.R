### convert to geno
pacman::p_load(tidyverse, data.table)

##### load vcfs #####
high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf")

### Calculate genetic load -- derived allele (da) based ####
source("scripts/7_deleterious_mutations_load/snpeff/2_calculate_load/function_load.R")

load_high_da <- calculate_load(high)

#additive instead of codominant
source("scripts/7_deleterious_mutations_load/snpeff/2_calculate_load/function_load.R")


#### Safe ####  
load_high_da <- load_high_da$load
names(load_high_da) <- c("id", "high_Lp", "high_Lr", "high_Lt_codom")

save(load_high_da, file="output/load/snpeff/genetic_load_da_high.RData")


