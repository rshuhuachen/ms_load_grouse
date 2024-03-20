### does gerp based on coding regions only affect LMS? ###

# packages
pacman::p_load(tidyverse, data.table)

# load annotated gerp data
load(file = "output/load/gerp/gerp_annotated_region.RData")
gerp_all$CHROM <- gsub("__", ";", gerp_all$CHROM)
gerp_all$CHROM <- gsub("HRSCAF_", "HRSCAF=", gerp_all$CHROM)

# load annotation
load(file = "output/load/snpeff/snpeff_high_annotated_region.RData")
snpef_all$CHROM <- gsub("__", ";", snpef_all$CHROM)
snpef_all$CHROM <- gsub("HRSCAF_", "HRSCAF=", snpef_all$CHROM)

## load functions to calculate load
source("scripts/7_deleterious_mutations_load/per_region/functions_snpefloads.R")
source("scripts/7_deleterious_mutations_load/per_region/functions_gerploads.R")

#### calculate load per region as defined below 

regions <- c("region_promoter", "region_gene","region_tss" ,"region_exon","region_intron" ,"region_down" ,
             "region_up")

#load existing combined load file
load("output/load/all_loads_combined_da_nosex_29scaf.RData") #loads no sex chr only 30 scaf

#first gerp5
load_per_region <- loads
for (region in regions){
  subset_locs <- subset(gerp_all, gerp_all[,region] == 1) #subset based on region name
  
  subset_locs$chr_pos <- paste0(subset_locs$CHROM, subset_locs$POS, sep = "_") #make a col for the snp position
  
  sub_genotypes <- subset(gerps_cat5_vcf, chr_pos %in% subset_locs$chr_pos) #subset genotypes based on subset
  
  sub_genotypes$chr_pos <- NULL #remove chr_pos again
  
  load_sub <- calculate_load_gerp(sub_genotypes) #calculate loads
  load_sub$n_genotyped <- NULL #don't need this anymore
  
  region_name <- gsub("region_", "", region)
  names(load_sub)[2] <- paste0("gerp5_n_total_", region_name) #add the region in the column name for easier merge
  names(load_sub)[3] <- paste0("gerp5_load_m_", region_name)
  names(load_sub)[4] <- paste0("gerp5_load_e_", region_name)
  names(load_sub)[5] <- paste0("gerp5_load_t_add_", region_name)
  
  load_per_region <- left_join(load_per_region, load_sub, by = "id")
}

#snpeff
for (region in regions){
  subset_locs <- subset(snpef_all, snpef_all[,region] == 1) #subset based on region name
  
  subset_locs$chr_pos <- paste0(subset_locs$CHROM, subset_locs$POS, sep = "_") #make a col for the snp position
  
  sub_genotypes <- subset(snpef, chr_pos %in% subset_locs$chr_pos) #subset genotypes based on subset
  
  sub_genotypes$chr_pos <- NULL #remove chr_pos again
  
  load_sub <- calculate_load_snpef(sub_genotypes) #calculate loads
  load_sub$n_genotyped <- NULL #don't need this anymore
  
  region_name <- gsub("region_", "", region)
  names(load_sub)[2] <- paste0("high_n_total_", region_name) #add the region in the column name for easier merge
  names(load_sub)[3] <- paste0("high_load_m_", region_name)
  names(load_sub)[4] <- paste0("high_load_e_", region_name)
  names(load_sub)[5] <- paste0("high_load_t_add_", region_name)
  
  load_per_region <- left_join(load_per_region, load_sub, by = "id")
}

save(load_per_region, file = "output/load/all_loads_combined_da_nosex_30scaf_plus_per_region.RData")
