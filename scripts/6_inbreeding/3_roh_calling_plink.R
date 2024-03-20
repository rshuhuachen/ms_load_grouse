##### ROH calling #####

# In this script, we will 
# 1. Estimate SNP density and number
# 2. Run ROH calling

### Libraries
pacman::p_load(tidyverse, data.table)

#### 1. Estimate SNP density and number ####

VCF= paste0(getwd(), "data/genomic/intermediate/ltet_snps_filtered.vcf")

system(paste0("vcftools --vcf ", VCF, " --SNPdensity 10000 --out output/genotyping/genome_density"))
dens <- fread("data/genomes/genome_density.snpden")

#  mean variants/kb = 7.041, min variants/kb = 0, max variants/kb = 39.900

#### 2. Run ROH ####

source("scripts/6_inbreeding/function_run_rohcalling.R")
# this function runs ROH using plink and visualizes the ROHs using a circoplot, plus calculates FROH
# please see the source code of this function for further details
# this function was created to allow easy exploration of the different parameters settings

## Note: VCF files are filtered but no LD pruning of MAF filtering (Meyermans et al 2020)

#### Best settings after experimenting -> this parameter setting combination was used for the manuscript and is saved in the repository ####

run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 125, 
        choose_homozyg_window_het = 5, 
        choose_homozyg_window_missing = 5,
        choose_homozyg_window_threshold = 0.032, 
        choose_homozyg_gap = 1000,
        choose_homozyg_kb = 1000,
        choose_homozyg_density = 50, 
        choose_homozyg_snp = 100,
        choose_homozyg_het = 100, 
        runname = "plink_froh")


### Other parameter settings:
#### Setting 0: best settings after experimenting with adjustment to bcftools ####

run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 125, 
        choose_homozyg_window_het = 5, 
        choose_homozyg_window_missing = 5,
        choose_homozyg_window_threshold = 0.032, 
        choose_homozyg_gap = 1000,
        choose_homozyg_kb = 100,
        choose_homozyg_density = 50, 
        choose_homozyg_snp = 100,
        choose_homozyg_het = 100, 
        runname = paste0("bestrun_allids_adjusted", format(Sys.time(), "_%d_%m_%y")))

#### Setting 1: default settings ####
run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 100, 
        choose_homozyg_window_het = 5, 
        choose_homozyg_window_missing = 5, 
        choose_homozyg_window_threshold = 0.05, 
        choose_homozyg_gap = 1000,
        choose_homozyg_kb = 1000,
        choose_homozyg_density = 50,
        choose_homozyg_snp = 100,
        choose_homozyg_het = 100, 
        runname = paste0("default_allids", format(Sys.time(), "_%d_%m_%y")))



#### Setting 3: relaxed settings ####
# https://doi.org/10.1038/s41467-021-23222-9
run_roh(vcf = VCF,  
        choose_homozyg_window_snp = 125, 
        choose_homozyg_window_het = 5, 
        choose_homozyg_window_missing = 10, #set to missing data % empirical
        choose_homozyg_window_threshold = 0.05, 
        choose_homozyg_gap = 1000,
        choose_homozyg_kb = 1000,
        choose_homozyg_density = 50, 
        choose_homozyg_snp = 100,
        choose_homozyg_het = 100, 
        runname = paste0("longroh_relaxed_allids", format(Sys.time(), "_%d_%m_%y")))

#### Setting 4: stoffel et al 2021 settings ####

run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 50, 
        choose_homozyg_window_het = 2, 
        choose_homozyg_window_missing = 2, 
        choose_homozyg_window_threshold = 0.05, 
        choose_homozyg_gap = 300,
        choose_homozyg_kb = 1200,
        choose_homozyg_density = 200, 
        choose_homozyg_snp = 50,
        choose_homozyg_het = 2, 
        runname = paste0("stoffelpars_allids", format(Sys.time(), "_%d_%m_%y")))

#### Setting 5: relaxed stoffel et al 2021 settings ####

run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 50, 
        choose_homozyg_window_het = 5, 
        choose_homozyg_window_missing = 5, 
        choose_homozyg_window_threshold = 0.05, 
        choose_homozyg_gap = 300,
        choose_homozyg_kb = 1200,
        choose_homozyg_density = 200, 
        choose_homozyg_snp = 50,
        choose_homozyg_het = 2, 
        runname = paste0("stoffelpars_relaxed_allids", format(Sys.time(), "_%d_%m_%y")))


#### Setting 6: cockerill et al 2022 settings ####
# https://doi.org/10.3390/genes13112124
run_roh(vcf = VCF, 
        choose_homozyg_window_snp = 100,
        choose_homozyg_window_het = 3, 
        choose_homozyg_window_missing = 5, 
        choose_homozyg_window_threshold = 0.05, 
        choose_homozyg_gap = 1000,
        choose_homozyg_kb = 100,
        choose_homozyg_density = 50,  
        choose_homozyg_snp = 100,
        choose_homozyg_het = 100, 
        runname = paste0("cockerill_allids", format(Sys.time(), "_%d_%m_%y")))

