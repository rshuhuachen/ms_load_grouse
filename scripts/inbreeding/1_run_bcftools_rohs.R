## bcftools roh funtion was used to estimate ROH using the script in fie /vol/cluster-data/rchen/wgr/analyses/rohcalling/mcmc_roh/run_mcmcroh.sh
pacman::p_load(data.table, dplyr, ggplot2)

### run BCFtools
### Set the directories of both your pre-processed VCF file as well as the absolute path of the github repository 
rawvcfrepo = "data/genomic/intermediate/" # Directory of the pre-processed VCF file (output of script inbreeding-sexual-traits/scripts/1_genotyping/snpfiltering/snpfilter.sh)
VCF= paste0(rawvcfrepo, "ltet_snps_filtered.vcf")
ROH_OUT = paste0(getwd(), "/output/inbreeding/bcftools_roh.txt")

system(paste0("bcftools roh -G30 --AF-dflt 0.4 ", VCF, " -o ", ROH_OUT))

system("grep \"^RG " , ROH_OUT, "> ", getwd(), "output/inbreeding/bcftools_roh_rg.txt") #select info on rohs only and exclude info on individual sites

roh_raw <- fread("output/inbreeding/bcftools_roh_rg.txt", 
             col.names=c("state", "file_id", "chr", "start", "end", "length", "nsnp", "qual"))

# add real id
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

roh_raw <- left_join(roh_raw, ids[,c(2, 9)], by = c("file_id" = "loc"))

# take out sex scaf
roh_autoscaf <- subset(roh_raw, chr != "ScEsiA3_16870;HRSCAF=19426")

summary(roh_autoscaf$length)
summary(roh_autoscaf$nsnp)

# different filtering thresholds
roh_clean <- subset(roh_autoscaf, qual > 30 & length >= 10000 & nsnp >= 100) 

# get froh
froh <- roh_clean %>% group_by(id) %>% summarise(total_length_bp = sum(length))
froh$froh <- froh$total_length_bp/1004266063
summary(froh$froh)

n_roh <- roh_clean %>% group_by(id) %>% count()
summary(n_roh$n)

ggplot(roh_clean, aes(x = length)) + geom_histogram() + scale_y_log10()
ggplot(froh, aes(x = froh)) + geom_histogram()

# export
save(froh, file = "output/inbreeding/froh.RData")
write.csv(roh_clean, "output/inbreeding/bcftools_roh_rg_clean.txt", quote=F, row.names=F)

# relationship with all scaf FROH
load("data/phenotypes/phenotypes_lifetime.RData") #LMS

froh_pheno <- left_join(froh[,c("id", "froh")], pheno_wide[,c("id", "LMS_min", "core", "site", "born")], by = "id")
summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + core + (1|site) , data = froh_pheno, family = "poisson", ziformula = ~1))


