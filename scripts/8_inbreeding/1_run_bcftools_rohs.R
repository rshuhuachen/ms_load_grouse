### load packages ####
pacman::p_load(data.table, dplyr, ggplot2, glmmTMB)

### bcftools ####
# to estimate Froh, we identify ROHs using BCFtools 
# within the R script, we run BCFtools through the shell (not very computationally intensive)

# set the parameters of where the vcf file can be found
rawvcfrepo = "data/genomic/intermediate/" 
VCF= paste0(rawvcfrepo, "ltet_snps_filtered.vcf")
ROH_OUT = paste0(getwd(), "/output/inbreeding/bcftools_roh.txt")

# run bcftools
system(paste0("bcftools roh -G30 --AF-dflt 0.4 ", VCF, " -o ", ROH_OUT))

# bcftools outputs lots of information, but we only care about the ROHs which are the lines that start with "RG"
system("grep \"^RG " , ROH_OUT, "> ", getwd(), "output/inbreeding/bcftools_roh_rg.txt") #select info on rohs only and exclude info on individual sites

#### Process output ####

# now, we can import the bcftools output 
roh_raw <- fread("output/inbreeding/bcftools_roh_rg.txt", 
             col.names=c("state", "file_id", "chr", "start", "end", "length", "nsnp", "qual"))

# the sample IDs currently refer to long file names, we can add the real ID for simplicity
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
roh_raw <- left_join(roh_raw, ids[,c(2, 9)], by = c("file_id" = "loc"))

# Froh is based on autosomes only, so here we take out the sex-linked scaffold
roh_autoscaf <- subset(roh_raw, chr != "ScEsiA3_16870;HRSCAF=19426")

# BCFtools does not have any filters or quality thresholds, meaning that a ROH can be defined as just 2 SNPs. 
# Here, we set some thresholds both for quality (minimum of 30), ROH length (minimum of 10kb) and number of SNPs within the ROH (minimum of 100) 
roh_clean <- subset(roh_autoscaf, qual > 30 & length >= 10000 & nsnp >= 100) 

# To now calculate Froh, we sum the total length of all ROHs for all individuals and divide it by the total autosomal genome length (which is 1004266063)
froh <- roh_clean %>% group_by(id) %>% summarise(total_length_bp = sum(length))
froh$froh <- froh$total_length_bp/1004266063

# We can also calculate the number of ROHs per individual, regardless of their length
n_roh <- roh_clean %>% group_by(id) %>% count()

# Quick visualisation
ggplot(roh_clean, aes(x = length)) + geom_histogram() + scale_y_log10()
ggplot(froh, aes(x = froh)) + geom_histogram()

# Now we export the Froh numbers for further analysis
save(froh, file = "output/inbreeding/froh.RData")
write.csv(roh_clean, "output/inbreeding/bcftools_roh_rg_clean.txt", quote=F, row.names=F)

# For a quick model to estimate the relationship between Froh and LMS, we check it using a frequentist framework here
load("data/phenotypes/phenotypes_lifetime.RData") # phenotypic and LMS data

froh_pheno <- left_join(froh[,c("id", "froh")], pheno_wide[,c("id", "LMS_min", "core", "site", "born")], by = "id")
summary(glmmTMB(LMS_min ~ scale(froh) + core + (1|site) , data = froh_pheno, family = "poisson", ziformula = ~1))


