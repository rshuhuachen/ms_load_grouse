## bcftools roh funtion was used to estimate ROH using the script in file scripts/6_inbreeding/1_run_mcmcroh.sh
pacman::p_load(data.table, dplyr, ggplot2)

### Run BCFtools
# Set the filename of the output of bcftools
ROH_OUT = "output/inbreeding/mcmc_roh.txt"

# Filter bcftools raw output file to only take ROHs
system("grep \"^RG" , ROH_OUT, "> output/inbreeding/mcmc_roh_rg.txt")

# Read in output
roh <- fread("output/inbreeding/mcmc_roh_rg.txt", 
             col.names=c("state", "file_id", "chr", "start", "end", "length", "nsnp", "qual"))

# Summarise
summary(roh$qual)

# Add real id
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
roh <- left_join(roh, ids[,c(2, 9)], by = c("file_id" = "loc"))

# Filter by minimum lenth and quality
roh <- subset(roh, qual > 30 & length >= 100000 & nsnp >= 100) 

write_tsv(roh, file="output/inbreeding/mcmc_roh_rg_filtered.txt", col.names=T)

# Calculate FROH
roh_bp_perid <- roh %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_bp_perid$froh <- roh_bp_perid$total_length_bp/1004266063
summary(roh_bp_perid$froh)
roh_n_perid <- roh %>% group_by(id) %>% count()
summary(roh_n_perid$n)

# Write out froh
froh_mcmc <- roh_bp_perid[,c("id", "total_length_bp", "froh")]
names(froh_mcmc)[3] <- "froh_mcmc"
write.csv(froh_mcmc, "output/inbreeding/bcftools_mcmc_froh.csv", row.names=F, quote=F)

### Check correlation with FROH based on only 29 autosomal scafs 
scaf <- fread("data/genomic/refgenome/30_largest.scafs.tsv")
scaf <- subset(scaf, scaf_no != 4) #exclude sex scaffold
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

roh_29scaf <- subset(roh1, chr %in% scaf$scaf)
roh_29scaf_bp_perid <- roh_29scaf %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_29scaf_bp_perid$froh <- roh_29scaf_bp_perid$total_length_bp/1004266063
summary(roh_29scaf_bp_perid$froh)
roh_29scaf_n_perid <- roh_29scaf %>% group_by(id) %>% count()
summary(roh_29scaf_n_perid$n)

plot(roh_bp_perid$froh, roh_29scaf_bp_perid$froh) #highly correlated