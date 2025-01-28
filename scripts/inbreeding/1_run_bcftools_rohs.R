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

#add real id
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

roh_raw <- left_join(roh_raw, ids[,c(2, 9)], by = c("file_id" = "loc"))
write.csv(roh_raw, "output/inbreeding/bcftools_roh_rg_clean.txt", quote=F, row.names=F)

# different filtering thresholds
roh1 <- subset(roh, qual > 30 & length >= 100000 & nsnp >= 100) #reassess?


roh1_bp_perid <- roh1 %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh1_bp_perid$froh <- roh1_bp_perid$total_length_bp/1004266063
summary(roh1_bp_perid$froh)
roh1_n_perid <- roh1 %>% group_by(id) %>% count()
summary(roh1_n_perid$n)

roh2_bp_perid <- roh2 %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh2_bp_perid$froh <- roh2_bp_perid$total_length_bp/1004266063
summary(roh2_bp_perid$froh)
roh2_n_perid <- roh2 %>% group_by(id) %>% count()
summary(roh2_n_perid$n)

roh3_bp_perid <- roh3 %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh3_bp_perid$froh <- roh3_bp_perid$total_length_bp/1004266063
summary(roh3_bp_perid$froh)
roh3_n_perid <- roh3 %>% group_by(id) %>% count()
summary(roh3_n_perid$n)

plot(roh1_bp_perid$froh, roh3_bp_perid$froh)
plot(roh2_bp_perid$froh, roh3_bp_perid$froh)
plot(roh1_bp_perid$froh, roh2_bp_perid$froh)

#### choose option 3
roh <- subset(roh, length >= 50000 & nsnp >= 50)
summary(roh$nsnp)

write.csv(roh, "output/inbreeding/roh_bcftools/roh_rg_clean_filtered.txt", quote=F, row.names=F)

## summarise
roh_bp_perid <- roh %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_n_perid <- roh %>% group_by(id) %>% count()

roh_snp_perid <- roh %>% group_by(id) %>% summarise(total_length_snp = sum(nsnp))
roh_mean_snp_perid <- roh %>% group_by(id) %>% summarise(mean_length_snp = mean(nsnp))

# summary

summary(roh_bp_perid$total_length_bp)
summary(roh_snp_perid$total_length_snp)
summary(roh_mean_snp_perid$mean_length_snp)
summary(roh_n_perid$n)

# froh
roh_bp_perid$froh <- roh_bp_perid$total_length_bp/1004266063
summary(roh_bp_perid$froh)

#write out froh
froh_mcmc <- roh_bp_perid[,c("id", "total_length_bp", "froh")]
names(froh_mcmc)[3] <- "froh_mcmc"
write.csv(froh_mcmc, "output/inbreeding/bcftools_mcmc_froh.csv", row.names=F, quote=F)

##### Only 29 autosomal scafs #####
scaf <- fread("/vol/cluster-data/rchen/geneticload/gerp/analyses/git/cactus_insert_ltet_take4/data/genomes/30_largest.scafs_nosex.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

roh_29scaf <- subset(roh1, chr %in% scaf$scaf)
roh_29scaf_bp_perid <- roh_29scaf %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_29scaf_bp_perid$froh <- roh_29scaf_bp_perid$total_length_bp/1004266063
summary(roh_29scaf_bp_perid$froh)
roh_29scaf_n_perid <- roh_29scaf %>% group_by(id) %>% count()
summary(roh_29scaf_n_perid$n)

#relationship with all scaf FROH
plot(roh1_bp_perid$froh, roh_29scaf_bp_perid$froh) #high correlation

load("/vol/cluster-data/rchen/git/genetic_load_ltet/data/phenotypes/phenotypes_wide_extra.RData") #LMS
pheno <- pheno %>% mutate(core = as.factor(case_when(is.na(LMS) ~ "no core", !is.na(LMS) ~ "core")))

froh_pheno <- left_join(roh_29scaf_bp_perid[,c("id", "froh")], pheno[,c("id", "LMS_min", "core", "site", "born")], by = "id")
summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + core + (1|born) + (1|site) , data = froh_pheno, family = "poisson", ziformula = ~1))

#### Take out sex scafs #####
scaf <- fread("/vol/cluster-data/rchen/geneticload/gerp/analyses/git/cactus_insert_ltet_take4/data/genomes/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

roh_autoscaf <- subset(roh1, chr != scaf$scaf[which(scaf$scaf_no==4)])
roh_autoscaf_bp_perid <- roh_autoscaf %>% group_by(id) %>% summarise(total_length_bp = sum(length))
roh_autoscaf_bp_perid$froh <- roh_autoscaf_bp_perid$total_length_bp/(1003452484-scaf$size_base[which(scaf$scaf_no==4)])
summary(roh_autoscaf_bp_perid$froh)
roh_autoscaf_n_perid <- roh_autoscaf %>% group_by(id) %>% count()
summary(roh_autoscaf_n_perid$n)

#relationship with all scaf FROH
plot(roh1_bp_perid$froh, roh_autoscaf_bp_perid$froh) #high correlation

load("/vol/cluster-data/rchen/git/genetic_load_ltet/data/phenotypes/phenotypes_wide_extra.RData") #LMS
pheno <- pheno %>% mutate(core = as.factor(case_when(is.na(LMS) ~ "no core", !is.na(LMS) ~ "core")))

froh_pheno <- left_join(roh_autoscaf_bp_perid[,c("id", "froh")], pheno[,c("id", "LMS_min", "core", "site", "born")], by = "id")
summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + core + (1|born) + (1|site) , data = froh_pheno, family = "poisson", ziformula = ~1))

roh_autoscaf_bp_perid <- roh_autoscaf_bp_perid %>% rename(froh_auto = froh)

save(roh_autoscaf_bp_perid, file = "output/inbreeding/froh_auto.RData")
save(roh_autoscaf, file = "output/inbreeding/rohs_auto.RData")
