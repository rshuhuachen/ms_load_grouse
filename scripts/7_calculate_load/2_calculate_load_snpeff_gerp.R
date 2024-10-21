#### Here we will calculate mutation load based on SnpEff and GERP #####

### load packages ###
pacman::p_load(dplyr, data.table)

### load function to calculate load ###
source("scripts/7_calculate_load/function_calculate_load.R")

##### snpeff categories #####

high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz")
moderate <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf.gz")
low <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf.gz")
lof <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_LOF.vcf.gz")
missense <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_missense.vcf.gz")

## calculate load for each
# in this function, we give the columns names, filter for only the largest 29 autosomal scaffolds and exclude annotations with warning messages

high_load <- calculate_load_snpeff(high, output_vcf = TRUE, loadtype = "high")
moderate_load <- calculate_load_snpeff(moderate, output_vcf = FALSE, loadtype = "moderate")
low_load <- calculate_load_snpeff(low, output_vcf = FALSE, loadtype = "low")
lof_load <- calculate_load_snpeff(lof, output_vcf = FALSE, loadtype = "lof")
missense_load <- calculate_load_snpeff(missense, output_vcf = FALSE, loadtype = "missense")

### save high impact mutations 
snpeff <- high_load$vcf
save(snpeff, file = "output/load/snpeff/snpeff_high.RData")

##### GERP #####

gerp_scafs <- list.files(path = "output/gerp", pattern = "count_mutations*", full.names = T)

gerp_raw <- data.frame()
for (i in 1:length(gerp_scafs)){
  scaf <- fread(gerp_scafs[i])
  gerp_raw <- rbind(gerp_raw, scaf)
}

gerp <- gerp_raw %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:hom_data, sum))

gerp <- gerp %>% mutate(het_load = het_data / n_genotyped,
                        hom_load = hom_data / n_genotyped,
                        total_load = (het_data * 0.5 + hom_data) / n_genotyped)

gerp_34_load <- subset(gerp, gerp_cat == "3-4")
gerp_45_load <- subset(gerp, gerp_cat == "4-5")

gerp_34_load$loadtype = "gerp34"
gerp_45_load$loadtype = "gerp45"

# extract count per category
gerp_count <- gerp %>% ungroup %>% select(c(gerp_cat, n_total)) %>% unique()
save(gerp_count, file = "output/load/gerp/gerps_count_per_cat.RData")

# for future: also read in all gerp >= 4 for subsetting etc
gerp_snp_scafs <- list.files(path = "output/gerp/beds", pattern = "gerp_overlapSNP*", full.names = T)
gerp_snp_scafs <- gerp_snp_scafs[-22] #empty

gerp_snp <- data.frame()
for (i in 1:length(gerp_snp_scafs)){
  scaf <- read.table(gerp_snp_scafs[i])
  scaf <- scaf %>% filter(V5 >= 4)
  gerp_snp <- rbind(gerp_snp, scaf)
}

gerp_45_load_check <- calculate_load_gerp(gerp_snp, output_vcf = TRUE, loadtype = "gerp45") #413489 = n
check <- left_join(gerp_45_load, gerp_45_load_check$load, by = "id")
View(check[,c("total_load.x", "total_load.y")])
gerp <- gerp_45_load_check$vcf
save(gerp, file = "output/load/gerp/gerp_over4.RData")

# for the time being, get load estimate from genetic_load dir
# 
# load("/vol/cluster-data/rchen/geneticload/gerp/analyses/git/cactus_insert_ltet_take4/output/gerp/scores/sum_gerp_cat_with0.RData")
# gerp_34_load <- load_cat_select %>% select(c(id, gerp_cat, n_total, n_genotyped, gerp_Lp, gerp_Lr, gerp_Lt_additive)) %>% 
#   filter(gerp_cat == "3-4")%>%
#     select(-c(gerp_cat)) 
# gerp_34_load$loadtype = "gerp34"
# names(gerp_34_load) <- names(moderate_load)  
# 
# gerp_45_load <- load_cat_select %>% select(c(id, gerp_cat, n_total, n_genotyped, gerp_Lp, gerp_Lr, gerp_Lt_additive)) %>% 
#   filter(gerp_cat == "4-5")%>%
#   select(-c(gerp_cat)) 
# gerp_45_load$loadtype = "gerp45"
# names(gerp_45_load) <- names(moderate_load)  

#### Save in a single file #####

load <- rbind(high_load$load[,c("id", "het_load", "hom_load", "total_load", "loadtype")] ,
              moderate_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")], 
              low_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")],
              lof_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")],
              missense_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")], 
              gerp_34_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")], 
              gerp_45_load[,c("id", "het_load", "hom_load", "total_load", "loadtype")])

save(load, file = "output/load/all_loads_combined_da_nosex_29scaf.RData")
write.table(load, file = "output/load/all_loads_combined_da_nosex_29scaf.tsv", sep="\t", row.names = F)

cor.test(load$total_load[which(load$loadtype == "gerp45")], load$total_load[which(load$loadtype == "high")])

### Test for lek effects ####
load("data/phenotypes/phenotypes_lifetime.RData")
pheno_load <- left_join(pheno_wide, load, by = "id")
summary(lm(total_load ~ site, data = subset(pheno_load, loadtype == "gerp45")))
summary(lm(total_load ~ site, data = subset(pheno_load, loadtype == "high")))

#### Overlap GERP > 4 and high ####
high_nowarning <- subset(high, !grepl("WARNING", V8))

# only select 29 largest autosomal scaffolds
scaf <- fread("data/genomic/refgenome/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)
scaf <- subset(scaf, scaf_no != 4)

high_nowarning <- subset(high_nowarning, V1 %in% scaf$scaf)
high_nowarning$chr_pos <- paste0(high_nowarning$V1, "_", high_nowarning$V2)
gerp$chr_pos <- paste0(gerp$chr, "_", gerp$pos)
nrow(subset(high_nowarning, chr_pos %in% gerp$chr_pos)) #274
