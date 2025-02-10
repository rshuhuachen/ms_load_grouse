#### Calculate total load per chromosome ####

### packages ####
pacman::p_load(tidyverse, data.table, glmmTMB)

### load GERP data ####
load(file = "output/load/gerp/gerp_over4.RData")

### load SnpEff data ####
high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz")

#### load functions #####
source("scripts/7_calculate_load/function_calculate_load.R")
source("scripts/theme_ggplot.R")

#### GERP ####
# split mutations by scaffold

#### scaf name ####
load("data/genomic/raw/metadata/scaffold_names_dovetail.RData")
gerp <- left_join(gerp, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

#### load per chr ####
list_gerps <- data.frame()
for (i in 1:length(unique(gerp$chr))){
  sub_chr <- subset(gerp, chr == unique(gerp$chr)[i])
  load_chr <- calculate_load_gerp(vcf = sub_chr[,-c(ncol(sub_chr))], output_vcf = F, loadtype = unique(gerp$scaf_nr)[i])
  list_gerps <- rbind(list_gerps, load_chr)
}

#### SnpEff ###
high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz")
## metadata on filenames and ids
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(high) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id)#rename columns

# only select 29 largest autosomal scaffolds
scaf <- fread("data/genomic/refgenome/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)
scaf <- subset(scaf, scaf_no != 4)

high <- subset(high, CHROM %in% scaf$scaf)

high <- left_join(high, genome[,c("contig", "scaf_nr")], by = c("CHROM" = "contig"))

list_high <- data.frame()
for (i in 1:length(unique(high$scaf_nr))){
  sub_chr <- subset(high, scaf_nr == unique(high$scaf_nr)[i])
  load_chr <- calculate_load_snpeff(vcf = sub_chr[,-c(ncol(sub_chr))], output_vcf = F, loadtype = unique(high$scaf_nr)[i])
  list_high <- rbind(list_gerps, load_chr)
}

### merge together
list_gerps$method = "gerp45"
list_high$method = "high"
per_chr_load <- rbind(list_gerps, list_high)
save(per_chr_load, file = "output/load/load_per_chr_snpeff_gerp.RData")

#### merge with pheno data ###
load("data/phenotypes/phenotypes_lifetime.RData")

per_chr_load <- left_join(per_chr_load, pheno_wide, by = "id")


###### BELOW IS OLD ######
# #### loop over chr to get output of model ####
# n_mut <- unique(list_gerps[,c("loadtype", "n_total")])
# summary_models <- data.frame()
# 
# for (i in 1:length(unique(gerp_load$loadtype))){
# 
#   ### total load
#   model_total <- glmmTMB(LMS_min ~ scale(total_load) + core + (1|site),
#                          family = "poisson", ziformula = ~1,
#                          data = subset(gerp_load, loadtype == unique(gerp_load$loadtype)[i]))
#   sum_total <- summary(model_total)
#   coef_total <- as.data.frame(sum_total$coefficients$cond)
# 
#   summary_total <- data.frame(chr = unique(gerp_load$loadtype)[i],
#                               est = coef_total$Estimate[2],
#                               se = coef_total$`Std. Error`[2],
#                               zval = coef_total$`z value`[2],
#                               pval = coef_total$`Pr(>|z|)`[2],
#                               load= "Total")
# 
#   ### hom/het load
#   model_hom <- glmmTMB(LMS_min ~ scale(het_load) + scale(hom_load) + core + (1|site),
#                        family = "poisson", ziformula = ~1,
#                        data = subset(gerp_load, loadtype == unique(gerp_load$loadtype)[i]))
#   sum_hom <- summary(model_hom)
#   coef_hom <- as.data.frame(sum_hom$coefficients$cond)
# 
#   summary_het <- data.frame(chr = unique(gerp_load$loadtype)[i],
#                             est = coef_hom$Estimate[2],
#                             se = coef_hom$`Std. Error`[2],
#                             zval = coef_hom$`z value`[2],
#                             pval = coef_hom$`Pr(>|z|)`[2],
#                             load= "Heterozygous")
# 
#   summary_hom <- data.frame(chr = unique(gerp_load$loadtype)[i],
#                             est = coef_hom$Estimate[3],
#                             se = coef_hom$`Std. Error`[3],
#                             zval = coef_hom$`z value`[3],
#                             pval = coef_hom$`Pr(>|z|)`[3],
#                             load= "Homozygous")
# 
#   summary <- rbind(summary_total, summary_het, summary_hom)
#   summary <- summary %>% mutate(sig = case_when(pval < 0.05 ~ "Significant", TRUE ~ "Non significant"))
#   summary_models <- rbind(summary_models, summary)
# }
# 
# summary_models <- left_join(summary_models, n_mut, by = c("chr" = "loadtype"))
# 
# #### plot ####
# ggplot(summary_models, aes(x = est, y = as.factor(chr), col = sig)) +
#   geom_point(size=4) +
#   geom_segment(aes(x = (est-se), xend = (est+se), yend = as.factor(chr)), col = "black", linewidth=1) +
#   geom_vline(xintercept = 0, col = "darkred", linetype = "dotted") +
#   facet_grid(~load, scales="free")+
#   scale_color_manual(values=c(clr_grey, clr_gerp)) +
#   labs(col = "Significant", x = "Beta estimate load* ", y = "Chromosome") -> gerp_per_chr
# 
# png("plots/gerp_per_chr_lms.png", width=1000, height = 800)
# gerp_per_chr
# dev.off()
# 
# ggplot(subset(summary_models, load == "Total"), aes(x = est, y = as.factor(chr), col = sig)) +
#   geom_point(size=4) +
#   geom_segment(aes(x = (est-se), xend = (est+se), yend = as.factor(chr)), col = "black", linewidth=1) +
#   geom_vline(xintercept = 0, col = "darkred", linetype = "dotted") +
#   scale_color_manual(values=c(clr_grey, clr_gerp)) +
#   labs(col = "Significant", x = "Beta estimate load", y = "Chromosome") -> gerp_per_chr_total
# 
# ## nr of snps per chr ##
# system('cat output/ancestral/ltet_filtered_ann_aa.vcf | grep -v "^#" | cut -f 1 | sort | uniq -c > output/n_snps_per_chr.txt')
# 
# n_snp <- read.table("output/n_snps_per_chr.txt")
# names(n_snp) <- c("n_snp", "contig")
# 
# genome <- left_join(genome, n_snp, by = "contig")
# 
# ggplot(subset(genome[c(1:30),]), aes(x = n_bases, y = n_snp)) + geom_point()
# 
# ## nr of gerp mutations per chr
# n_gerp <- gerp %>% group_by(chr) %>% count()
# genome <- left_join(genome, n_gerp, by = c("contig" = "chr"))
# 
# ggplot(subset(genome[c(1:30),]), aes(x = n_bases, y = n)) + geom_point()
# 
# summary_models <- left_join(summary_models, genome, by = c("chr" = "scaf_nr"))
# 
# save(summary_models, file = "output/load_per_chr/summary_frequentist_load_per_chr_gerp.RData")
