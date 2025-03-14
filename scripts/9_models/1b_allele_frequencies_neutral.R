#### calculated allele frequencies of low impact mutations and GERP < 0 (both expected to be neutral) ####
### here: load in, assign ancestral and derived alleles and calculate mean ####

pacman::p_load(tidyverse, data.table)

### First: get RData frame with just neutral GERP snps

gerp_snp_scafs <- list.files(path = "output/gerp/beds", pattern = "gerp_overlapSNP*", full.names = T)
gerp_snp_scafs <- gerp_snp_scafs[-22] #empty

gerp_snp <- data.frame()
for (i in 1:length(gerp_snp_scafs)){
  scaf <- read.table(gerp_snp_scafs[i])
  scaf <- scaf %>% filter(V5 < 0)
  gerp_snp <- rbind(gerp_snp, scaf)
}

## metadata on filenames and ids
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(gerp_snp) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns

save(gerp_snp, file = "output/load/gerp/gerps_neutral.RData")
# reformat gerp to .vcf
# merge with clean ID names

gerps_vcf <- gerp_snp %>% select(c(chr, pos, ancestral, derived, qual, info, format, D154259:D229777)) #exclude gerp scores
gerps_vcf <- gerps_vcf %>% mutate(ID = ".", .after = pos) %>% mutate(FILTER = ".", .after = qual)
names(gerps_vcf)[1:9] <- c("#CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT") #same names as .vcf

write_tsv(gerps_vcf, file = "output/load/gerp/gerps_neutral.vcf", col_names = TRUE)

### calculate allele frequencies with vcftools
system("vcftools --gzvcf data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf.gz --freq --out output/af/allele_freq_low")

system("vcftools --vcf output/load/gerp/gerps_neutral.vcf --freq --out output/af/allele_freq_gerp0")

##
vcf_low <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf.gz")
vcf_gerp<-gerps_vcf
#vcf_gerp <- read.table("output/load/gerp/gerps_neutral.vcf")

allele_freq_low <- fread("output/af/allele_freq_low.frq", skip = 1)
allele_freq_gerp <- fread("output/af/allele_freq_gerp0.frq", skip = 1)

names(allele_freq_low) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")
names(allele_freq_gerp) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")

#### make long with one row one allele freq
allele_freq_low_long <- gather(allele_freq_low, allele, allele_freq_low, allele_freq_1:allele_freq_2, factor_key=TRUE)
allele_freq_gerp0_long <- gather(allele_freq_gerp, allele, allele_freq_gerp, allele_freq_1:allele_freq_2, factor_key=TRUE)

### merge derived allele
allele_freq_low_long <- left_join(allele_freq_low_long, vcf_low[,c("V1", "V2", "V5")], by = c("scaf" = "V1", "pos" = "V2")) #chr, pos, derived
allele_freq_gerp0_long <- left_join(allele_freq_gerp0_long, vcf_gerp[,c("#CHROM", "POS", "ALT")], by = c("scaf" = "#CHROM", "pos" = "POS")) #chr, pos, derived

## subset where allele = derived
allele_freq_low_long_derived <- subset(allele_freq_low_long, substr(allele_freq_low,1,1) == V5)
allele_freq_gerp0_long_derived <- subset(allele_freq_gerp0_long, substr(allele_freq_gerp,1,1) == ALT)

## only take frequency
allele_freq_low_long_derived$frequency <- as.numeric(substr(allele_freq_low_long_derived$allele_freq_low, 3, nchar(allele_freq_low_long_derived$allele_freq_low)))
allele_freq_gerp0_long_derived$frequency <- as.numeric(substr(allele_freq_gerp0_long_derived$allele_freq_gerp, 3, nchar(allele_freq_gerp0_long_derived$allele_freq_gerp)))

## summary allele freq
summary(allele_freq_low_long_derived$frequency)
summary(allele_freq_gerp0_long_derived$frequency)

## how many are fixed?
nrow(subset(allele_freq_low_long_derived, frequency == 1))
nrow(subset(allele_freq_gerp0_long_derived, frequency == 1))

save(allele_freq_low_long_derived, file = "output/load/snpeff/allelefreq_snpeff_low.derived.RData")
save(allele_freq_gerp0_long_derived, file = "output/load/gerp/allelefreq_gerp0.derived.RData")

## histograms
source("scripts/theme_ggplot.R")

ggplot(allele_freq_low_long_derived, aes(x = frequency)) + geom_histogram(fill = clrs[5], color=clrs[1]) + 
  labs(x = "Allele frequencies derived low impact mutations SnpEff", y = "Count") -> hist_af_low

ggplot(allele_freq_gerp0_long_derived, aes(x = frequency)) + geom_histogram(fill = clrs[5], color=clrs[1])+ 
  labs(x = "Allele frequencies derived mutations with GERP score < 0", y = "Count")-> hist_af_gerp

