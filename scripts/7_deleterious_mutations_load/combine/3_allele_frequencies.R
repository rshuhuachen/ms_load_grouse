#### calculated allele frequencies of high impact mutations and GERP >= 4 ####
### here: load in, assign ancestral and derived alleles and calculate mean ####
### not sure how vcftool assigns frist and second allele (is it major/minor or ancestral/derived?)

pacman::p_load(tidyverse, data.table)

### First: reformat gerp to .vcf

#load in all gerps
gerps_vcf <- data.frame()
for (i in c(1:3,5:30)){
  scaf <- read.table(paste0("output/gerp/beds/gerp_overlapSNP_scaf_",i, ".tsv.gz"))
  gerps_vcf <- rbind(gerps_vcf, scaf)
}

# reformat to vcf
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

# merge with clean ID names
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(gerps_vcf) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns

gerps_vcf_vcf <- gerps_vcf %>% select(c(chr, pos, ancestral, derived, qual, info, format, D154259:D229777)) #exclude gerp scores
gerps_vcf_vcf <- gerps_vcf_vcf %>% mutate(ID = ".", .after = pos) %>% mutate(FILTER = ".", .after = qual)
names(gerps_vcf_vcf)[1:9] <- c("#CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT") #same names as .vcf

write_tsv(gerps_vcf_vcf, file = "output/load/gerp/gerps_cat5_vcfformat.vcf", col_names = TRUE)

### calculate allele frequencies with vcftools
system(paste0("vcftools --vcf data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf --freq --out output/snpeff/allele_freq_high")

system(paste0("vcftools --vcf output/load/gerp/gerps_cat5_vcfformat.vcf --freq --out output/gerp/allele_freq_gerp5")

##
vcf_high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf")
vcf_gerp <- read.table("output/load/gerp/gerps_cat5_vcfformat.vcf")

allele_freq_high <- fread("output/allele_freq_high.frq", skip = 1)
allele_freq_gerp <- fread("output/allele_freq_gerp5.frq", skip = 1)

names(allele_freq_high) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")
names(allele_freq_gerp) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")

#### make long with one row one allele freq
allele_freq_high_long <- gather(allele_freq_high, allele, allele_freq_high, allele_freq_1:allele_freq_2, factor_key=TRUE)
allele_freq_gerp_long <- gather(allele_freq_gerp, allele, allele_freq_gerp, allele_freq_1:allele_freq_2, factor_key=TRUE)

### merge derived allele
allele_freq_high_long <- left_join(allele_freq_high_long, vcf_high[,c("V1", "V2", "V5")], by = c("scaf" = "V1", "pos" = "V2")) #chr, pos, derived
allele_freq_gerp_long <- left_join(allele_freq_gerp_long, vcf_gerp[,c("V1", "V2", "V5")], by = c("scaf" = "V1", "pos" = "V2")) #chr, pos, derived

## subset where allele = derived
allele_freq_high_long_derived <- subset(allele_freq_high_long, substr(allele_freq_high,1,1) == V5)
allele_freq_gerp_long_derived <- subset(allele_freq_gerp_long, substr(allele_freq_gerp,1,1) == V5)

## only take frequency
allele_freq_high_long_derived$frequency <- as.numeric(substr(allele_freq_high_long_derived$allele_freq_high, 3, nchar(allele_freq_high_long_derived$allele_freq_high)))
allele_freq_gerp_long_derived$frequency <- as.numeric(substr(allele_freq_gerp_long_derived$allele_freq_gerp, 3, nchar(allele_freq_gerp_long_derived$allele_freq_gerp)))

## summary allele freq
summary(allele_freq_high_long_derived$frequency)
summary(allele_freq_gerp_long_derived$frequency)

## how many are fixed?
nrow(subset(allele_freq_high_long_derived, frequency == 1))
nrow(subset(allele_freq_gerp_long_derived, frequency == 1))

save(allele_freq_high_long_derived, file = "output/load/snpeff/allelefreq_snpeff_high.derived.RData")
save(allele_freq_gerp_long_derived, file = "output/load/gerp/allelefreq_gerp5.derived.RData")

## histograms
source("scripts/theme_ggplot.R")

ggplot(allele_freq_high_long_derived, aes(x = frequency)) + geom_histogram(fill = clrs[5], color=clrs[1]) + 
  labs(x = "Allele frequencies derived high impact mutations SnpEff", y = "Count") -> hist_af_high

ggplot(allele_freq_gerp_long_derived, aes(x = frequency)) + geom_histogram(fill = clrs[5], color=clrs[1])+ 
  labs(x = "Allele frequencies derived mutations with GERP score >= 4", y = "Count")-> hist_af_gerp

