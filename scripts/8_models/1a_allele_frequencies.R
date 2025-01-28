#### calculated allele frequencies of high impact mutations and GERP >= 4 ####
### here: load in, assign ancestral and derived alleles and calculate mean ####

pacman::p_load(tidyverse, data.table)

### First: reformat gerp to .vcf

load(file = "output/load/gerp/gerp_over4.RData")

# merge with clean ID names

gerps_vcf <- gerp %>% select(c(chr, pos, ancestral, derived, qual, info, format, D154259:D229777)) #exclude gerp scores
gerps_vcf <- gerps_vcf %>% mutate(ID = ".", .after = pos) %>% mutate(FILTER = ".", .after = qual)
names(gerps_vcf)[1:9] <- c("#CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT") #same names as .vcf

write_tsv(gerps_vcf, file = "output/load/gerp/gerps_over4.vcf", col_names = TRUE)

### calculate allele frequencies with vcftools
system("vcftools --gzvcf data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz --freq --out output/af/allele_freq_high")

system("vcftools --vcf output/load/gerp/gerps_over4.vcf --freq --out output/af/allele_freq_gerp45")

##
vcf_high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz")
vcf_gerp<-gerps_vcf
#vcf_gerp <- read.table("output/load/gerp/gerps_over4.vcf")

allele_freq_high <- fread("output/af/allele_freq_high.frq", skip = 1)
allele_freq_gerp <- fread("output/af/allele_freq_gerp45.frq", skip = 1)

names(allele_freq_high) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")
names(allele_freq_gerp) <- c("scaf",  "pos",  "n_alleles",  "n_chr",  "allele_freq_1", "allele_freq_2")

#### make long with one row one allele freq
allele_freq_high_long <- gather(allele_freq_high, allele, allele_freq_high, allele_freq_1:allele_freq_2, factor_key=TRUE)
allele_freq_gerp_long <- gather(allele_freq_gerp, allele, allele_freq_gerp, allele_freq_1:allele_freq_2, factor_key=TRUE)

### merge derived allele
allele_freq_high_long <- left_join(allele_freq_high_long, vcf_high[,c("V1", "V2", "V5")], by = c("scaf" = "V1", "pos" = "V2")) #chr, pos, derived
allele_freq_gerp_long <- left_join(allele_freq_gerp_long, vcf_gerp[,c("#CHROM", "POS", "ALT")], by = c("scaf" = "#CHROM", "pos" = "POS")) #chr, pos, derived

## subset where allele = derived
allele_freq_high_long_derived <- subset(allele_freq_high_long, substr(allele_freq_high,1,1) == V5)
allele_freq_gerp_long_derived <- subset(allele_freq_gerp_long, substr(allele_freq_gerp,1,1) == ALT)

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

