#### Calculating sMLH and comparing sMLH with froh ####
#### Another method to estimate inbreeding is by calculating standardised multilocus heterozygosity
## Here we calculate sMLH to compare them to FROH

### Set up ####
### load packages
pacman::p_load(inbreedR, tidyverse,data.table,vcfR,reshape2)

### load in vcf file
VCF= paste0(getwd(), "data/genomic/intermediate/ltet_snps_filtered.vcf")

#### Calculate sMLH with inbreedR ####
### read vcf
vcf <- read.vcfR(VCF, verbose = FALSE )
# extract genotypes
gt <- extract.gt(vcf)
# transpose and transform to data.frame
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE)
# NA handling
gt[gt == "."] <- NA
# split columns
snp_geno <- do.call(cbind, apply(gt, 2, function(x) colsplit(x, "/", c("a","b"))))
# convert to raw inbreedR format
snp_genotypes <- inbreedR::convert_raw(snp_geno)
# check data
check_data(snp_genotypes)

# calculate smlh
het <- sMLH(snp_genotypes)
smlh <- data.frame(smlh = het)

### Add animal ID ####
# get id from vcf
# system("/vol/biotools/bin/bcftools query -l $VCF > data/genomic/raw/metadata/ids_vcf_ltet.txt")

file_id <- read.table("/prj/blackgrouse/githubs/rohcalling_newest/data/ids_vcf_ltet.txt")
smlh$file_id <- file_id[,1]
smlh$file_id <- gsub("/vol/cluster-data/rchen/wgr/data/processed/", "", smlh$file_id)

#add read id
idfile <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
smlh <- left_join(smlh, unique(idfile[,c(3,9)]), by = c("file_id" = "file"))
smlh <- smlh[,-2]

write.csv(smlh, "output/inbreeding/sMLH.csv", row.names=F, quote=F)

