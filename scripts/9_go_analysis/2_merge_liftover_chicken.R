
###### Based on function: chicken annotation #####
pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, tidyverse, GenomicRanges)

### load annotation data

chicken_ann_txt <- unique(fread("data/genomic/annotation/gallus_gallus/liftoff_gallus_ltet.gff"))
names(chicken_ann_txt) <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "info")

chicken_ann_txt$seqnames <- gsub(";", "__", chicken_ann_txt$seqnames)
chicken_ann_txt$seqnames <- gsub("=", "_", chicken_ann_txt$seqnames)

# extract just the 'gene' part in the info column for later GO analyses
chicken_ann_txt$gene <- gsub(".*gene=([^;]+).*", "\\1", chicken_ann_txt$info)

chicken_ann_gr <- as(chicken_ann_txt, "GRanges")

## annotate:snpeff, gerp

### load data ###

snpef <- read.table("output/snpeff/annotation/ltet_ann_aa_snp_output_HIGH.vcf")
gerp <- read.table("output/gerp/gerps_cat5_vcfformat.vcf")

#### rename columns #####
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(snpef) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns
names(gerp) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns

### change scaf names
snpef$CHROM <- gsub(";", "__", snpef$CHROM)
snpef$CHROM <- gsub("=", "_", snpef$CHROM)

gerp$CHROM <- gsub(";", "__", gerp$CHROM)
gerp$CHROM <- gsub("=", "_", gerp$CHROM)

### remove genotypes: not necessary here
snpef <- snpef[,c(1:9)]
gerp <- gerp[,c(1:9)]

snpef$end <- snpef$POS
snpef$start <- snpef$POS
snpef_gr <- as(snpef, "GRanges")

gerp$end <- gerp$POS
gerp$start <- gerp$POS
gerp_gr <- as(gerp, "GRanges")

#merge by overlap

snpef_annotated_chicken <- as.data.frame(mergeByOverlaps(snpef_gr, chicken_ann_gr, ignore.strand=TRUE))

gerp_annotated_chicken <- as.data.frame(mergeByOverlaps(gerp_gr, chicken_ann_gr))

save(gerp_annotated_chicken, file = "output/load/gerp/chicken_annotated_gerp.RData")
save(snpef_annotated_chicken, file = "output/load/snpeff/chicken_annotated_high.RData")
