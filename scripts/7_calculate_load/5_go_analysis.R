### based on chicken gene annotation, what genes do we find enriched? ###

pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, dplyr, GenomicRanges)

### load data ###

load("output/load/snpeff/snpeff_high.RData") # snpeff mutations
load("output/load/gerp/gerp_over4.RData") # gerp mutations

all <- read.table("output/ancestral/ltet_filtered_ann_aa.vcf")

#### rename columns #####
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

gerp <- gerp %>% dplyr::select(c(chr, pos, ancestral, derived, qual, info, format, D154259:D229777))
gerp <- gerp %>% mutate(ID = NA, .after =  pos) %>% mutate(FILTER = NA, .after = qual)
names(gerp) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns
names(all) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns

### change scaf names

all$CHROM <- gsub(";", "__", all$CHROM)
all$CHROM <- gsub("=", "_", all$CHROM)

snpeff$CHROM <- gsub(";", "__", snpeff$CHROM)
snpeff$CHROM <- gsub("=", "_", snpeff$CHROM)

gerp$CHROM <- gsub(";", "__", gerp$CHROM)
gerp$CHROM <- gsub("=", "_", gerp$CHROM)

### remove genotypes: not necessary here
snpef <- snpeff[,c(1:9)]
gerp <- gerp[,c(1:9)]
all <- all[,c(1:9)]

####### Gene regions ###########

### load annotation data
annotation_dir <- "/vol/cluster-data/rchen/git/grouse-annotation/output"
#annotation_dir <- "/Users/vistor/Documents/Work/GitHub/PhD/grouse-annotation/output"

promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
genes=unique(gffToGRanges(paste0(annotation_dir, "/genes.gff3")))
TSS=unique(gffToGRanges(paste0(annotation_dir, "/TSS.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))
downstream=unique(gffToGRanges(paste0(annotation_dir, "/downstream.gff3")))
upstream=unique(gffToGRanges(paste0(annotation_dir, "/upstream.gff3")))
threeUTR =unique(gffToGRanges(paste0(annotation_dir, "/threeUTRs.gff3")))
fiveUTR=unique(gffToGRanges(paste0(annotation_dir, "/fiveUTRs.gff3")))

#### Annotate SNPeff regions ####
snpeff$end <- snpeff$POS
snpeff$start <- snpeff$POS
snpef_gr <- as(snpeff, "GRanges")

sig_promoter_snpef <- mergeByOverlaps(promoter, snpef_gr)
sig_promoter_snpef <- as.data.frame(sig_promoter_snpef@listData)
sig_promoter_snpef <- sig_promoter_snpef %>% add_column("region" = "promoter") %>% 
  dplyr::select(c(`snpef_gr.seqnames`, POS, region, gene_id)) 

sig_gene_snpef <- mergeByOverlaps(snpef_gr, genes)
sig_gene_snpef <- as.data.frame(sig_gene_snpef@listData)
sig_gene_snpef <- sig_gene_snpef %>% add_column("region" = "gene") %>% 
  dplyr::select(c(`snpef_gr.seqnames`, POS, region, gene_id)) 

sig_tss_snpef <- mergeByOverlaps(snpef_gr, TSS)
sig_tss_snpef <- as.data.frame(sig_tss_snpef@listData)
sig_tss_snpef <- sig_tss_snpef %>% add_column("region" = "tss") %>% 
  dplyr::select(c(`snpef_gr.seqnames`, POS, region, gene_id)) 

sig_exon_snpef <- mergeByOverlaps(snpef_gr, exons_gene)
sig_exon_snpef <- as.data.frame(sig_exon_snpef@listData)
sig_exon_snpef <- sig_exon_snpef %>% add_column("region" = "exon") %>% 
  dplyr::select(c(`snpef_gr.seqnames`, POS, region, exons_gene.ID)) 
sig_exon_snpef <- sig_exon_snpef %>% rename("gene_id" = "exons_gene.ID")

sig_intron_snpef <- mergeByOverlaps(snpef_gr, introns)
sig_intron_snpef <- as.data.frame(sig_intron_snpef@listData)
sig_intron_snpef <- sig_intron_snpef %>% add_column("region" = "intron") %>% 
  dplyr::select(c(`snpef_gr.seqnames`, POS, region, introns.ID)) 
sig_intron_snpef <- sig_intron_snpef %>% rename("gene_id" = "introns.ID")

snpef_all <- rbind(sig_promoter_snpef, sig_gene_snpef, sig_tss_snpef, sig_exon_snpef, sig_intron_snpef)

save(snpef_all, file = "output/load/snpeff/snpeff_high_annotated_region_gene_ids.RData")

#### Annotate gerp regions ####
gerp$end <- gerp$POS
gerp$start <- gerp$POS
gerp_gr <- as(gerp, "GRanges")

sig_promoter_gerp <- mergeByOverlaps(promoter, gerp_gr)
sig_promoter_gerp <- as.data.frame(sig_promoter_gerp@listData)
sig_promoter_gerp <- sig_promoter_gerp %>% add_column("region" = "promoter") %>% 
  dplyr::select(c(`gerp_gr.seqnames`, POS, region, gene_id)) 

sig_gene_gerp <- mergeByOverlaps(gerp_gr, genes)
sig_gene_gerp <- as.data.frame(sig_gene_gerp@listData)
sig_gene_gerp <- sig_gene_gerp %>% add_column("region" = "gene") %>% 
  dplyr::select(c(`gerp_gr.seqnames`, POS, region, gene_id)) 

sig_tss_gerp <- mergeByOverlaps(gerp_gr, TSS)
sig_tss_gerp <- as.data.frame(sig_tss_gerp@listData)
sig_tss_gerp <- sig_tss_gerp %>% add_column("region" = "tss") %>% 
  dplyr::select(c(`gerp_gr.seqnames`, POS, region, gene_id)) 

sig_exon_gerp <- mergeByOverlaps(gerp_gr, exons_gene)
sig_exon_gerp <- as.data.frame(sig_exon_gerp@listData)
sig_exon_gerp <- sig_exon_gerp %>% add_column("region" = "exon") %>% 
  dplyr::select(c(`gerp_gr.seqnames`, POS, region, exons_gene.ID)) 
sig_exon_gerp <- sig_exon_gerp %>% rename("gene_id"= "exons_gene.ID")

sig_intron_gerp <- mergeByOverlaps(gerp_gr, introns)
sig_intron_gerp <- as.data.frame(sig_intron_gerp@listData)
sig_intron_gerp <- sig_intron_gerp %>% add_column("region" = "intron") %>% 
  dplyr::select(c(`gerp_gr.seqnames`, POS, region, introns.ID)) 
sig_intron_gerp <- sig_intron_gerp %>% rename("gene_id" = "introns.ID")

gerp_all <- rbind(sig_promoter_gerp, sig_gene_gerp, sig_tss_gerp, sig_exon_gerp, sig_intron_gerp)

save(gerp_all, file = "output/load/gerp/gerp_high_annotated_region_gene_ids.RData")

### Annotate ALL mutations ####
all$end <- all$POS
all$start <- all$POS
all_gr <- as(all, "GRanges")

sig_promoter_all <- mergeByOverlaps(promoter, all_gr)
sig_promoter_all <- as.data.frame(sig_promoter_all@listData)
sig_promoter_all <- sig_promoter_all %>% add_column("region" = "promoter") %>% 
  dplyr::select(c(`all_gr.seqnames`, POS, region, gene_id)) 

sig_gene_all <- mergeByOverlaps(all_gr, genes)
sig_gene_all <- as.data.frame(sig_gene_all@listData)
sig_gene_all <- sig_gene_all %>% add_column("region" = "gene") %>% 
  dplyr::select(c(`all_gr.seqnames`, POS, region, gene_id)) 

sig_tss_all <- mergeByOverlaps(all_gr, TSS)
sig_tss_all <- as.data.frame(sig_tss_all@listData)
sig_tss_all <- sig_tss_all %>% add_column("region" = "tss") %>% 
  dplyr::select(c(`all_gr.seqnames`, POS, region, gene_id)) 

sig_exon_all <- mergeByOverlaps(all_gr, exons_gene)
sig_exon_all <- as.data.frame(sig_exon_all@listData)
sig_exon_all <- sig_exon_all %>% add_column("region" = "exon") %>% 
  dplyr::select(c(`all_gr.seqnames`, POS, region, exons_gene.ID)) 
sig_exon_all <- sig_exon_all %>% rename("gene_id" = "exons_gene.ID")

sig_intron_all <- mergeByOverlaps(all_gr, introns)
sig_intron_all <- as.data.frame(sig_intron_all@listData)
sig_intron_all <- sig_intron_all %>% add_column("region" = "intron") %>% 
  dplyr::select(c(`all_gr.seqnames`, POS, region, introns.ID)) 
sig_intron_all <- sig_intron_all %>% rename("gene_id" = "introns.ID")

all_all <- rbind(sig_promoter_all, sig_gene_all, sig_tss_all, sig_exon_all, sig_intron_all)

save(all_all, file = "output/load/all_snps_annotated_region_gene_ids.RData")

#### merge with 'similar to' column to get gene IDs ####
lookup <- fread("data/genomic/refgenome/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

snpef_all <- left_join(snpef_all, lookup, by = c("gene_id" = "original"))
gerp_all <- left_join(gerp_all, lookup, by = c("gene_id" = "original"))
all_all <- left_join(all_all, lookup, by = c("gene_id" = "original"))

snpef_all$similar <- toupper(snpef_all$similar)
gerp_all$similar <- toupper(gerp_all$similar)
all_all$similar <- toupper(all_all$similar)

library(readr)
snpef_all %>% dplyr::select(similar) %>% unique() %>% write_tsv(file = "output/go_analysis/gene_id_snpef_high.tsv")
gerp_all %>% dplyr::select(similar) %>% unique() %>% write_tsv(file = "output/go_analysis/gene_id_gerp.tsv")
all_all %>% dplyr::select(similar) %>% unique() %>% write_tsv(file = "output/go_analysis/gene_id_all.tsv")


