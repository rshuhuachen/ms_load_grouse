### annotate gene regions per SNP for high impact and gerp mutations ###

pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, tidyverse, GenomicRanges)

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

####### Gene regions ###########

### load annotation data
annotation_dir <- "data/genomic/annotation"

promoter=unique(gffToGRanges(paste0(annotation_dir, "/promoters.gff3")))
genes=unique(gffToGRanges(paste0(annotation_dir, "/genes.gff3")))
TSS=unique(gffToGRanges(paste0(annotation_dir, "/TSS.gff3")))
exons_gene=unique(gffToGRanges(paste0(annotation_dir, "/exons_gene.gff3")))
introns=unique(gffToGRanges(paste0(annotation_dir, "/introns_transcripts.gff3")))
downstream=unique(gffToGRanges(paste0(annotation_dir, "/downstream.gff3")))
upstream=unique(gffToGRanges(paste0(annotation_dir, "/upstream.gff3")))

#### Annotate SNPeff regions ####
snpef$end <- snpef$POS
snpef$start <- snpef$POS
snpef_gr <- as(snpef, "GRanges")

snpef_promoter <- as.data.frame(subsetByOverlaps(snpef_gr, promoter)) %>%
  add_column("region_promoter" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_genes <- as.data.frame(subsetByOverlaps(snpef_gr, genes)) %>%
  add_column("region_gene" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_TSS <- as.data.frame(subsetByOverlaps(snpef_gr, TSS))%>%
  add_column("region_tss" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_exons <- as.data.frame(subsetByOverlaps(snpef_gr, exons_gene))%>%
  add_column("region_exon" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_introns <- as.data.frame(subsetByOverlaps(snpef_gr, introns))%>%
  add_column("region_intron" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_down <- as.data.frame(subsetByOverlaps(snpef_gr, downstream))%>%
  add_column("region_down" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_up <- as.data.frame(subsetByOverlaps(snpef_gr, upstream))%>%
  add_column("region_up" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_all <- left_join(snpef, snpef_promoter[,c("chr", "start", "region_promoter")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_genes[,c("chr", "start", "region_gene")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_TSS[,c("chr", "start", "region_tss")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_exons[,c("chr", "start", "region_exon")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_introns[,c("chr", "start", "region_intron")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_down[,c("chr", "start", "region_down")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_up[,c("chr", "start", "region_up")], by = c("CHROM" = "chr", "start"))

snpef_all <- snpef_all %>% mutate(region_noncoding = as.factor(case_when(
  is.na(region_promoter) & is.na(region_gene) & is.na(region_tss) & is.na(region_exon) & 
    is.na(region_intron) & is.na(region_down) & is.na(region_up)  ~ 1)))

summary(snpef_all$region_noncoding)
save(snpef_all, file = "output/load/snpeff/snpeff_high_annotated_region.RData")

#### Annotate gerp regions ####
gerp$end <- gerp$POS
gerp$start <- gerp$POS
gerp_gr <- as(gerp, "GRanges")

gerp_promoter <- as.data.frame(subsetByOverlaps(gerp_gr, promoter)) %>%
  add_column("region_promoter" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_genes <- as.data.frame(subsetByOverlaps(gerp_gr, genes)) %>%
  add_column("region_gene" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_TSS <- as.data.frame(subsetByOverlaps(gerp_gr, TSS))%>%
  add_column("region_tss" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_exons <- as.data.frame(subsetByOverlaps(gerp_gr, exons_gene))%>%
  add_column("region_exon" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_introns <- as.data.frame(subsetByOverlaps(gerp_gr, introns))%>%
  add_column("region_intron" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_down <- as.data.frame(subsetByOverlaps(gerp_gr, downstream))%>%
  add_column("region_down" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_up <- as.data.frame(subsetByOverlaps(gerp_gr, upstream))%>%
  add_column("region_up" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_all <- left_join(gerp, gerp_promoter[,c("chr", "start", "region_promoter")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_genes[,c("chr", "start", "region_gene")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_TSS[,c("chr", "start", "region_tss")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_exons[,c("chr", "start", "region_exon")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_introns[,c("chr", "start", "region_intron")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_down[,c("chr", "start", "region_down")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_up[,c("chr", "start", "region_up")], by = c("CHROM" = "chr", "start")))

gerp_all <- gerp_all %>% mutate(region_noncoding = as.factor(case_when(
  is.na(region_promoter) & is.na(region_gene) & is.na(region_tss) & is.na(region_exon) & 
    is.na(region_intron) & is.na(region_down) & is.na(region_up)  1)))

summary(gerp_all$region_noncoding)

save(gerp_all, file = "output/load/gerp/gerp_annotated_region.RData")

#### Add in dataframe ###
snpef_snps <- nrow(snpef)
gerp_snps <- nrow(gerp)

data.frame("method" = c(rep("High impact SnpEff",8), rep("GERP scores >= 4",8)),
           "region" = c(rep(c("Downstream", "Exon",  "Gene body", "Intron", "Promoter",  "TSS", "Upstream", "Non-coding"),times=2)),
           "n_snps" = c(nrow(snpef_down), nrow(snpef_exons), nrow(snpef_genes), nrow(snpef_introns),
                        nrow(snpef_promoter),  nrow(snpef_TSS), nrow(snpef_up), nrow(subset(snpef_all,region_noncoding == 1)),
                        nrow(gerp_down), nrow(gerp_exons),  nrow(gerp_genes), nrow(gerp_introns),
                        nrow(gerp_promoter), nrow(gerp_TSS), nrow(gerp_up), nrow(subset(gerp_all,region_noncoding == 1))),
           "perc_snps" = c(nrow(snpef_down)/snpef_snps*100, nrow(snpef_exons)/snpef_snps*100, 
                           nrow(snpef_genes)/snpef_snps*100, nrow(snpef_introns)/snpef_snps*100,
                           nrow(snpef_promoter)/snpef_snps*100, 
                           nrow(snpef_TSS)/snpef_snps*100, nrow(snpef_up)/snpef_snps*100,
                           nrow(subset(snpef_all,region_noncoding == 1))/snpef_snps*100,
                           nrow(gerp_down)/gerp_snps*100, nrow(gerp_exons)/gerp_snps*100, 
                           nrow(gerp_genes)/gerp_snps*100, nrow(gerp_introns)/gerp_snps*100,
                           nrow(gerp_promoter)/gerp_snps*100, nrow(gerp_TSS)/gerp_snps*100, nrow(gerp_up)/gerp_snps*100,
                           nrow(subset(gerp_all,region_noncoding == 1))/gerp_snps*100)) -> snps_by_region

snps_by_region$region <- factor(snps_by_region$region, levels = c("Intron", 
                                                                  "Exon", "Upstream", "Downstream",
                                                                  "Gene body","TSS",
                                                                  "Promoter", "Non-coding"))

source("scripts/theme_ggplot.R")

save(snps_by_region, file = "output/snps_by_region_plot.RData")
load(file = "output/snps_by_region_plot.RData")

# number of snps, separate per method
ggplot(subset(snps_by_region, method == "GERP scores >= 4" & region != "Non-coding"), aes(y = region, x = n_snps)) + 
  geom_col(position="dodge", fill = clr_gerp) +
  labs(x = "Number of SNPs", y = "Region")+
  geom_text(aes(label = paste0("N = ", prettyNum(round(n_snps, 2), big.mark=",")), y = region, hjust=-0.3),
            position=position_dodge(width=0.9), size = 6) + 
  scale_x_continuous(labels = c("0", "50k", "100k", "150k", "200k"), 
                     breaks = c(0,50000,100000,150000, 200000), limits = c(0,250000))+
  theme(legend.position = "bottom") -> snps_per_region_plot_gerp

ggplot(subset(snps_by_region, method == "High impact SnpEff"& region != "Non-coding"), aes(y = region, x = n_snps)) + 
  geom_col(position="dodge", fill = clr_high) +
  labs(x = "Number of SNPs", y = "Region")+
  geom_text(aes(label = paste0("N = ", prettyNum(round(n_snps, 2), big.mark=",")), 
                y = region, hjust=-0.3),
            position=position_dodge(width=0.9), size = 6) + 
  theme(legend.position = "bottom") -> snps_per_region_plot_high

cowplot::plot_grid(snps_per_region_plot_gerp, snps_per_region_plot_high,
                   ncol = 2, 
                   align = "hv", axis = "lb", 
                   labels = "auto", label_fontface = "plain", label_size = 22)
