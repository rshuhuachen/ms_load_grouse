### based on chicken gene annotation, what genes do we find enriched? ###

pacman::p_load(genomation, data.table, GenomicFeatures, rtracklayer, fuzzyjoin, tibble, dplyr, GenomicRanges)

### load data ###

load("output/load/snpeff/snpeff_high.RData") # snpeff mutations
load("output/load/gerp/gerp_over4.RData") # gerp mutations

all <- read.table("output/ancestral/ltet_filtered_ann_aa.vcf")

#### rename columns #####
filenames <- fread("data/metadata/ltet_snps_noindel_minDP20_maxmis0.7_maxdp60_minQC30_pruned_hwe.fam")
ids <- fread("data/metadata/file_list_all_bgi_clean.csv")
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(snpef) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns
names(gerp) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns
names(all) <- c(c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT"), idnames$id) #rename columns

### change scaf names
snpef$CHROM <- gsub(";", "__", snpef$CHROM)
snpef$CHROM <- gsub("=", "_", snpef$CHROM)

gerp$CHROM <- gsub(";", "__", gerp$CHROM)
gerp$CHROM <- gsub("=", "_", gerp$CHROM)

### remove genotypes: not necessary here
snpef <- snpef[,c(1:9)]
gerp <- gerp[,c(1:9)]
all <- all[,c(1:9)]

##### Two types of annotation: regions vs functions #####

# Since we don't have gene names just gene IDs for grouse, use functional based one on chicken instead
# So first: annotate gene regions, then annotate with genes for GO 

####### Gene regions ###########
###### based on grouse

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

snpef_threeUTR <- as.data.frame(subsetByOverlaps(snpef_gr, threeUTR))%>%
  add_column("region_threeutr" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

snpef_fiveUTR <- as.data.frame(subsetByOverlaps(snpef_gr, fiveUTR))%>%
  add_column("region_fiveutr" = 1) %>% dplyr::rename(chr = seqnames) %>% unique()

snpef_all <- left_join(snpef, snpef_promoter[,c("chr", "start", "region_promoter")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_genes[,c("chr", "start", "region_gene")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_TSS[,c("chr", "start", "region_tss")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_exons[,c("chr", "start", "region_exon")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_introns[,c("chr", "start", "region_intron")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_down[,c("chr", "start", "region_down")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_up[,c("chr", "start", "region_up")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_threeUTR[,c("chr", "start", "region_threeutr")], by = c("CHROM" = "chr", "start"))%>%
  left_join(snpef_fiveUTR[,c("chr", "start", "region_fiveutr")],by = c("CHROM" = "chr", "start"))

snpef_all <- snpef_all %>% mutate(region_noncoding = as.factor(case_when(
  is.na(region_promoter) & is.na(region_gene) & is.na(region_tss) & is.na(region_exon) & 
    is.na(region_intron) & is.na(region_down) & is.na(region_up) & 
    is.na(region_threeutr) & is.na(region_fiveutr) ~ 1)))

summary(snpef_all$region_noncoding)
save(snpef_all, file = "results/snpeff/snpeff_high_annotated_region.RData")

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

gerp_threeUTR <- as.data.frame(subsetByOverlaps(gerp_gr, threeUTR))%>%
  add_column("region_threeutr" = 1) %>% dplyr::rename(chr = seqnames)%>% unique()

gerp_fiveUTR <- as.data.frame(subsetByOverlaps(gerp_gr, fiveUTR))%>%
  add_column("region_fiveutr" = 1) %>% dplyr::rename(chr = seqnames) %>% unique()

gerp_all <- left_join(gerp, gerp_promoter[,c("chr", "start", "region_promoter")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_genes[,c("chr", "start", "region_gene")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_TSS[,c("chr", "start", "region_tss")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_exons[,c("chr", "start", "region_exon")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_introns[,c("chr", "start", "region_intron")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_down[,c("chr", "start", "region_down")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_up[,c("chr", "start", "region_up")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_threeUTR[,c("chr", "start", "region_threeutr")], by = c("CHROM" = "chr", "start"))%>%
  left_join(gerp_fiveUTR[,c("chr", "start", "region_fiveutr")],by = c("CHROM" = "chr", "start"))

gerp_all <- gerp_all %>% mutate(region_noncoding = as.factor(case_when(
  is.na(region_promoter) & is.na(region_gene) & is.na(region_tss) & is.na(region_exon) & 
    is.na(region_intron) & is.na(region_down) & is.na(region_up) & 
    is.na(region_threeutr) & is.na(region_fiveutr) ~ 1)))

summary(gerp_all$region_noncoding)

save(gerp_all, file = "results/gerp/gerp_annotated_region.RData")

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

save(snps_by_region, file = "results/snps_by_region_plot.RData")
load(file = "results/snps_by_region_plot.RData")

ggplot(snps_by_region, aes(y = region, x = perc_snps, fill = method)) + 
  geom_col(position="dodge") +
  scale_fill_manual(values = c(clr_gerp, clr_high)) + 
  labs(x = "Percentage of SNPs", y = "Genomic region", fill = "Approach used")+
  geom_text(aes(label = paste0(prettyNum(round(perc_snps, 2), big.mark=","), "%"), group = method, y = region, hjust=-0.3),
                position=position_dodge(width=0.9)) + 
  theme(legend.position = "bottom") +
  xlim(0, 110) -> snps_per_region_plot

ggsave(snps_per_region_plot, file = "plots/final/snps_deleterious_per_gene_region.png", width=10, height=10)

# number instead, separate per method
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

###### Based on function: chicken annotation #####

### load annotation data
#chicken_ann_gr <- unique(gffToGRanges("results/go_analysis/liftoff_gallus_ltet.gff")) -> not working like this, first convert to correct format in txt
#do this manually....
chicken_ann_txt <- unique(fread("results/go_analysis/liftoff_gallus_ltet.gff"))
names(chicken_ann_txt) <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "info")

chicken_ann_txt$seqnames <- gsub(";", "__", chicken_ann_txt$seqnames)
chicken_ann_txt$seqnames <- gsub("=", "_", chicken_ann_txt$seqnames)

# extract just the 'gene' part in the info column for later GO analyses
chicken_ann_txt$gene <- gsub(".*gene=([^;]+).*", "\\1", chicken_ann_txt$info)

chicken_ann_gr <- as(chicken_ann_txt, "GRanges")

## annotate: entire vcf, snpeff, gerp
all$end <- all$POS
all$start <- all$POS

all$CHROM <- gsub(";", "__", all$CHROM)
all$CHROM <- gsub("=", "_", all$CHROM)

all_gr <- as(all, "GRanges")

#merge by overlap

snpef_annotated_chicken <- as.data.frame(mergeByOverlaps(snpef_gr, chicken_ann_gr, ignore.strand=TRUE))

gerp_annotated_chicken <- as.data.frame(mergeByOverlaps(gerp_gr, chicken_ann_gr))

all_annotated_chicken <- as.data.frame(mergeByOverlaps(all_gr, chicken_ann_gr))

save(gerp_annotated_chicken, file = "results/go_analysis/chicken_annotated_gerp.RData")
save(snpef_annotated_chicken, file = "results/go_analysis/chicken_annotated_high.RData")

### export gene names

all_genes <- data.frame(genes = unique(all_annotated_chicken$chicken_ann_gr.gene))
snpef_genes <- data.frame(genes = unique(snpef_annotated_chicken$chicken_ann_gr.gene))
gerp_genes <- data.frame(genes = unique(gerp_annotated_chicken$chicken_ann_gr.gene))

write.csv(all_genes, file = "results/go_analysis/genes_all.csv", row.names = F, quote=F)
write.csv(snpef_genes, file = "results/go_analysis/genes_snpef_high.csv", row.names = F, quote=F)
write.csv(gerp_genes, file = "results/go_analysis/genes_gerp_cat5.csv", row.names = F, quote=F)


#### merge with 'similar to' column to get gene IDs ####
lookup <- fread("../grouse-annotation/data/lookup_ANN_gene_id.txt")
names(lookup) <- c("original", "similar")

gerp_all <- left_join(gerp_all, lookup, by = c("ID" = "original"))

names(gerp_all) <- c("chr", "pos", "region", "ann", "gene_id")

gerp_all$gene_id <- toupper(gerp_all$gene_id)
