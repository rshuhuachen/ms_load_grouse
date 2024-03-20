##### In this script we will separate the gff file into regions of interest based on the great tit ########

#### Packages #####
#remotes::install_github("jasongraf1/JGmisc")
pacman::p_load(BiocManager, rtracklayer, GenomicFeatures, BiocGenerics, data.table)

#### Genome data ####
# first change scaffold names
gff_raw <- fread("data/genomic/annotation/PO2979_Lyrurus_tetrix_black_grouse.annotation.gff") 
gff_raw$V1 <- gsub(";", "__", gff_raw$V1)
gff_raw$V1 <- gsub("=", "_", gff_raw$V1)
write.table(gff_raw, file = "data/genomic/PO2979_Lyrurus_tetrix_black_grouse.annotation_editedscafnames.gff", sep = "\t", col.names = FALSE, quote=F, row.names = FALSE)

#read in new file
gff <- makeTxDbFromGFF("data/genomic/PO2979_Lyrurus_tetrix_black_grouse.annotation_editedscafnames.gff", format="gff3", organism="Lyrurus tetrix") 

## divide up ###

promoters <- promoters(gff, upstream=2000, downstream=200, columns=c("tx_name", "gene_id")) # From NIOO
TSS <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) # TSS as in Laine et al., 2016. Nature Communications
downstream <- flank(genes(gff), 10000, start=FALSE, both=FALSE, use.names=TRUE)
upstream <- promoters(genes(gff), upstream=10000, downstream=0)

exons_gene <- unlist(exonsBy(gff, "gene")) # group exons by genes
introns <- unlist(intronsByTranscript(gff, use.names=TRUE))

### write out files

export(genes(gff), "data/genomic/annotation/genes.gff3")
export(promoters, "data/genomic/annotation/promoters.gff3")
export(TSS, "data/genomic/annotation/TSS.gff3")
export(downstream, "data/genomic/annotation/downstream.gff3")
export(upstream, "data/genomic/annotation/upstream.gff3")

export(exons_gene, "data/genomic/annotation/exons_gene.gff3")

introns@ranges@NAMES[is.na(introns@ranges@NAMES)]<-"Unknown"
export(introns, "data/genomic/annotation/introns_transcripts.gff3")

