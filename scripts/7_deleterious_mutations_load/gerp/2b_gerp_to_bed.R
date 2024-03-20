library(dplyr); library(readr); library(stringr); library(data.table)
args <- commandArgs(trailingOnly = TRUE)
tsv_in <- args[[1]]
bed_out <- args[[2]]
path_out <- args[[3]]

scafnames <- fread("data/genomic/refgenome/30_largest.scafs.tsv")
scafnames <- scafnames[,c("scaf", "scaf_no")]
scafnames$scaf <- gsub(":", ";", scafnames$scaf)
scafnames$scaf <- gsub("\\.", "=", scafnames$scaf)

read_gerp <- function(file){
    scafnr <- gsub(".*maf_", "", file)
    scafnr <- gsub(".maf.rates", "", scafnr)
    read_tsv(file, col_names = c("neutral_rate_n", "rs_score")) %>%
    mutate(scaf_no = as.integer(scafnr),
           pos = row_number(),
           start = pos -1) -> out
    out <- right_join(scafnames, out, by = c("scaf_no"))
    out <- out %>% rename(seq=scaf)
    return(out)
}

dir.create(path = path_out,
           showWarnings = FALSE)

read_gerp(tsv_in) %>%
  select(seq, start, pos, neutral_rate_n, rs_score) %>%
  write_tsv(file = bed_out, col_names = FALSE)