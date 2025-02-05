library(dplyr); library(data.table); library(readr)
args <- commandArgs(trailingOnly = TRUE)
genome <- args[[1]]
scafname <- args[[2]]
windows <- args[[3]]
winlength <- args[[4]]

#genome="data/genomes/30_largest.scafs_reduced.tsv"
gen <- fread(genome, col.names = c("scaf", "length"))
gen$scaf <- gsub(":", ";", gen$scaf)
gen$scaf <- gsub("\\.", "=", gen$scaf)

generate_windows <- function(scaf, length, window_length) {
  # Create a sequence of start coordinates
  length <- as.numeric(length)
  window_length <- as.numeric(window_length)
  starts <- seq(1, length, by = window_length)
  # Calculate the end coordinates
  ends <- starts + window_length - 1
  # Create a dataframe with start and end coordinates
  window_df <- data.frame(scaf = scaf, start = starts, end = ends)
  # Change last end to actual end
  window_df$end[nrow(window_df)] <- length
  return(window_df)
}

scafname<- gsub(":", ";", scafname)
scafname <- gsub("\\.", "=", scafname)

gen_scaf <- subset(gen, scaf == scafname)
windows_df <- generate_windows(scaf = scafname, 
                               length = gen_scaf$length, 
                               window_length = winlength)


write_tsv(windows_df, file = windows, col_names = F)


