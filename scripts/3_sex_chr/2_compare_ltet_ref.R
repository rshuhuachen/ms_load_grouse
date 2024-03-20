pacman::p_load(tidyverse, data.table)

files <- list.files(path = "output/genotyping/chicken_blast", 
                    pattern = ".out", full.names = T) 

list_out <- NULL
for (i in 1:length(files)){
    one_file <- fread(files[i])
    list_out[[i]] <- one_file}

for (i in 1:length(list_out)){
  names(list_out[[i]]) <- c("query", "subject", "perc_match", "length", "mismatch", "gap", "start_query", "end_query", "start_subject", "end_subject",
                            "evalue", "bitscore")
}

n_alignment <- lapply(list_out, nrow)
n_alignment <- as.data.frame(do.call(rbind, n_alignment))

mean_length <- lapply(list_out, function(x){mean(x$length)})
mean_length <- as.data.frame(do.call(rbind, mean_length))

max_length <- lapply(list_out, function(x){max(x$length)})
max_length <- as.data.frame(do.call(rbind, max_length))

sum_length <- lapply(list_out, function(x){sum(x$length)})
sum_length <- as.data.frame(do.call(rbind, sum_length))

alignment <- data.frame(name = files,
                        n_row = n_alignment$V1,
                        mean_length = mean_length$V1, 
                        max_length = max_length$V1,
                        sum_length = sum_length$V1)

alignment$name <- gsub("output/genotyping/chicken_blast/out_scaf_scaf_", "", alignment$name)
alignment$name <- gsub(".fa.out", "", alignment$name)
alignment$name <- gsub(";", ":", alignment$name)
alignment$name <- gsub("=", ".", alignment$name)

scafnames <- fread("data/genomic/refgenome/30_largest.scafs.tsv")

alignment <- left_join(alignment, scafnames[,c("scaf", "size_base", "n")], by = c("name" = "scaf"))

#### TBC #####
#also check with the alignment in the hal file to see the overlap
halstats <- fread("output/cactus/stats_reduced_363_hal.txt")
halstats <- subset(halstats, !grepl("bird|Anc|root", halstats$GenomeName))

