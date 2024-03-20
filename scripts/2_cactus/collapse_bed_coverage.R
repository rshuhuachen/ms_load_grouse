library(dplyr); library(readr)

args <- commandArgs(trailingOnly = TRUE)
bed_in <- args[[1]]
bed_out <- args[[2]]

bed <- read_tsv(bed_in, col_names = c("seq", "start", "end", "id", "cov"))

bed %>%
  group_by(seq) %>%
  mutate(cov_bin = cumsum(!(cov == lag(cov, default = 0) & start == lag(end, default = 0) ))) %>%
  group_by(seq, cov_bin) %>%
  summarize(start = min(start),
            end = max(end),
            cov = cov[[1]]) %>%
  ungroup() %>%
  select(-cov_bin) %>%
  write_tsv( bed_out, col_names = FALSE )