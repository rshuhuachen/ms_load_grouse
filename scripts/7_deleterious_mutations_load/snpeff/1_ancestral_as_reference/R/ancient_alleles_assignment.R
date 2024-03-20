library(tidyverse)
library(here)
library(dplyr)
setwd("/vol/cluster-data/rchen/git/genetic_load_ltet/")
# import sites diverging in ancestral genome and reference genome
data_anc <- read_tsv(here("scripts/snpeff/1_ancestral_as_reference/processed_data/ancestral_66_ltet_snps.tsv.gz"), skip = 1, col_names = F) 
names(data_anc) <- c("refSequence", "refPosition", "Lyrurustetrix", "Anc66")
data_anc <- data_anc %>% mutate(POS  = refPosition + 1)

#rename the scaffolds, for some reason cactus renamed them within the hal which was used to get the snps, now convert back to the actual names
data_anc$refSequence <- gsub(":", ";", data_anc$refSequence)
data_anc$refSequence <- gsub("\\.", "=", data_anc$refSequence)

data_anc <- unique(data_anc) #some duplicates weirdly

# import variable sites in population (SNPs)
data_snps <- read_tsv(here("scripts/snpeff/1_ancestral_as_reference/processed_data/ltet_snps_vcf.tsv.gz"))
data_snps <- unique(data_snps)

# combine SNPs and ancestral data
data_joined <- left_join(data_snps,
                         data_anc,
                         by = c(`#CHROM` = "refSequence", POS = "POS")) %>% 
  filter(complete.cases(refPosition))   # cases where ltet and anc are identical don't need updating

# allele complement assignment (and switch to upper case)
compl <- c(a = "T", A = "T", c = "G", C = "G", g = "C", G = "C", t = "A", T = "A")

# check if the original or its complementary allele is needed, based on REF:
# -> adjust Anc accordingly, THEN check if Anc == ALT or Anc %in% ALT...
data_check_anc <- data_joined %>% 
  mutate(across(c(REF, ALT, Lyrurustetrix, Anc66), str_to_upper)) %>% 
  mutate(ref_correct = case_when(
    REF == Lyrurustetrix ~ TRUE,           # ltet IS REF (no adjustment needed)
    REF == compl[Lyrurustetrix] ~ FALSE,   # ltet is complement of REF (need complement of Anc)
    .default = NA ),                # sanity check (ltet is neither REF nor its complement, should not occur -> check below)
  anc_corrected = if_else(ref_correct, Anc66, compl[Anc66]),
  anc_in_alt = map2_lgl(anc_corrected, ALT, grepl))# %>% 

# check how often REF is neither matching nor complement (should be 0)
sum(is.na(data_check_anc$ref_correct))
# check how often Anc in not listed in Alt
n_alt_non_match <- sum(1 - data_check_anc$anc_in_alt)

# fraction of diverging ancestral alleles not within the ltet population (relative to all diverging alleles)
f_miss_div <- n_alt_non_match / nrow(data_check_anc)
# fraction of diverging ancestral alleles not within the ltet population (relative to all variable sites)
f_miss_total <- n_alt_non_match / nrow(data_snps)
# fraction of variable sites where REF != Anc 
f_div_total <- nrow(data_check_anc) / nrow(data_snps)

tibble( type = c("div Anc / total SNPs", "miss Anc / div Anc", "miss Anc / total SNPs"),
        `%` = sprintf("%.2f", c(f_div_total, f_miss_div, f_miss_total) * 100),
        `n (10^6)` = str_c(sprintf("%.2f", c(nrow(data_check_anc), n_alt_non_match, n_alt_non_match)* 1e-6), 
                  "/",
                  sprintf("%.2f", c( nrow(data_snps), nrow(data_check_anc), nrow(data_snps)) * 1e-6))) %>% 
  knitr::kable(format = "latex") %>% 
  write_lines(here("scripts/snpeff/1_ancestral_as_reference/results/ancestral/ancestral_allele_mismatches.tex"))

# assign new ref (available ancestral allele) for vcf updating
data_new_ref <- data_snps %>% 
  left_join(data_check_anc %>% 
              # remove diverging alleles that are not in the population
              filter(anc_in_alt) %>%  
              select(`#CHROM`, POS, AA = anc_corrected),
            by = c("#CHROM", "POS")) %>% 
  # if allele is not diverging in Anc, keep REF
  mutate(AA = if_else(is.na(AA), REF, AA))

data_new_ref %>% 
  write_tsv(here("scripts/snpeff/1_ancestral_as_reference/results/ancestral/anc_allele_assignment.tsv.gz"))

data_new_ref %>% 
  select(`#CHROM`, FROM = POS, TO = POS, AA) %>%
  write_tsv(here("scripts/snpeff/1_ancestral_as_reference/results/ancestral/anc_allele_assignment.bed"))