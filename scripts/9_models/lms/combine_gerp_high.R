### Packages ####
pacman::p_load(tidyverse, brms, bayesplot, data.table)

# function to calculate load
source("scripts/7_calculate_load/function_calculate_load.R")

### Approach 1: combining GERP and high mutations in one total load parameter ####
load(file = "output/load/snpeff/snpeff_high.RData")
load(file = "output/load/gerp/gerp_over4.RData")

names(gerp)[1:12]
names(snpeff)[1:12]

snpeff_c <- data.frame(chr = snpeff$CHROM,
                       start = snpeff$POS -1,
                       pos = snpeff$POS,
                       neutral_rate_n = NA,
                       rs_score = NA,
                       ancestral = snpeff$REF,
                       derived = snpeff$ALT,
                       qual = snpeff$QUAL,
                       info = snpeff$INFO,
                       format = snpeff$FORMAT)

snpeff_c <- cbind(snpeff_c, snpeff[c(10:199)])

combined <- rbind(gerp, snpeff_c)

load_combined <- calculate_load_gerp(vcf = combined, output_vcf = FALSE, loadtype= "gerp45_plus_high")

# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# combine
pheno_wide_load <- left_join(pheno_wide, load_combined, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load)) #some ids without genotypes, excluded for wgr

#### Model ####
summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                         family = "poisson", ziformula = ~1))


# brms model
brm_load_t <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_t, file = "output/models/total_hom_het/lms_total_gerp45_plus_high.RData")

brm_load_gerp_high <- brm_load_t
bayesplot::mcmc_intervals_data(brm_load_gerp_high, prob = 0.8, prob_outer = 0.95) %>%
  write.csv(file = "output/models/total_hom_het/lms_total_gerp45_plus_high_intervals.csv")

### Approach 2: both loads in the same model ####
# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load_high <- load_per_region %>% filter(loadtype == "high")
load_gerp45 <- load_per_region %>% filter(loadtype == "gerp45")

load <- left_join(load_high, load_gerp45, by = "id", suffix = c("_high", "_gerp45"))

# combine
pheno_wide_load <- left_join(pheno_wide, load, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load_high)) #some ids without genotypes, excluded for wgr

#### Model ####
brm_load_t <- brm(LMS_min ~ scale(total_load_high) + scale(total_load_gerp45) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_t, file = "output/models/total_hom_het/lms_total_gerp45_high_sep.RData")


bayesplot::mcmc_intervals_data(brm_load_gerp_high_sep, prob = 0.8, prob_outer = 0.95) %>%
  write.csv(file = "output/models/total_hom_het/lms_total_gerp45_high_sep_intervals.csv")

