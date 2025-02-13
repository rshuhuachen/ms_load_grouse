#### Simplify LMS parameter here for modelling ####

# load packages
pacman::p_load(tidyverse, brms, bayesplot, data.table, ggridges)

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

### simplify LMS into binary operator ####
pheno_wide_load <- pheno_wide_load %>% mutate(LMS_binary = case_when(LMS_min > 0 ~ 1,
                                                                     LMS_min == 0 ~ 0))

#### model ####
brm_load_binary_high <- brm(LMS_binary ~ scale(total_load_high) + core + (1|site), data = pheno_wide_load,
                            family = "bernoulli",
                            prior = prior(normal(0,1), class = b),
                            cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                            iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_binary_high, file = "output/models/total_hom_het/lms_total_binary_high.RData")

brm_load_binary_gerp45 <- brm(LMS_binary ~ scale(total_load_gerp45) + core + (1|site), data = pheno_wide_load,
                              family = "bernoulli",
                              prior = prior(normal(0,1), class = b),
                              cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                            iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_binary_gerp45, file = "output/models/total_hom_het/lms_total_binary_gerp45.RData")


### simplify LMS into three categories ####
pheno_wide_load <- pheno_wide_load %>% mutate(LMS_cat = case_when(LMS_min > 5 ~ "High",
                                                                  LMS_min > 0 & LMS_min < 5 ~ "Moderate",
                                                                  LMS_min == 0 ~ "Low"))

pheno_wide_load$LMS_cat <- factor(pheno_wide_load$LMS_cat, levels = c("Low", "Moderate", "High"))

#### model ####
brm_load_cat_high <- brm(LMS_cat ~ scale(total_load_high) + core + (1|site), data = pheno_wide_load,
                            family = "categorical",
                            cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                            iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_cat_high, file = "output/models/total_hom_het/lms_total_cat_high.RData")

brm_load_cat_gerp45 <- brm(LMS_cat ~ scale(total_load_gerp45) + core + (1|site), data = pheno_wide_load,
                              family = "categorical",
                              cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                              iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_cat_gerp45, file = "output/models/total_hom_het/lms_total_cat_gerp45.RData")

