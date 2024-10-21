##### Bayesian models for LMS: total load #####

pacman::p_load(bayesplot, brms, dplyr, data.table)


args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
out <- args[[2]]
iter <- args[[3]]
warm <- args[[4]]
thin <- args[[5]]

method <- as.character(method)
# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
print(method)
load <- load_per_region %>% filter(loadtype == method)

# combine
pheno_wide_load <- left_join(pheno_wide, load, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load)) #some ids without genotypes, excluded for wgr

#### model ####
brm_load_t <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = iter, thin = thin, warmup = warm, seed = 1908)

save(brm_load_t, file = out)

