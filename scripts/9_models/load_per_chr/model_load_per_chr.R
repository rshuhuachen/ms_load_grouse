##### Bayesian models for LMS: total load #####

pacman::p_load(bayesplot, brms, dplyr, data.table)

args <- commandArgs(trailingOnly = TRUE)
method_arg <- args[[1]]
chr <- args[[2]]
out <- args[[3]]
iter <- args[[4]]
warm <- args[[5]]
thin <- args[[6]]

method_arg <- as.character(method_arg)
chr <- as.character(chr)

# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# load mutation load measures
load(file = "output/load/load_per_chr_snpeff_gerp.RData")

# subset only the relevant method/loadtype
print(method_arg)
print(chr)
load <- per_chr_load %>% filter(method == method_arg & loadtype == chr)

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

