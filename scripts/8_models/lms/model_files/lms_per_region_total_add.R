pacman::p_load(bayesplot, brms, dplyr, data.table, readr)
args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
region <- args[[2]]
out <- args[[3]]

##### Bayesian models for MS #####

#load data
file = "output/load/pheno_loads_lifetime.RData")

#### brms settings ####
iter = 1000000
thin = 1000
warm = 1/2*iter

#### Model
parameter <- paste0(gsub("_region", "", method), "_load_t_add_", gsub("region_", "", region))
formula <- formula(paste0("LMS_min ~ scale(", parameter, ") + core + (1|site) + (1|born)"))

brm_region <- brm(formula, data = pheno_load,
                               family = "zero_inflated_poisson",
                               prior = prior(normal(0,1), class = b),
                               cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                               iter = iter, thin = thin, warmup = warm)

save(brm_region, file = out)

