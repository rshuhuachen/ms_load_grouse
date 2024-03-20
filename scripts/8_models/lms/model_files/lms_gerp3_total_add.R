##### Bayesian models for MS #####
pacman::p_load(bayesplot, brms, dplyr, data.table)

iter <- args[[1]]
warm <- args[[2]]
thin <- args[[3]]

#load data
file = "output/load/pheno_loads_lifetime.RData")

#### model ####

brm_load_t_add_gerp3_lms <- brm(LMS_min ~ scale(gerp_Lt_additive_cat3) + core + (1|site) + (1|born), data = pheno_load,
                               family = "zero_inflated_poisson",
                               prior = prior(normal(0,1), class = b),
                               cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                               iter = iter, thin = thin, warmup = warm)

save(brm_load_t_add_gerp3_lms, file = "results/bayes_models/lms_gerp3_total_add_zi_29scaf.RData")
