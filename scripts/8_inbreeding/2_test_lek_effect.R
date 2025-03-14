### Test for lek effects on Froh ####
pacman::p_load(dplyr, data.table, brms)

load(file = "output/inbreeding/froh.RData")
load("data/phenotypes/phenotypes_lifetime.RData")

pheno_froh <- left_join(pheno_wide, froh, by = "id")

brm_lek_froh <- brm(scale(froh) ~ site, data = pheno_froh,
                    family = "gaussian",
                    prior = prior(normal(0,1), class = b),
                    cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                    iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

interval_lek_froh <- mcmc_intervals_data(brm_lek_froh, prob = 0.8, prob_outer = 0.95)
write.csv(interval_lek_froh, file = "output/load/interval_gerp_load_site_effect.csv", quote=F, row.names = F)
