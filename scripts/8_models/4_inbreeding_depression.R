### packages ####
pacman::p_load(tidyverse, data.table, brms, bayesplot)

### load Froh values ####
load(file = "output/inbreeding/froh.RData")

### load pheno data ####
load("data/phenotypes/phenotypes_lifetime.RData") #LMS

### merge data ####
froh_pheno <- left_join(froh[,c("id", "froh")], pheno_wide[,c("id", "LMS_min", "core", "site", "born")], by = "id")

summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + core + (1|site) , data = froh_pheno, family = "poisson", ziformula = ~1))

### run brms model ####
iter = 1000000
burn = 500000
thin = 1000

#### model ####
fit <- brm(LMS_min ~ scale(froh) + core + (1|site), data = froh_pheno,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = iter, thin = thin, warmup = burn, seed = 1908)

save(fit, file = "output/models/from_lms.RData")