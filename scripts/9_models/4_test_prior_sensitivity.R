##### Testing prior sensitivity with 2 different sets of priors, lower iterations #####

### load packages ####
pacman::p_load(bayesplot, brms, dplyr, data.table)
set.seed(1908)

### load data ####
# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load_high <- load_per_region %>% filter(loadtype == "high")
load_gerp45 <- load_per_region %>% filter(loadtype == "gerp45")

# combine
load_high <- left_join(pheno_wide, load_high, by = "id")
load_high <- subset(load_high, !is.na(total_load)) #some ids without genotypes, excluded for wgr

load_gerp45 <- left_join(pheno_wide, load_gerp45, by = "id")
load_gerp45 <- subset(load_gerp45, !is.na(total_load)) #some ids without genotypes, excluded for wgr

#### model with different priors ####

### set pars
iter = 10000
thin = 100
warm = 5000

### GERP total ####
### prior 1: default
brm_load_gerp_p1 <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = load_gerp45,
                        family = "zero_inflated_poisson",
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_gerp_p1) # -0.21
intervals_gerp_p1 <- mcmc_intervals_data(brm_load_gerp_p1, prob =0.8, prob_outer = 0.95)
write.csv(intervals_gerp_p1, file = "output/models/prior_sensitivity/brm_total_load_gerp_p1.csv", quote=F, row.names = F)

### prior 2: positive
prior2 <- c(prior(normal(20,10), class = b),
            prior(normal(50,10), class = Intercept))

brm_load_gerp_p2 <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = load_gerp45,
                        family = "zero_inflated_poisson",
                        prior = prior2,
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_gerp_p2) # -0.20
intervals_gerp_p2 <- mcmc_intervals_data(brm_load_gerp_p2, prob =0.8, prob_outer = 0.95)
write.csv(intervals_gerp_p2, file = "output/models/prior_sensitivity/brm_total_load_gerp_p2.csv", quote=F, row.names = F)

### GERP hom het ####
### prior 1: default
brm_load_homhet_gerp_p1 <- brm(LMS_min ~ scale(hom_load) + scale(het_load) + core + (1|site), data = load_gerp45,
                        family = "zero_inflated_poisson",
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_homhet_gerp_p1) # -0.55 and -0.58 for hom and het
intervals_gerp_hom_het_p1 <- mcmc_intervals_data(brm_load_homhet_gerp_p1, prob =0.8, prob_outer = 0.95)
write.csv(intervals_gerp_hom_het_p1, file = "output/models/prior_sensitivity/brm_hom_het_load_gerp_p1.csv", quote=F, row.names = F)

### prior 2: positive
brm_load_homhet_gerp_p2 <- brm(LMS_min ~ scale(hom_load) + scale(het_load) + core + (1|site), data = load_gerp45,
                               family = "zero_inflated_poisson",
                               prior = prior2,
                               cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                               iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_homhet_gerp_p2) # -0.57 and -0.60 for hom and het
intervals_gerp_hom_het_p2 <- mcmc_intervals_data(brm_load_homhet_gerp_p2, prob =0.8, prob_outer = 0.95)
write.csv(intervals_gerp_hom_het_p2, file = "output/models/prior_sensitivity/brm_hom_het_load_gerp_p2.csv", quote=F, row.names = F)

#### High impact total ####
### prior 1: default
brm_load_high_p1 <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = load_high,
                  family = "zero_inflated_poisson",
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_high_p1) # -0.11
intervals_high_p1 <- mcmc_intervals_data(brm_load_high_p1, prob =0.8, prob_outer = 0.95)
write.csv(intervals_high_p1, file = "output/models/prior_sensitivity/brm_total_load_high_p1.csv", quote=F, row.names = F)

### prior 2: positive
brm_load_high_p2 <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = load_high,
                        family = "zero_inflated_poisson",
                        prior = prior2,
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_high_p2) # -0.20
intervals_high_p2 <- mcmc_intervals_data(brm_load_high_p2, prob =0.8, prob_outer = 0.95)
write.csv(intervals_high_p2, file = "output/models/prior_sensitivity/brm_total_load_high_p2.csv", quote=F, row.names = F)

### High hom het ####
### prior 1: default
brm_load_homhet_high_p1 <- brm(LMS_min ~ scale(hom_load) + scale(het_load) + core + (1|site), data = load_high,
                               family = "zero_inflated_poisson",
                               cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                               iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_homhet_high_p1) # -0.55 and -0.58 for hom and het
intervals_high_hom_het_p1 <- mcmc_intervals_data(brm_load_homhet_high_p1, prob =0.8, prob_outer = 0.95)
write.csv(intervals_high_hom_het_p1, file = "output/models/prior_sensitivity/brm_hom_het_load_high_p1.csv", quote=F, row.names = F)

### prior 2: positive
brm_load_homhet_high_p2 <- brm(LMS_min ~ scale(hom_load) + scale(het_load) + core + (1|site), data = load_high,
                               family = "zero_inflated_poisson",
                               prior = prior2,
                               cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                               iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_homhet_high_p2) # -0.21
intervals_high_hom_het_p2 <- mcmc_intervals_data(brm_load_homhet_high_p2, prob =0.8, prob_outer = 0.95)
write.csv(intervals_high_hom_het_p2, file = "output/models/prior_sensitivity/brm_hom_het_load_high_p2.csv", quote=F, row.names = F)

### test in one model ####
test <- left_join(load_gerp45[,c("id", "LMS_min", "core", "site", "total_load", "hom_load", "het_load")], 
                  load_high[,c("id", "total_load", "hom_load", "het_load")], by = "id", suffix = c("_gerp", "_high"))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load_gerp) + scale(total_load_high) + core + (1|site), family = "poisson", ziformula = ~1,
                         data = test))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load_gerp) + scale(total_load_high) + core + (1|site), family = "poisson", ziformula = ~1,
                         data = test))

ggplot(test, aes(x = hom_load_gerp, y = hom_load_high)) + geom_point() + geom_smooth(method='lm')
cor.test(test$hom_load_gerp, test$hom_load_high) # sig positive

ggplot(test, aes(x = het_load_gerp, y = het_load_high)) + geom_point() + geom_smooth(method='lm')
cor.test(test$het_load_gerp, test$het_load_high) # sig positive
