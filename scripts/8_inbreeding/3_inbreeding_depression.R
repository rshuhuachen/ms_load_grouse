### Here, we run the formal brms model to test for inbreeding depression ####
# Could be done on a server, but it's just one model. If you don't have many CPU's, this can be intensive to run within R(Studio)!

### packages ####
pacman::p_load(tidyverse, data.table, brms, bayesplot)

### load Froh values ####
load(file = "output/inbreeding/froh.RData")

### load pheno data ####
load("data/phenotypes/phenotypes_lifetime.RData") #LMS

### merge data ####
froh_pheno <- left_join(froh[,c("id", "froh")], pheno_wide[,c("id", "LMS_min", "core", "site", "born")], by = "id")

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

#### Model diagnosis ####
#get posteriors
posterior <- as.array(fit)
log_ps <- log_posterior(fit)
nuts <- nuts_params(fit) #divergence

#get only beta and sd
betas <- variables(fit)[grep("b_", variables(fit))]
sd <- variables(fit)[grep("sd_", variables(fit))]

#global patterns in divergence
diverge_beta <- mcmc_parcoord(posterior, np = nuts, pars= betas)

#identify collinearity between parameters
collin_beta <- mcmc_pairs(posterior, np = nuts, pars= betas)

#traceplot
trace_beta <- mcmc_trace(posterior, pars = betas, np = nuts)

#rhat
rhat <- mcmc_rhat(brms::rhat(fit))

#effective sample size
neff <- mcmc_neff(neff_ratio(fit))

#autocorrelation
autocor_beta <- mcmc_acf(posterior, pars = betas)
autocor_sd <- mcmc_acf(posterior, pars=sd)

#quick glance results
areas <- mcmc_areas(fit, pars=betas)

#combine in list
diagnosis <- list(diverge_beta = diverge_beta, 
                collin_beta = collin_beta, 
                trace_beta = trace_beta, 
                rhat = rhat, 
                neff = neff, 
                autocor_beta = autocor_beta, 
                autocor_sd = autocor_sd,
                areas = areas)



# save output
pdf(file=paste0("output/models/diagnosis/froh_lms.pdf"))
print(diverge_beta)
print(collin_beta)
print(trace_beta)
print(rhat)
print(neff)
print(autocor_beta)
print(autocor_sd)
print(areas)
dev.off()

#### Testing for prior sensitivity #####
### prior 1: default
brm_froh_p1 <- brm(LMS_min ~ scale(froh) + core + (1|site), data = froh_pheno,
                        family = "zero_inflated_poisson",
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_load_gerp_p1) # similar

### prior 2: positive
prior2 <- c(prior(normal(20,10), class = b),
            prior(normal(50,10), class = Intercept))

brm_froh_p2 <- brm(LMS_min ~ scale(froh) + core + (1|site), data = froh_pheno,
                        family = "zero_inflated_poisson",
                        prior = prior2,
                        cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                        iter = iter, thin = thin, warmup = warm, seed = 1908)

summary(brm_froh_p2) # similar
