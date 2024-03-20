library(dplyr); library(data.table); library(readr); library(brms)

args <- commandArgs(trailingOnly = TRUE)
response <- args[[1]]
iterations <- args[[2]]
burn <- args[[3]]
thin <- args[[4]]

### load data ###

load(file = "output/load/pheno_loads_annual.RData")
pheno_core <- subset(pheno_annual_load, id %in% core$id)

#annual measurements for core males only

### modelling ####
formula <- formula(paste0("scale(", response, ") ~ scale(froh_bcftools) + age + sqrt(age) + (1|year) + (1|site/id)"))

fit <- brm(formula,
             family = "gaussian",
           data = pheno_core, 
           cores =8,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,1), class = b),
           iter = iterations, 
           thin = thin, warmup = burn)

save(fit, file = paste0("results/bayes_models/trait_models/model_", response, "_froh.RData"))

