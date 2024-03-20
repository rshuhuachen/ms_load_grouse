library(dplyr); library(data.table); library(readr); library(brms)

args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
loadtype <- args[[2]]
region <- args[[3]]
response <- args[[4]]
iterations <- args[[5]]
burn <- args[[6]]
thin <- args[[7]]

predictor <- paste0(method, "_load_", loadtype, "_", region)


### load data ###

load(file = "output/load/pheno_loads_annual.RData")
pheno_core <- subset(pheno_annual_load, id %in% core$id)

### rename some cols
pheno_core <- pheno_core %>%
  rename(gerp5_load_t_add_all = gerp_Lt_additive_cat5,
         gerp5_load_e_all = gerp_Lr_cat5,
         gerp5_load_m_all = gerp_Lp_cat5,
         high_load_t_add_all = high_Lt_additive,
         high_load_e_all = high_Lr,
         high_load_m_all = high_Lp,)

names(pheno_core) <- gsub("t_add", "tadd", names(pheno_core))
#annual measurements for core males only

### modelling ####
formula <- formula(paste0("scale(", response, ") ~ scale(", predictor, ") + age + sqrt(age) + (1|year) + (1|site/id)"))

fit <- brm(formula,
             family = "gaussian",
           data = pheno_core, 
           cores =8,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,1), class = b),
           iter = iterations, 
           thin = thin, warmup = burn)

save(fit, file = paste0("results/bayes_models/trait_models/model_", response, "_", predictor, ".RData"))

