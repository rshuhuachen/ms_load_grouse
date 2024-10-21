library(dplyr); library(data.table); library(readr); library(brms)

args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
response <- args[[2]]
out <- args[[3]]
iterations <- args[[4]]
burn <- args[[5]]
thin <- args[[6]]

### load data ###

load(file = "data/phenotypes/phenotypes_annual.RData")

# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load <- load_per_region %>% filter(loadtype == method)

# merge files
pheno <- left_join(pheno_long, load, by = "id")
pheno$born <- pheno$year - pheno$age

pheno <- pheno %>% mutate(age_cat = as.factor(case_when(age == 1 ~ "yearling", age > 1  ~ "adult")))

### modelling ####
formula <- formula(paste0("scale(", response, ") ~ scale(total_load) + age_cat + (1|year) + (1|site/id)"))

fit <- brm(formula,
            family = "gaussian",
           data = pheno, 
           cores =8,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,1), class = b),
           iter = iterations, 
           thin = thin, warmup = burn)

save(fit, file = out)

