### load packages ###
library(dplyr); library(data.table); library(readr); library(brms)

### args ###
args <- commandArgs(trailingOnly = TRUE)
go <- args[[1]]
rg <- args[[2]]
out <- args[[3]]
iterations <- args[[4]]
burn <- args[[5]]
thin <- args[[6]]

#### model effects on sexual traits ####
### load phenotypes ####
load(file="data/phenotypes/phenotypes_lifetime.RData")
load(file = "output/biological_pathways/loads_per_go_term_per_region.RData")

data <- left_join(pheno_wide[,c("id", "site", "LMS_min", "core")], loads, by = "id")

### fit model ####
set.seed(1908)

### modelling ####
formula <- formula(paste0("LMS_min ~ scale(total_load) + core + (1|site)"))

print(go)
nrow(subset(data, loadtype == go & region == rg & method == "gerp"))
# 
# rg = "promo"
# loadtype = "immune"
# 
# summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), 
#                          family = "poisson", ziformula = ~1, data = subset(data, region == rg & loadtype == go & method == "gerp")))

fit <- brm(formula,
           family = "zero_inflated_poisson",
           data = subset(data, loadtype == go & method == "gerp" & region == rg), 
           cores =8,
           control = list(adapt_delta = 0.99, max_treedepth = 15),
           prior = prior(normal(0,1), class = b),
           iter = iterations, 
           thin = thin, warmup = burn)

save(fit, file = out)
