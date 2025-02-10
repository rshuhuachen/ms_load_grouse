##### Bayesian models for LMS: total load for other gerp cats #####

pacman::p_load(bayesplot, brms, dplyr, data.table)

args <- commandArgs(trailingOnly = TRUE)
gerpcat <- args[[1]]
out <- args[[2]]
iter <- args[[3]]
warm <- args[[4]]
thin <- args[[5]]

if(gerpcat != "0"){
  gerpcat = paste0(substr(gerpcat, 1, 1), "-", substr(gerpcat, 2, 2))
}


if(gerpcat == "0"){
  gerpcat = "< 0"
}

#### load in all gerps ####

gerp_scafs <- list.files(path = "output/gerp", pattern = "count_mutations*", full.names = T)

gerp_raw <- data.frame()
for (i in 1:length(gerp_scafs)){
  scaf <- fread(gerp_scafs[i])
  gerp_raw <- rbind(gerp_raw, scaf)
}

gerp <- gerp_raw %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:hom_data, sum))

gerp <- gerp %>% mutate(het_load = het_data / n_genotyped,
                        hom_load = hom_data / n_genotyped,
                        total_load = (het_data * 0.5 + hom_data) / n_genotyped)


#### extract the gerp category ####

gerp_cat <- subset(gerp, gerp_cat == gerpcat)

# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# combine
pheno_wide_load <- left_join(pheno_wide, gerp_cat, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load)) #some ids without genotypes, excluded for wgr

#### model ####
brm_load_t <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = iter, thin = thin, warmup = warm, seed = 1908)

save(brm_load_t, file = out)


