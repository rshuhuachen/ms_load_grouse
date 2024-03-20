### in this script we'll sum all the gerp scores across the different scaffolds (minus sex chr scaffold 4)

pacman::p_load(dplyr, data.table, lme4)

theme_set(theme_classic())

### load in data based on different gerp categories
gerps_data <- list.files(path="output/gerp/scores", pattern = "load_per_id_cat", full.names=T)

gerps_cat <- data.frame()
for (i in 1:length(gerps_data)){
  scaf <- fread(gerps_data[i])
  gerps_cat <- rbind(gerps_cat, scaf)
}

sum_gerp_cat <- gerps_cat %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:gerp_count_t, sum))

sum_gerp_cat <- sum_gerp_cat %>% mutate(
  gerp_Lp = gerp_count_p / n_genotyped,
  gerp_Lr = gerp_count_r / n_genotyped,
  gerp_Lt = gerp_count_t / n_genotyped)

### load in data based on different gerp categories but only additive total load
### load in data based on different gerp categories
gerps_data_additive <- list.files(path="output/gerp/scores", pattern = "load_per_id_cat_total_additive", full.names=T)

gerps_cat_additive <- data.frame()
for (i in 1:length(gerps_data_additive)){
  scaf <- fread(gerps_data_additive[i])
  gerps_cat_additive <- rbind(gerps_cat_additive, scaf)
}

sum_gerp_cat_additive <- gerps_cat_additive %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:gerp_count_t_v2, sum))

sum_gerp_cat_additive <- sum_gerp_cat_additive %>% mutate(
  gerp_Lt_additive = gerp_count_t_v2 / n_genotyped)

write_tsv(sum_gerp_cat, file = "output/load/gerp/sum_gerp_cat.tsv", col_names=TRUE)
write_tsv(sum_gerp_cat_additive, file = "output/load/gerp/sum_gerp_cat_additive.tsv", col_names=TRUE)
