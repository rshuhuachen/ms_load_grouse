### Here we will combine genetic loads calculated with GERP and SnpEff ###

# load packages
pacman::p_load(tidyverse)

# load data

sum_gerp_cat <- read_tsv("output/load/gerp/sum_gerp_cat.tsv")
sum_gerp_cat_additive <- read_tsv("output/load/gerp/sum_gerp_cat_additive.tsv")
high <- read_tsv("output/load/snpeff/genetic_load_da_nosex_30scaf_HIGH.tsv")
high_add <- read_tsv("output/load/snpeff/genetic_load_da_nosex_total_additive_30scaf_HIGH.tsv")


#select relevant columns and rename before merging

# first make gerp scores into a wide dataframe

cat1 <- filter(sum_gerp_cat, gerp_cat == "0-1")
cat2 <- filter(sum_gerp_cat, gerp_cat == "1-2")
cat3 <- filter(sum_gerp_cat, gerp_cat == "2-3")
cat4 <- filter(sum_gerp_cat, gerp_cat == "3-4")
cat5 <- filter(sum_gerp_cat, gerp_cat == "4-5")

sum_gerp_cat_wide <- cat1[,c("id", "gerp_Lp", "gerp_Lr", "gerp_Lt")] %>% left_join(cat2[,c("id", "gerp_Lp", "gerp_Lr", "gerp_Lt")], by = "id", suffix = c("", "_cat2")) %>%
  left_join(cat3[,c("id", "gerp_Lp", "gerp_Lr", "gerp_Lt")], by = "id", suffix = c("", "_cat3")) %>%
  left_join(cat4[,c("id", "gerp_Lp", "gerp_Lr", "gerp_Lt")], by = "id", suffix = c("", "_cat4")) %>%
  left_join(cat5[,c("id", "gerp_Lp", "gerp_Lr", "gerp_Lt")], by = "id", suffix = c("", "_cat5"))

names(sum_gerp_cat_wide)[2] <- "gerp_Lp_cat1"
names(sum_gerp_cat_wide)[3] <- "gerp_Lr_cat1"
names(sum_gerp_cat_wide)[4] <- "gerp_Lt_cat1"

## make sum_gerp_cat_additive wide ###
cat1_add <- filter(sum_gerp_cat_additive, gerp_cat == "0-1")
cat2_add <- filter(sum_gerp_cat_additive, gerp_cat == "1-2")
cat3_add <- filter(sum_gerp_cat_additive, gerp_cat == "2-3")
cat4_add <- filter(sum_gerp_cat_additive, gerp_cat == "3-4")
cat5_add <- filter(sum_gerp_cat_additive, gerp_cat == "4-5")

sum_gerp_cat_additive_wide <- cat1_add[,c("id", "gerp_Lt_additive")] %>% 
  left_join(cat2_add[,c("id", "gerp_Lt_additive")], by = "id", suffix = c("", "_cat2")) %>%
  left_join(cat3_add[,c("id", "gerp_Lt_additive")], by = "id", suffix = c("", "_cat3")) %>%
  left_join(cat4_add[,c("id", "gerp_Lt_additive")], by = "id", suffix = c("", "_cat4")) %>%
  left_join(cat5_add[,c("id", "gerp_Lt_additive")], by = "id", suffix = c("", "_cat5"))

names(sum_gerp_cat_additive_wide)[2] <- "gerp_Lt_additive_cat1"

#### combine all

loads <- left_join(high[,c("id", "load_p", "load_r", "load_t")], high_add[,c("id", "load_t_v2")], by = "id") %>%
  left_join(sum_gerp_cat_wide, by = "id") %>%
  left_join(sum_gerp_cat_additive_wide, by = "id")

loads <- loads %>% rename(
  high_Lp = load_p,
  high_Lr = load_r,
  high_Lt = load_t,
  high_Lt_add = load_t_v2)

save(loads, file = "output/load/all_loads_combined_da_nosex_29scaf.RData")
write_tsv(loads, file = "output/load/all_loads_combined_da_nosex_29scaf.tsv", col_names=T)
