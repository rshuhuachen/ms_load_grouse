#combine lifetime fitness data and load
pacman::p_load(dplyr)

#lifetime
#load data
load("data/phenotypes/phenotypes_lifetime.RData") #lifetime phenotypes
load("output/load/all_loads_combined_da_nosex_30scaf_plus_per_region.RData") #loads no sex chr only 30 scaf
load("output/inbreeding/froh_combined.RData")
# combine
pheno_wide_load <- left_join(pheno_wide, load_per_region, by = "id")
pheno_wide_load <- left_join(pheno_wide_load, combined[c("id", "froh_bcftools")], by = "id")

save(pheno_wide_load, file = "output/load/pheno_loads_lifetime.RData")

annual
#load data
load("data/phenotypes/phenotypes_annual.RData") #annual

# combine
pheno_annual_load <- left_join(pheno_long, load_per_region, by = "id")
pheno_annual_load <- left_join(pheno_annual_load, combined[c("id", "froh_bcftools")], by = "id")

save(pheno_annual_load, file = "output/load/pheno_loads_annual.RData")
