### Here, we will quantify whether the total load based on GERP >= 4 is indeed biologically more deleterious compared to those with GERP <0
### however, to estimate this while keeping the number of mutations equal, we will randomly draw a subset of neutral mutations


### load packages ####
pacman::p_load(tidyverse, data.table, lme4)


### load neutral gerp mutations ####
load(file = "output/load/gerp/gerps_neutral.RData")

### load function to calculate load ####
source("scripts/7_calculate_load/function_calculate_load.R")

load_0 <- calculate_load_gerp(gerp_snp, output_vcf = F, loadtype = "gerp0")
model_out <- model_load(load_0, 1)

random_draws <- function(geno, n_draws, n_mutations, file, method, emperical_beta){
  source("scripts/theme_ggplot.R")
  all_draws <- data.frame()
  
  if(method == "GERP"){
    for (i in 1:n_draws){
      draw <- geno[sample(nrow(geno), n_mutations),] #randomly draw snps
      ## load functions
      source("scripts/7_calculate_load/function_calculate_load.R")
      load <- calculate_load_gerp(draw, output_vcf = F, loadtype = "random_draw") #calculate load
      model_out <- model_load(load, i)
      model_out$method <- method
      
      all_draws <- rbind(all_draws, model_out)
    }}
  
  if(method == "High impact"){
    for (i in 1:n_draws){
      draw <- geno[sample(nrow(geno), n_mutations),] #randomly draw snps
      ## load functions
      source("scripts/7_calculate_load/function_calculate_load.R")
      load <- calculate_load_snpeff(draw, output_vcf = F, loadtype = "random_draw") #calculate load
      model_out <- model_load(load, i)
      model_out$method <- method
      
      all_draws <- rbind(all_draws, model_out)
    }}
  
  ### conclusion
  all_draws <- all_draws %>% mutate(conclusion = as.factor(case_when(
    beta < 0 & pval < 0.05 ~ "Significantly negative",
    beta > 0 & pval < 0.05 ~ "Significantly positive",
    TRUE ~ "Insignificant"
  )))
  
  return(all_draws)
  save(all_draws, file = file)
}

# run functions

# all mutations
draws_gerp0 <- random_draws(geno = gerp_snp, n_draws = 5000, n_mutations = 413489, file = "output/random_draws/all_gerp0.RData", method="GERP", emperical_beta = -0.21)
save(draws_gerp0, file = "output/random_draws/all_gerp0.RData")