### packages ####
pacman::p_load(tidyverse, data.table)

### load in summaries per window ####

sum_windows <- list.files(path = "output/load_per_window", pattern = "summary_load_per_window_scaf*", full.names = T)

scafs_windows <- list()
for (i in 1:length(sum_windows)){
  load(file = sum_windows[i])
  scafs_windows[[i]] <- loads_per_scaf
}

### loop over loads to output model results
loads_per_scaf_window <- data.frame()
for (i in 1:length(scafs_windows)){
  for (j in 1:length(scafs_windows[[i]])){
    load <- as.data.frame(scafs_windows[[i]][[j]])
    

    # only if there were snps in the window
    if (nrow(load) > 1 & sum(load$n_total) > 0 & sum(load$total_load, na.rm=T) > 0){
      # merge with pheno
      load("data/phenotypes/phenotypes_lifetime.RData")
      load_pheno <- left_join(load, pheno_wide, by = "id")

      # model with glmmtmb
      model_total <- glmmTMB(LMS_min ~ scale(total_load) + core + (1|site),
                            family = "poisson", ziformula = ~1,
                            data = load_pheno)
      
      sum_total <- summary(model_total)
      coef_total <- as.data.frame(sum_total$coefficients$cond)
      summary_total <- data.frame(loadtype = load$loadtype[1],
                                  scaf = unlist(strsplit(load$loadtype[1], "_", fixed=T))[2],
                                  window_nr = unlist(strsplit(load$loadtype[1], "_", fixed=T))[4],
                                   n_snps_window = load$n_total[1],
                                   est = coef_total$Estimate[2],
                                   se = coef_total$`Std. Error`[2],
                                   zval = coef_total$`z value`[2],
                                   pval = coef_total$`Pr(>|z|)`[2],
                                   load= "Total")
      
      loads_per_scaf_window <- rbind(loads_per_scaf_window, summary_total)
    }
  }
}
