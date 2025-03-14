### packages ####
pacman::p_load(tidyverse, data.table, ggforce)

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

save(loads_per_scaf_window, file = "output/load_per_window/loads_per_scaf_window.RData")


### plot ####
source("scripts/theme_ggplot.R")

loads_per_scaf_window$window_nr <- as.numeric(loads_per_scaf_window$window_nr)
loads_per_scaf_window$scaf <- as.numeric(loads_per_scaf_window$scaf)
loads_per_scaf_window <- loads_per_scaf_window %>% mutate(sig = case_when(
  pval < 0.05 & est > 0 ~ "sig pos",
  pval < 0.05 & est < 0 ~ "sig neg",
  TRUE ~ "nonsig"
))

loads_per_scaf_window$qval <- p.adjust(loads_per_scaf_window$pval, method="fdr", n = nrow(loads_per_scaf_window))
loads_per_scaf_window <- loads_per_scaf_window %>% mutate(sig_q = case_when(
  qval < 0.05 & est > 0 ~ "sig pos",
  qval < 0.05 & est < 0 ~ "sig neg",
  TRUE ~ "nonsig"
))

loads_per_scaf_window %>% 
  filter(scaf < 10) %>%
  ggplot(aes(x = window_nr, y = est)) + 
  geom_line()+
  geom_hline(yintercept = 0, col = "darkred", linetype = "dotted")+
  geom_point(aes(col = sig_q), size = 2) + 
 # geom_segment(aes(y = est - se, yend = est + se, x = window_nr)) +
  scale_color_manual(values = c(clr_grey, "#703D57",  "#FFCD70")) +
  facet_wrap(scaf~., scales="free_x", ncol = 1) -> plot_scaf_10
ggsave(plot_scaf_10, file = "plots/sup/est_across_windows_gerp.png", width=14, height=18)
