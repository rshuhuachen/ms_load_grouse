#### Simplify LMS parameter here for modelling ####

# load packages
pacman::p_load(tidyverse, brms, bayesplot, data.table, ggridges)

source("scripts/theme_ggplot.R")

# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load_high <- load_per_region %>% filter(loadtype == "high")
load_gerp45 <- load_per_region %>% filter(loadtype == "gerp45")

load <- left_join(load_high, load_gerp45, by = "id", suffix = c("_high", "_gerp45"))

# combine
pheno_wide_load <- left_join(pheno_wide, load, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load_high)) #some ids without genotypes, excluded for wgr


##### model ####
brm_load_core_high <- brm(LMS_min ~ scale(total_load_high) + (1|born) + (1|site), data = subset(pheno_wide_load, core == "core"),
                            family = "zero_inflated_poisson",
                            prior = prior(normal(0,1), class = b),
                            cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                            iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_core_high, file = "output/models/total_hom_het/lms_total_core_only_high.RData")

brm_load_core_gerp45 <- brm(LMS_min ~ scale(total_load_gerp45) + (1|born) + (1|site), data = subset(pheno_wide_load, core == "core"),
                            family = "zero_inflated_poisson",
                              prior = prior(normal(0,1), class = b),
                              cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                              iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_core_gerp45, file = "output/models/total_hom_het/lms_total_core_only_gerp45.RData")

##### plot ####
load(file = "output/models/total_hom_het/lms_total_core_only_high.RData")
load(file = "output/models/total_hom_het/lms_total_core_only_gerp45.RData")

# get intervals
brm_load_core_gerp45_interval <- mcmc_intervals_data(brm_load_core_gerp45, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load_gerp45")
brm_load_core_high_interval <- mcmc_intervals_data(brm_load_core_high, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load_high")

brms_plota_interval <- rbind(brm_load_core_gerp45_interval, brm_load_core_high_interval)

brms_plota_interval$model <- c("Total GERP load", "Total SnpEff load")

data.frame(parameter = brms_plota_interval$parameter,
           median = round(brms_plota_interval$m, 2),
           ci_95 = paste0(round(brms_plota_interval$ll, 2), ", ", round(brms_plota_interval$hh, 2)),
           ci_80 = paste0(round(brms_plota_interval$l, 2), ", ", round(brms_plota_interval$h, 2)))

interval_core_gerp <- mcmc_intervals_data(brm_load_core_gerp45, prob =0.8, prob_outer = 0.95)
interval_core_gerp <- data.frame(parameter = interval_core_gerp$parameter,
                                   median = round(interval_core_gerp$m, 2),
                                   ci_95 = paste0(round(interval_core_gerp$ll, 2), ", ", round(interval_core_gerp$hh, 2)),
                                   ci_80 = paste0(round(interval_core_gerp$l, 2), ", ", round(interval_core_gerp$h, 2)))

write_tsv(interval_core_gerp, file = "output/models/intervals/total_gerp_core_full.tsv")

interval_core_high <- mcmc_intervals_data(brm_load_core_high, prob =0.8, prob_outer = 0.95)
interval_core_high <- data.frame(parameter = interval_core_high$parameter,
                                 median = round(interval_core_high$m, 2),
                                 ci_95 = paste0(round(interval_core_high$ll, 2), ", ", round(interval_core_high$hh, 2)),
                                 ci_80 = paste0(round(interval_core_high$l, 2), ", ", round(interval_core_high$h, 2)))

write_tsv(interval_core_high, file = "output/models/intervals/total_high_core_full.tsv")

# get areas
brm_load_core_gerp45_area <- mcmc_areas_data(brm_load_core_gerp45, pars = "b_scaletotal_load_gerp45")
brm_load_core_high_area <- mcmc_areas_data(brm_load_core_high, pars = "b_scaletotal_load_high")

brms_plota_areas <- rbind(brm_load_core_gerp45_area,
                          brm_load_core_high_area)

brms_plota_areas$model <- c( rep("Total GERP load", nrow(brm_load_core_gerp45_area)),
                             rep("Total SnpEff load", nrow(brm_load_core_high_area)))


#rearrange order for visualization
brms_plota_interval$model  <- factor(as.factor(brms_plota_interval$model),
                                     levels= c("Total SnpEff load", "Total GERP load"))

brms_plota_areas$model  <- factor(as.factor(brms_plota_areas$model),
                                  levels= c("Total SnpEff load", "Total GERP load"))

### plot

# split by interval
brms_plota <- split(brms_plota_areas, brms_plota_areas$interval)

brms_plota$bottom <- brms_plota$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_plota$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_plota_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_plota_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_plota_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "Parameter")+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp, clr_froh), 0.7)) +
  scale_y_discrete(labels = c("Total SnpEff load", "Total GERP load", 
                              expression(italic(F)[ROH])))+
  scale_color_manual(values =c(clr_high, clr_gerp, clr_froh)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> loads_lms_core

loads_lms_core

### simplify LMS into three categories ####
pheno_wide_load <- pheno_wide_load %>% mutate(LMS_cat = case_when(LMS_min > 5 ~ "High",
                                                                  LMS_min > 0 & LMS_min <= 5 ~ "Moderate",
                                                                  LMS_min == 0 ~ "Low"))

pheno_wide_load$LMS_cat <- factor(pheno_wide_load$LMS_cat, levels = c("Low", "Moderate", "High"))

##### model ####
brm_load_cat_high <- brm(LMS_cat ~ scale(total_load_high) + core + (1|site), data = pheno_wide_load,
                         family = "categorical",
                         cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                         iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_cat_high, file = "output/models/total_hom_het/lms_total_cat_high.RData")

brm_load_cat_gerp45 <- brm(LMS_cat ~ scale(total_load_gerp45) + core + (1|site), data = pheno_wide_load,
                           family = "categorical",
                           cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                           iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_cat_gerp45, file = "output/models/total_hom_het/lms_total_cat_gerp45.RData")

##### plot ####
load(file = "output/models/total_hom_het/lms_total_cat_high.RData")
summary(brm_load_cat_high)

# get intervals
brm_load_cat_gerp45_interval <- mcmc_intervals_data(brm_load_cat_gerp45, prob =0.8, prob_outer = 0.95)
brm_load_cat_gerp45_interval<- subset(brm_load_cat_gerp45_interval, grepl("*scaletotal_load*", parameter))
brm_load_cat_high_interval <- mcmc_intervals_data(brm_load_cat_high, prob =0.8, prob_outer = 0.95)
brm_load_cat_high_interval<- subset(brm_load_cat_high_interval, grepl("*scaletotal_load*", parameter))

brms_plota_interval <- rbind(brm_load_cat_gerp45_interval,
                             brm_load_cat_high_interval)

brms_plota_interval$model <- c("Total GERP load", "Total GERP load", "Total SnpEff load", "Total SnpEff load")
brms_plota_interval$lms_cat <- c("Moderate LMS", "High LMS", "Moderate LMS", "High LMS")

# export full intervals for gerp load

# get areas
brm_load_cat_gerp45_area <- mcmc_areas_data(brm_load_cat_gerp45)
brm_load_cat_gerp45_area<- subset(brm_load_cat_gerp45_area, grepl("*scaletotal_load*", parameter))
brm_load_cat_high_area <- mcmc_areas_data(brm_load_cat_high)
brm_load_cat_high_area<- subset(brm_load_cat_high_area, grepl("*scaletotal_load*", parameter))

brms_plota_areas <- rbind(brm_load_cat_gerp45_area,
                          brm_load_cat_high_area)

brms_plota_areas$model <- c(rep("Total GERP load", nrow(brm_load_cat_gerp45_area)),
                            rep("Total SnpEff load", nrow(brm_load_cat_high_area)))

brms_plota_areas$lms_cat <- c(rep("Moderate LMS", 0.5*nrow(brm_load_cat_gerp45_area)),
                              rep("High LMS", 0.5*nrow(brm_load_cat_gerp45_area)),
                              rep("Moderate LMS", 0.5*nrow(brm_load_cat_high_area)),
                              rep("High LMS", 0.5*nrow(brm_load_cat_high_area)))


# rearrange order for visualization
brms_plota_interval$model  <- factor(as.factor(brms_plota_interval$model),
                                     levels= c("Total SnpEff load", "Total GERP load"))

brms_plota_areas$model  <- factor(as.factor(brms_plota_areas$model),
                                  levels= c("Total SnpEff load", "Total GERP load"))

brms_plota_interval$lms_cat  <- factor(as.factor(brms_plota_interval$lms_cat),
                                       levels= c("Moderate LMS", "High LMS"))

brms_plota_areas$lms_cat  <- factor(as.factor(brms_plota_areas$lms_cat),
                                    levels= c("Moderate LMS", "High LMS"))

### plot

# split by interval
brms_plota <- split(brms_plota_areas, brms_plota_areas$interval)

brms_plota$bottom <- brms_plota$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_plota$outer) +  
  aes(x = .data$x, y = .data$lms_cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_plota_interval, aes(x = l, xend = h, yend = lms_cat), col = "black", linewidth=3)+
  geom_segment(data=brms_plota_interval, aes(x = ll, xend = hh, yend = lms_cat), col = "black")+
  geom_point(data=brms_plota_interval, aes(x = m, y = lms_cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "Parameter")+
  facet_grid(~model)+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp), 0.7)) +
  # scale_y_discrete(labels = c("Total SnpEff load", "Total GERP load"))+
  scale_color_manual(values =c(clr_high, clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> loads_lms_cat

#### combine in plot ####
cowplot::plot_grid(loads_lms_core,  loads_lms_cat, 
                   ncol = 1, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig_lms_simple

ggsave(fig_lms_simple, file = "plots/sup/lms_simplified.png", width=12, height=14)

