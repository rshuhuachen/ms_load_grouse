### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, 
               extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

#### plot sup snpeff lof and missense ####
load("output/models/total_hom_het/lms_total_lof.RData")
brm_load_t_add_lof_lms <- brm_load_t
load("output/models/total_hom_het/lms_total_missense.RData")
brm_load_t_add_missense_lms <- brm_load_t

# get intervals
brms_lof_interval <- mcmc_intervals_data(brm_load_t_add_lof_lms, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
brms_mis_interval <- mcmc_intervals_data(brm_load_t_add_missense_lms, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")

brms_othersnpeff_interval <- rbind(brms_lof_interval,
                                   brms_mis_interval)

brms_othersnpeff_interval$model <- c("Loss of Function", "Missense")

# get areas
brms_lof_area <- mcmc_areas_data(brm_load_t_add_lof_lms, pars = "b_scaletotal_load")
brms_mis_area <- mcmc_areas_data(brm_load_t_add_missense_lms, pars = "b_scaletotal_load")

brms_othersnpeff_areas <- rbind(brms_lof_area,
                                brms_mis_area)

brms_othersnpeff_areas$model <- c(rep("Loss of Function", nrow(brms_lof_area)),
                                  rep("Missense", nrow(brms_mis_area)))

#rearrange order for visualization
brms_othersnpeff_interval$model  <- factor(as.factor(brms_othersnpeff_interval$model),
                                           levels= c("Missense", "Loss of Function"))

brms_othersnpeff_areas$model  <- factor(as.factor(brms_othersnpeff_areas$model),
                                        levels= c("Missense", "Loss of Function"))

#plot

# split by interval
brms_othersnpeff <- split(brms_othersnpeff_areas, brms_othersnpeff_areas$interval)

brms_othersnpeff$bottom <- brms_othersnpeff$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_othersnpeff$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_othersnpeff_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_othersnpeff_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_othersnpeff_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "SnpEff category")+
  scale_fill_manual(values =alpha(c(clr_high, clr_high), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_high)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> other_snpeff


png(file = "plots/sup/other_snpeff.png", height = 600, width = 600)
other_snpeff
dev.off()
brms_othersnpeff_interval <- brms_othersnpeff_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_othersnpeff_interval, file = "output/models/intervals/other_snpeff.csv", quote=F, row.names = F)

### GERP categories ####
load(file = "output/models/other_gerp_cats/lms_total_gerp_0.RData")
brm_load_gerp0 <- brm_load_t
load(file = "output/models/other_gerp_cats/lms_total_gerp_01.RData")
brm_load_gerp01 <- brm_load_t
load(file = "output/models/other_gerp_cats/lms_total_gerp_12.RData")
brm_load_gerp12 <- brm_load_t
load(file = "output/models/other_gerp_cats/lms_total_gerp_23.RData")
brm_load_gerp23 <- brm_load_t
load(file = "output/models/other_gerp_cats/lms_total_gerp_34.RData")
brm_load_gerp34 <- brm_load_t

# load(file = "output/models/total_hom_het/lms_total_gerp34.RData")
# brm_load_t_add_gerp4_lms <- brm_load_t
# r2_bayes(brm_load_t_add_gerp4_lms)
load(file = "output/models/total_hom_het/lms_total_gerp45.RData")
brm_load_gerp4 <- brm_load_t

rm(brm_load_t)

#extract intervals and areas
#total
gerp0_t_interval <- mcmc_intervals_data(brm_load_gerp0, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp01_t_interval <- mcmc_intervals_data(brm_load_gerp01, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp12_t_interval <- mcmc_intervals_data(brm_load_gerp12, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp23_t_interval <- mcmc_intervals_data(brm_load_gerp23, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp34_t_interval <- mcmc_intervals_data(brm_load_gerp34, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp4_t_interval <- mcmc_intervals_data(brm_load_gerp4, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scaletotal_load", parameter))

brms_gerps_interval <- rbind(gerp0_t_interval,
                             gerp01_t_interval,
                             gerp12_t_interval,
                             gerp23_t_interval, 
                             gerp34_t_interval, 
                             gerp4_t_interval)

brms_gerps_interval$cat <- c("< 0", "0-1","1-2","2-3","3-4", "≥ 4")

#total

gerp0_t_area <- mcmc_areas_data(brm_load_gerp0) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp01_t_area <- mcmc_areas_data(brm_load_gerp01) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp12_t_area <- mcmc_areas_data(brm_load_gerp12) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp23_t_area <- mcmc_areas_data(brm_load_gerp23) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp34_t_area <- mcmc_areas_data(brm_load_gerp34) %>%
  subset(grepl("b_scaletotal_load", parameter))

gerp4_t_area <- mcmc_areas_data(brm_load_gerp4) %>%
  subset(grepl("b_scaletotal_load", parameter))

brms_gerps_area <- rbind(gerp0_t_area,
                         gerp01_t_area,
                         gerp12_t_area,
                         gerp23_t_area,
                         gerp34_t_area,
                         gerp4_t_area)

obs <- nrow(gerp01_t_area)
brms_gerps_area$cat <- rep(c("< 0", "0-1","1-2","2-3","3-4", "≥ 4"), each = obs)

#rearrange order for visualization
brms_gerps_area$cat  <- factor(as.factor(brms_gerps_area$cat),
                               levels= c("< 0", "0-1","1-2","2-3","3-4", "≥ 4"))

### plot

# split by interval
brms_gerps <- split(brms_gerps_area, brms_gerps_area$interval)

brms_gerps$bottom <- brms_gerps$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_gerps$outer) +  
  aes(x = .data$x, y = .data$cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = cat, col = cat))+
  geom_segment(data=brms_gerps_interval, aes(x = l, xend = h, yend = cat), col = "black", linewidth=3)+
  geom_segment(data=brms_gerps_interval, aes(x = ll, xend = hh, yend = cat), col = "black")+
  geom_point(data=brms_gerps_interval, aes(x = m, y = cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "GERP score category")+
  scale_fill_manual(values =alpha(c("#E5988A", "#E5988A","#E5988A","#E5988A","#E5988A",clr_gerp), 0.7)) +
  scale_color_manual(values =c("#E5988A", "#E5988A","#E5988A","#E5988A","#E5988A",clr_gerp)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_gerpcats
supp_posterior_gerpcats

png(file = "plots/sup/other_gerp.png", height = 600, width = 600)
supp_posterior_gerpcats
dev.off()
brms_gerps_interval <- brms_gerps_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_gerps_interval, file = "output/models/intervals/other_gerp.csv", quote=F, row.names = F)

