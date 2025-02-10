### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

#### Check result 1 million iterations normal models ####
### GERP categories total load ####

#### only total #####
load(file = "../genetic_load_ltet/results/bayes_models/lms_gerp0_total_add_zi_29scaf.RData")
load(file = "results/bayes_models/lms_gerp1_total_add_zi_29scaf.RData")
load(file = "results/bayes_models/lms_gerp2_total_add_zi_29scaf.RData")
load(file = "results/bayes_models/lms_gerp3_total_add_zi_29scaf.RData")
load(file = "results/bayes_models/lms_gerp4_total_add_zi_29scaf.RData")
load(file = "results/bayes_models/lms_gerp5_total_add_zi_29scaf.RData")
brm_load_t_add_gerp5_lms <- brm_load_t_add_high_lms
brm_load_rp_gerp5_lms <- brm_load_rp_high_lms
rm(brm_load_t_add_high_lms)
rm(brm_load_rp_high_lms)

#extract intervals and areas
#total
gerp0_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp0_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat0", parameter))

gerp1_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp1_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat1", parameter))

gerp2_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp2_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat2", parameter))

gerp3_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp3_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat3", parameter))

gerp4_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp4_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat4", parameter))

gerp5_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_gerps_interval <- rbind(gerp0_t_interval, 
                             gerp1_t_interval,
                             gerp2_t_interval,
                             gerp3_t_interval,
                             gerp4_t_interval,
                             gerp5_t_interval)

brms_gerps_interval$cat <- c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4")

#rearrange order for visualization
#brms_gerps_interval$cat  <- factor(as.factor(brms_gerps_interval$cat),
#                                   levels= c("≥ 4", "3-4", "2-3", "1-2", "0-1", "< 0"))

write_tsv(brms_gerps_interval, file="results/bayes_models/tables/allgerpcats_tem_lms_intervals.tsv")

#total
gerp0_t_area <- mcmc_areas_data(brm_load_t_add_gerp0_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat0", parameter))

gerp1_t_area <- mcmc_areas_data(brm_load_t_add_gerp1_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat1", parameter))

gerp2_t_area <- mcmc_areas_data(brm_load_t_add_gerp2_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat2", parameter))

gerp3_t_area <- mcmc_areas_data(brm_load_t_add_gerp3_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat3", parameter))

gerp4_t_area <- mcmc_areas_data(brm_load_t_add_gerp4_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat4", parameter))

gerp5_t_area <- mcmc_areas_data(brm_load_t_add_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))


brms_gerps_area <- rbind(gerp0_t_area, 
                         gerp1_t_area,
                         gerp2_t_area,
                         gerp3_t_area,
                         gerp4_t_area,
                         gerp5_t_area)

obs <- nrow(gerp0_t_area)
brms_gerps_area$cat <- rep(c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"), each = obs)

#rearrange order for visualization
brms_gerps_area$cat  <- factor(as.factor(brms_gerps_area$cat),
                               levels= c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"))

### plot

# split by interval
brms_gerps <- split(brms_gerps_area, brms_gerps_area$interval)

brms_gerps$bottom <- brms_gerps$outer %>%
  group_by(!!! groups) %>%
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
  labs(x = "Standardised beta coefficient", y = "GERP score category")+
  scale_fill_manual(values =alpha(c("#F7E2DE", "#F0C4BC","#E5988A","#DA6B58", clr_gerp,"#B83E28"), 0.7)) +
  scale_color_manual(values =c("#F7E2DE", "#F0C4BC","#E5988A","#DA6B58", clr_gerp,"#B83E28")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_gerpcats
supp_posterior_gerpcats


### Same for SnpEff low and moderate #####

load(file = "results/bayes_models/lms_snpef_low_total_add_region_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_moderate_total_add_region_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_high_total_add_region_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_modifier_total_add_region_1mil_zi_29scaf.RData")

r2_bayes(brm_load_t_add_low_lms)
r2_bayes(brm_load_t_add_high_lms)
r2_bayes(brm_load_t_add_mod_lms)
r2_bayes(brm_load_t_add_modifier_lms)

#extract intervals and areas
#total
low_t_interval <- mcmc_intervals_data(brm_load_t_add_low_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalelow_Lt_additive", parameter))

mod_t_interval <- mcmc_intervals_data(brm_load_t_add_mod_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalemod_Lt_additive", parameter))

high_t_interval <- mcmc_intervals_data(brm_load_t_add_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

modifier_t_interval <- mcmc_intervals_data(brm_load_t_add_modifier_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalemodifier_Lt_additive", parameter))

brms_snpef_interval <- rbind(modifier_t_interval,
                             low_t_interval, 
                             mod_t_interval,
                             high_t_interval)

brms_snpef_interval$cat <- c("Modifier", "Low", "Moderate", "High")

#rearrange order for visualization
#brms_gerps_interval$cat  <- factor(as.factor(brms_gerps_interval$cat),
#                                   levels= c("≥ 4", "3-4", "2-3", "1-2", "0-1", "< 0"))

write_tsv(brms_snpef_interval, file="results/bayes_models/tables/allsnpef_t_lms_intervals.tsv")

#total
low_t_area <- mcmc_areas_data(brm_load_t_add_low_lms) %>%
  subset(grepl("b_scalelow_Lt_additive", parameter))

mod_t_area <- mcmc_areas_data(brm_load_t_add_mod_lms) %>%
  subset(grepl("b_scalemod_Lt_additive", parameter))

high_t_area <- mcmc_areas_data(brm_load_t_add_high_lms) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

modifier_t_area <- mcmc_areas_data(brm_load_t_add_modifier_lms) %>%
  subset(grepl("b_scalemodifier_Lt_additive", parameter))

brms_snpef_area <- rbind(modifier_t_area,
                         low_t_area, 
                         mod_t_area,
                         high_t_area)

obs <- nrow(low_t_area)
brms_snpef_area$cat <- rep(c("Modifier", "Low", "Moderate", "High"), each = obs)

#rearrange order for visualization
brms_snpef_area$cat  <- factor(as.factor(brms_snpef_area$cat),
                               levels= c("Modifier", "Low", "Moderate", "High"))

### plot

# split by interval
brms_snpef <- split(brms_snpef_area, brms_snpef_area$interval)

brms_snpef$bottom <- brms_snpef$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_snpef$outer) +  
  aes(x = .data$x, y = .data$cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = cat, col = cat))+
  geom_segment(data=brms_snpef_interval, aes(x = l, xend = h, yend = cat), col = "black", linewidth=3)+
  geom_segment(data=brms_snpef_interval, aes(x = ll, xend = hh, yend = cat), col = "black")+
  geom_point(data=brms_snpef_interval, aes(x = m, y = cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "SnpEff impact category")+
  scale_fill_manual(values =alpha(c("#CBE8F6", "#8EA4CC", "#549F93", clr_high ), 0.7)) +
  scale_color_manual(values =c("#CBE8F6", "#8EA4CC", "#549F93", clr_high)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_snpef
supp_posterior_snpef

### Also add GERP quantitative measure ####
load(file = "results/bayes_models/lms_snpef_gerp_q_min0_total_add_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_gerp_q_min1_total_add_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_gerp_q_min2_total_add_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_gerp_q_min3_total_add_1mil_zi_29scaf.RData")
load(file = "results/bayes_models/lms_snpef_gerp_q_min4_total_add_1mil_zi_29scaf.RData")

r2_bayes(brm_load_t_add_gerp_q_min0_lms)
r2_bayes(brm_load_t_add_gerp_q_min1_lms)
r2_bayes(brm_load_t_add_gerp_q_min2_lms)
r2_bayes(brm_load_t_add_gerp_q_min3_lms)
r2_bayes(brm_load_t_add_gerp_q_min4_lms)

#extract intervals and areas
#total
gerp_q_0_interval <- mcmc_intervals_data(brm_load_t_add_gerp_q_min0_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_q_min0", parameter))

gerp_q_1_interval <- mcmc_intervals_data(brm_load_t_add_gerp_q_min1_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_q_min1", parameter))

gerp_q_2_interval <- mcmc_intervals_data(brm_load_t_add_gerp_q_min2_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_q_min2", parameter))

gerp_q_3_interval <- mcmc_intervals_data(brm_load_t_add_gerp_q_min3_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_q_min3", parameter))

gerp_q_4_interval <- mcmc_intervals_data(brm_load_t_add_gerp_q_min4_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_q_min4", parameter))


brms_gerp_q_interval <- rbind(gerp_q_0_interval,
                              gerp_q_1_interval, 
                              gerp_q_2_interval,
                              gerp_q_3_interval,
                              gerp_q_4_interval
)

brms_gerp_q_interval$cat <- c("GERP > 0", "GERP > 1", "GERP > 2", "GERP > 3", "GERP > 4")

write_tsv(brms_gerp_q_interval, file="results/bayes_models/tables/gerp_q_t_lms_intervals.tsv")

#total
gerp_q_0_area <- mcmc_areas_data(brm_load_t_add_gerp_q_min0_lms) %>%
  subset(grepl("b_scalegerp_Lt_q_min0", parameter))

gerp_q_1_area <- mcmc_areas_data(brm_load_t_add_gerp_q_min1_lms) %>%
  subset(grepl("b_scalegerp_Lt_q_min1", parameter))

gerp_q_2_area <- mcmc_areas_data(brm_load_t_add_gerp_q_min2_lms) %>%
  subset(grepl("b_scalegerp_Lt_q_min2", parameter))

gerp_q_3_area <- mcmc_areas_data(brm_load_t_add_gerp_q_min3_lms) %>%
  subset(grepl("b_scalegerp_Lt_q_min3", parameter))

gerp_q_4_area <- mcmc_areas_data(brm_load_t_add_gerp_q_min4_lms) %>%
  subset(grepl("b_scalegerp_Lt_q_min4", parameter))

brms_gerp_q_area <- rbind(gerp_q_0_area,
                          gerp_q_1_area, 
                          gerp_q_2_area,
                          gerp_q_3_area,
                          gerp_q_4_area
)

obs <- nrow(gerp_q_0_area)
brms_gerp_q_area$cat <- rep(c("GERP > 0", "GERP > 1", "GERP > 2", "GERP > 3",  "GERP > 4"), each = obs)

#rearrange order for visualization
# brms_gerp_q_area$cat  <- factor(as.factor(brms_gerp_q_area$cat),
#                                levels= c("GERP > 0", "GERP > 1", "GERP > 2", "GERP > 3", "GERP > 4"))

### plot

# split by interval
brms_gerp_q <- split(brms_gerp_q_area, brms_gerp_q_area$interval)

brms_gerp_q$bottom <- brms_gerp_q$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_gerp_q$outer) +  
  aes(x = .data$x, y = .data$cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = cat, col = cat))+
  geom_segment(data=brms_gerp_q_interval, aes(x = l, xend = h, yend = cat), col = "black", linewidth=3)+
  geom_segment(data=brms_gerp_q_interval, aes(x = ll, xend = hh, yend = cat), col = "black")+
  geom_point(data=brms_gerp_q_interval, aes(x = m, y = cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "GERP score threshold")+
  scale_fill_manual(values =alpha(c("#F0C4BC","#E5988A","#DA6B58", clr_gerp,"#B83E28"), 0.7)) +
  scale_color_manual(values =c("#F0C4BC","#E5988A","#DA6B58", clr_gerp,"#B83E28")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_gerp_quant


#### Add in one plot all alternative models gerp and snpeff ####
plot_grid(supp_posterior_gerpcats, supp_posterior_gerp_quant, supp_posterior_snpef,
          labels=c("a", "b", "c"),
          align = "hv", axis = "lb",
          ncol=2, label_fontface = "plain", label_size = 22) -> supp_othercats

ggsave(supp_othercats, file = "plots/final/sup_posteriors_other_cats.png", height=12, width = 12)


### Total expressed realized  for all gerp cats ####
#combine total, expressed and realized in one


#extract intervals and areas
#total
gerp1_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp1_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat1", parameter))

gerp2_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp2_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat2", parameter))

gerp3_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp3_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat3", parameter))

gerp4_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp4_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat4", parameter))

gerp5_t_interval <- mcmc_intervals_data(brm_load_t_add_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

#realised
gerp1_r_interval <- mcmc_intervals_data(brm_load_rp_gerp1_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat1", parameter))

gerp2_r_interval <- mcmc_intervals_data(brm_load_rp_gerp2_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat2", parameter))

gerp3_r_interval <- mcmc_intervals_data(brm_load_rp_gerp3_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat3", parameter))

gerp4_r_interval <- mcmc_intervals_data(brm_load_rp_gerp4_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat4", parameter))

gerp5_r_interval <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat5", parameter))

#potential
gerp1_p_interval <- mcmc_intervals_data(brm_load_rp_gerp1_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat1", parameter))

gerp2_p_interval <- mcmc_intervals_data(brm_load_rp_gerp2_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat2", parameter))

gerp3_p_interval <- mcmc_intervals_data(brm_load_rp_gerp3_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat3", parameter))

gerp4_p_interval <- mcmc_intervals_data(brm_load_rp_gerp4_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat4", parameter))

gerp5_p_interval <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat5", parameter))

brms_gerps_interval <- rbind(gerp1_t_interval,
                             gerp2_t_interval,
                             gerp3_t_interval,
                             gerp4_t_interval,
                             gerp5_t_interval,
                             gerp1_r_interval,
                             gerp2_r_interval,
                             gerp3_r_interval,
                             gerp4_r_interval,
                             gerp5_r_interval,
                             gerp1_p_interval,
                             gerp2_p_interval,
                             gerp3_p_interval,
                             gerp4_p_interval,
                             gerp5_p_interval)


brms_gerps_interval$cat <- c(rep(c("0-1", "1-2", "2-3", "3-4", "4-5"), times = 3))
brms_gerps_interval$loadtype <- c(rep(c("Total", "Expressed", "Masked"), each = 5))

#rearrange order for visualization
brms_gerps_interval$cat  <- factor(as.factor(brms_gerps_interval$cat),
                                   levels= c("4-5", "3-4", "2-3", "1-2", "0-1"))

brms_gerps_interval$loadtype  <- factor(as.factor(brms_gerps_interval$loadtype),
                                        levels= c( "Total", "Expressed","Masked" ))

write_tsv(brms_gerps_interval, file="results/bayes_models/tables/allgerpcats_tem_lms_intervals.tsv")

#areas

#total
gerp1_t_area <- mcmc_areas_data(brm_load_t_add_gerp1_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat1", parameter))

gerp2_t_area <- mcmc_areas_data(brm_load_t_add_gerp2_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat2", parameter))

gerp3_t_area <- mcmc_areas_data(brm_load_t_add_gerp3_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat3", parameter))

gerp4_t_area <- mcmc_areas_data(brm_load_t_add_gerp4_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat4", parameter))

gerp5_t_area <- mcmc_areas_data(brm_load_t_add_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

#realised
gerp1_r_area <- mcmc_areas_data(brm_load_rp_gerp1_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat1", parameter))

gerp2_r_area <- mcmc_areas_data(brm_load_rp_gerp2_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat2", parameter))

gerp3_r_area <- mcmc_areas_data(brm_load_rp_gerp3_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat3", parameter))

gerp4_r_area <- mcmc_areas_data(brm_load_rp_gerp4_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat4", parameter))

gerp5_r_area <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat5", parameter))

#potential
gerp1_p_area <- mcmc_areas_data(brm_load_rp_gerp1_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat1", parameter))

gerp2_p_area <- mcmc_areas_data(brm_load_rp_gerp2_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat2", parameter))

gerp3_p_area <- mcmc_areas_data(brm_load_rp_gerp3_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat3", parameter))

gerp4_p_area <- mcmc_areas_data(brm_load_rp_gerp4_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat4", parameter))

gerp5_p_area <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat5", parameter))

brms_gerps_area <- rbind(gerp1_t_area,
                         gerp2_t_area,
                         gerp3_t_area,
                         gerp4_t_area,
                         gerp5_t_area,
                         gerp1_r_area,
                         gerp2_r_area,
                         gerp3_r_area,
                         gerp4_r_area,
                         gerp5_r_area,
                         gerp1_p_area,
                         gerp2_p_area,
                         gerp3_p_area,
                         gerp4_p_area,
                         gerp5_p_area)

obs <- nrow(gerp1_p_area)
brms_gerps_area$cat <- c(rep(rep(c("0-1", "1-2", "2-3", "3-4", "4-5"), each = obs), times = 3))
brms_gerps_area$loadtype <- c(rep(rep(c("Total", "Expressed", "Masked"), each =obs), each = 5))

#rearrange order for visualization
brms_gerps_area$cat  <- factor(as.factor(brms_gerps_area$cat),
                               levels= c("4-5", "3-4", "2-3", "1-2", "0-1"))

brms_gerps_area$loadtype  <- factor(as.factor(brms_gerps_area$loadtype),
                                    levels= c( "Total", "Expressed","Masked" ))

### plot

# split by interval
brms_gerps <- split(brms_gerps_area, brms_gerps_area$interval)

brms_gerps$bottom <- brms_gerps$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = subset(brms_gerps$outer, cat != "0-1" & cat != "1-2")) +  
  aes(x = .data$x, y = .data$cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = cat, col = cat))+
  geom_segment(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"), aes(x = l, xend = h, yend = cat), col = "black", linewidth=3)+
  geom_segment(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"), aes(x = ll, xend = hh, yend = cat), col = "black")+
  geom_point(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"), aes(x = m, y = cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "GERP score category")+
  facet_wrap(~loadtype, ncol = 3, scales="free_x") +
  scale_fill_manual(values =alpha(c(clr_gerp, "#DE7A68", "#E9A69B"), 0.7)) +
  scale_color_manual(values =c(clr_gerp, "#DE7A68", "#E9A69B")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_gerpcats

ggsave(supp_posterior_gerpcats, file = "plots/final/supp_posteriors_gerpcats.png", width=12, height=10)


#### GERP per different region ####
#combine total, expressed and realized in one
load(file = "results/bayes_models/lms_total_gerp5_region_promoter_zi_29scaf.RData")
gerp5_t_promoter <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_tss_zi_29scaf.RData")
gerp5_t_tss <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_gene_zi_29scaf.RData")
gerp5_t_gene <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_down_zi_29scaf.RData")
gerp5_t_down <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_up_zi_29scaf.RData")
gerp5_t_up <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_exon_zi_29scaf.RData")
gerp5_t_exon <- brm_region
load(file = "results/bayes_models/lms_total_gerp5_region_intron_zi_29scaf.RData")
gerp5_t_intron <- brm_region

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_promoter_zi_29scaf.RData")
gerp5_rp_promoter <- brm_region
r2_bayes(gerp5_rp_promoter)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_tss_zi_29scaf.RData")
gerp5_rp_tss <- brm_region
r2_bayes(gerp5_rp_tss)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_gene_zi_29scaf.RData")
gerp5_rp_gene <- brm_region
r2_bayes(gerp5_rp_gene)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_down_zi_29scaf.RData")
gerp5_rp_down <- brm_region
r2_bayes(gerp5_rp_down)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_up_zi_29scaf.RData")
gerp5_rp_up <- brm_region
r2_bayes(gerp5_rp_up)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_exon_zi_29scaf.RData")
gerp5_rp_exon <- brm_region
r2_bayes(gerp5_rp_exon)

load(file = "results/bayes_models/lms_expressed_masked_gerp5_region_intron_zi_29scaf.RData")
gerp5_rp_intron <- brm_region
r2_bayes(gerp5_rp_intron)

rm(brm_region)

#extract intervals and areas
#total
gerp5_t_promoter_interval <- mcmc_intervals_data(gerp5_t_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_promoter", parameter))

gerp5_t_tss_interval <- mcmc_intervals_data(gerp5_t_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_tss", parameter))

gerp5_t_gene_interval <- mcmc_intervals_data(gerp5_t_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_gene", parameter))

gerp5_t_down_interval <- mcmc_intervals_data(gerp5_t_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_down", parameter))

gerp5_t_up_interval <- mcmc_intervals_data(gerp5_t_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_up", parameter))

gerp5_t_intron_interval <- mcmc_intervals_data(gerp5_t_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_intron", parameter))

gerp5_t_exon_interval <- mcmc_intervals_data(gerp5_t_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_t_add_exon", parameter))

#realised
gerp5_r_promoter_interval <- mcmc_intervals_data(gerp5_rp_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_promoter", parameter))

gerp5_r_tss_interval <- mcmc_intervals_data(gerp5_rp_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_tss", parameter))

gerp5_r_gene_interval <- mcmc_intervals_data(gerp5_rp_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_gene", parameter))

gerp5_r_down_interval <- mcmc_intervals_data(gerp5_rp_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_down", parameter))

gerp5_r_up_interval <- mcmc_intervals_data(gerp5_rp_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_up", parameter))

gerp5_r_intron_interval <- mcmc_intervals_data(gerp5_rp_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_intron", parameter))

gerp5_r_exon_interval <- mcmc_intervals_data(gerp5_rp_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_e_exon", parameter))


#potential
#realised
gerp5_p_promoter_interval <- mcmc_intervals_data(gerp5_rp_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_promoter", parameter))

gerp5_p_tss_interval <- mcmc_intervals_data(gerp5_rp_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_tss", parameter))

gerp5_p_gene_interval <- mcmc_intervals_data(gerp5_rp_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_gene", parameter))

gerp5_p_down_interval <- mcmc_intervals_data(gerp5_rp_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_down", parameter))

gerp5_p_up_interval <- mcmc_intervals_data(gerp5_rp_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_up", parameter))

gerp5_p_intron_interval <- mcmc_intervals_data(gerp5_rp_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_intron", parameter))

gerp5_p_exon_interval <- mcmc_intervals_data(gerp5_rp_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_m_exon", parameter))

brms_gerps_regions_interval <- rbind(gerp5_t_promoter_interval, gerp5_t_tss_interval, gerp5_t_gene_interval,
                                     gerp5_t_down_interval, gerp5_t_up_interval, 
                                     gerp5_t_intron_interval, gerp5_t_exon_interval,
                                     gerp5_r_promoter_interval, gerp5_r_tss_interval, gerp5_r_gene_interval,
                                     gerp5_r_down_interval, gerp5_r_up_interval, 
                                     gerp5_r_intron_interval, gerp5_r_exon_interval,
                                     gerp5_p_promoter_interval, gerp5_p_tss_interval, gerp5_p_gene_interval,
                                     gerp5_p_down_interval, gerp5_p_up_interval, 
                                     gerp5_p_intron_interval, gerp5_p_exon_interval)


brms_gerps_regions_interval$region <- c(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"),
                                            times = 3))
brms_gerps_regions_interval$loadtype <- c(rep(c("Additive", "Homozygous", "Heterozygous"), each = 7))
#brms_gerps_regions_interval$loadtype <- c(rep(c("Total", "Expressed", "Masked"), each = 7))

#rearrange order for visualization
brms_gerps_regions_interval$region  <- factor(as.factor(brms_gerps_regions_interval$region),
                                              levels= c("Exon", "Intron", "Upstream", "Downstream", "Gene body",
                                                        "TSS", "Promoter"))

#brms_gerps_regions_interval$loadtype  <- factor(as.factor(brms_gerps_regions_interval$loadtype),
#                                        levels= c( "Total", "Expressed","Masked" ))

brms_gerps_regions_interval$loadtype  <- factor(as.factor(brms_gerps_regions_interval$loadtype),
                                                levels= c( "Additive", "Homozygous","Heterozygous" ))

write_tsv(brms_gerps_regions_interval, file="results/bayes_models/tables/gerp_tem_regions_lms_intervals.tsv")

#areas

#total
gerp5_t_promoter_area <- mcmc_areas_data(gerp5_t_promoter) %>%
  subset(grepl("b_scalegerp5_load_t_add_promoter", parameter))

gerp5_t_tss_area <- mcmc_areas_data(gerp5_t_tss) %>%
  subset(grepl("b_scalegerp5_load_t_add_tss", parameter))

gerp5_t_gene_area <- mcmc_areas_data(gerp5_t_gene) %>%
  subset(grepl("b_scalegerp5_load_t_add_gene", parameter))

gerp5_t_down_area <- mcmc_areas_data(gerp5_t_down) %>%
  subset(grepl("b_scalegerp5_load_t_add_down", parameter))

gerp5_t_up_area <- mcmc_areas_data(gerp5_t_up) %>%
  subset(grepl("b_scalegerp5_load_t_add_up", parameter))

gerp5_t_intron_area <- mcmc_areas_data(gerp5_t_intron) %>%
  subset(grepl("b_scalegerp5_load_t_add_intron", parameter))

gerp5_t_exon_area <- mcmc_areas_data(gerp5_t_exon) %>%
  subset(grepl("b_scalegerp5_load_t_add_exon", parameter))

#realised
gerp5_r_promoter_area <- mcmc_areas_data(gerp5_rp_promoter) %>%
  subset(grepl("b_scalegerp5_load_e_promoter", parameter))

gerp5_r_tss_area <- mcmc_areas_data(gerp5_rp_tss) %>%
  subset(grepl("b_scalegerp5_load_e_tss", parameter))

gerp5_r_gene_area <- mcmc_areas_data(gerp5_rp_gene) %>%
  subset(grepl("b_scalegerp5_load_e_gene", parameter))

gerp5_r_down_area <- mcmc_areas_data(gerp5_rp_down) %>%
  subset(grepl("b_scalegerp5_load_e_down", parameter))

gerp5_r_up_area <- mcmc_areas_data(gerp5_rp_up) %>%
  subset(grepl("b_scalegerp5_load_e_up", parameter))

gerp5_r_intron_area <- mcmc_areas_data(gerp5_rp_intron) %>%
  subset(grepl("b_scalegerp5_load_e_intron", parameter))

gerp5_r_exon_area <- mcmc_areas_data(gerp5_rp_exon) %>%
  subset(grepl("b_scalegerp5_load_e_exon", parameter))


#potential
#realised
gerp5_p_promoter_area <- mcmc_areas_data(gerp5_rp_promoter) %>%
  subset(grepl("b_scalegerp5_load_m_promoter", parameter))

gerp5_p_tss_area <- mcmc_areas_data(gerp5_rp_tss) %>%
  subset(grepl("b_scalegerp5_load_m_tss", parameter))

gerp5_p_gene_area <- mcmc_areas_data(gerp5_rp_gene) %>%
  subset(grepl("b_scalegerp5_load_m_gene", parameter))

gerp5_p_down_area <- mcmc_areas_data(gerp5_rp_down) %>%
  subset(grepl("b_scalegerp5_load_m_down", parameter))

gerp5_p_up_area <- mcmc_areas_data(gerp5_rp_up) %>%
  subset(grepl("b_scalegerp5_load_m_up", parameter))

gerp5_p_intron_area <- mcmc_areas_data(gerp5_rp_intron) %>%
  subset(grepl("b_scalegerp5_load_m_intron", parameter))

gerp5_p_exon_area <- mcmc_areas_data(gerp5_rp_exon) %>%
  subset(grepl("b_scalegerp5_load_m_exon", parameter))

brms_gerps_regions_area <- rbind(gerp5_t_promoter_area, gerp5_t_tss_area, gerp5_t_gene_area,
                                 gerp5_t_down_area, gerp5_t_up_area, 
                                 gerp5_t_intron_area, gerp5_t_exon_area,
                                 gerp5_r_promoter_area, gerp5_r_tss_area, gerp5_r_gene_area,
                                 gerp5_r_down_area, gerp5_r_up_area, 
                                 gerp5_r_intron_area, gerp5_r_exon_area,
                                 gerp5_p_promoter_area, gerp5_p_tss_area, gerp5_p_gene_area,
                                 gerp5_p_down_area, gerp5_p_up_area, 
                                 gerp5_p_intron_area, gerp5_p_exon_area)


obs <- nrow(gerp5_t_promoter_area)
brms_gerps_regions_area$region <- c(rep(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"), each = obs), times = 3))
#brms_gerps_regions_area$loadtype <- c(rep(rep(c("Total", "Expressed", "Masked"), each =obs), each = 7))
brms_gerps_regions_area$loadtype <- c(rep(rep(c("Additive", "Homozygous", "Heterozygous"), each =obs), each = 7))

#rearrange order for visualization
brms_gerps_regions_area$region  <- factor(as.factor(brms_gerps_regions_area$region),
                                          levels= c( "Intron", "Exon","Upstream", "Downstream", "Gene body",
                                                     "TSS", "Promoter"))

# brms_gerps_regions_area$loadtype  <- factor(as.factor(brms_gerps_regions_area$loadtype),
#                                                 levels= c( "Total", "Expressed","Masked" ))

brms_gerps_regions_area$loadtype  <- factor(as.factor(brms_gerps_regions_area$loadtype),
                                            levels= c( "Additive", "Homozygous","Heterozygous" ))
### plot

# split by interval
brms_gerps_regions <- split(brms_gerps_regions_area, brms_gerps_regions_area$interval)

brms_gerps_regions$bottom <- brms_gerps_regions$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_gerps_regions$outer) +  
  aes(x = .data$x, y = .data$region) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = loadtype, col = loadtype))+
  geom_segment(data=brms_gerps_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_gerps_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_gerps_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Region", title = "GERP score ≥ 4")+
  facet_wrap(~loadtype, ncol = 3, scales="free_x") +
  scale_fill_manual(values =alpha(c(clr_gerp, "#DE7A68", "#E9A69B"), 0.7)) +
  scale_color_manual(values =c(clr_gerp, "#DE7A68", "#E9A69B")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> supp_posterior_gerpregions

ggsave(supp_posterior_gerpregions, file = "plots/final/supp_posteriors_gerp_regions.png", width=12, height=16)

#### SNPeff per different region ####
#combine total, expressed and realized in one
load(file = "results/bayes_models/lms_total_high_region_promoter_withregion_zi_29scaf.RData")
high_t_promoter <- brm_region
load(file = "results/bayes_models/lms_total_high_region_tss_withregion_zi_29scaf.RData")
high_t_tss <- brm_region
load(file = "results/bayes_models/lms_total_high_region_gene_withregion_zi_29scaf.RData")
high_t_gene <- brm_region
load(file = "results/bayes_models/lms_total_high_region_down_withregion_zi_29scaf.RData")
high_t_down <- brm_region
load(file = "results/bayes_models/lms_total_high_region_up_withregion_zi_29scaf.RData")
high_t_up <- brm_region
load(file = "results/bayes_models/lms_total_high_region_exon_withregion_zi_29scaf.RData")
high_t_exon <- brm_region
load(file = "results/bayes_models/lms_total_high_region_intron_withregion_zi_29scaf.RData")
high_t_intron <- brm_region

load(file = "results/bayes_models/lms_expressed_masked_high_region_promoter_withregion_zi_29scaf.RData")
high_rp_promoter <- brm_region
r2_bayes(high_rp_promoter)

load(file = "results/bayes_models/lms_expressed_masked_high_region_tss_withregion_zi_29scaf.RData")
high_rp_tss <- brm_region
r2_bayes(high_rp_tss)

load(file = "results/bayes_models/lms_expressed_masked_high_region_gene_withregion_zi_29scaf.RData")
high_rp_gene <- brm_region
r2_bayes(high_rp_gene)

load(file = "results/bayes_models/lms_expressed_masked_high_region_down_withregion_zi_29scaf.RData")
high_rp_down <- brm_region
r2_bayes(high_rp_down)

load(file = "results/bayes_models/lms_expressed_masked_high_region_up_withregion_zi_29scaf.RData")
high_rp_up <- brm_region
r2_bayes(high_rp_up)

load(file = "results/bayes_models/lms_expressed_masked_high_region_exon_withregion_zi_29scaf.RData")
high_rp_exon <- brm_region
r2_bayes(high_rp_exon)

load(file = "results/bayes_models/lms_expressed_masked_high_region_intron_withregion_zi_29scaf.RData")
high_rp_intron <- brm_region
r2_bayes(high_rp_intron)

rm(brm_region)

#extract intervals and areas
#total
high_t_promoter_interval <- mcmc_intervals_data(high_t_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_promoter", parameter))

high_t_tss_interval <- mcmc_intervals_data(high_t_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_tss", parameter))

high_t_gene_interval <- mcmc_intervals_data(high_t_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_gene", parameter))

high_t_down_interval <- mcmc_intervals_data(high_t_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_down", parameter))

high_t_up_interval <- mcmc_intervals_data(high_t_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_up", parameter))

high_t_intron_interval <- mcmc_intervals_data(high_t_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_intron", parameter))

high_t_exon_interval <- mcmc_intervals_data(high_t_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_t_add_exon", parameter))

#realised
high_r_promoter_interval <- mcmc_intervals_data(high_rp_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_promoter", parameter))

high_r_tss_interval <- mcmc_intervals_data(high_rp_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_tss", parameter))

high_r_gene_interval <- mcmc_intervals_data(high_rp_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_gene", parameter))

high_r_down_interval <- mcmc_intervals_data(high_rp_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_down", parameter))

high_r_up_interval <- mcmc_intervals_data(high_rp_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_up", parameter))

high_r_intron_interval <- mcmc_intervals_data(high_rp_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_intron", parameter))

high_r_exon_interval <- mcmc_intervals_data(high_rp_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_e_exon", parameter))


#potential
#realised
high_p_promoter_interval <- mcmc_intervals_data(high_rp_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_promoter", parameter))

high_p_tss_interval <- mcmc_intervals_data(high_rp_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_tss", parameter))

high_p_gene_interval <- mcmc_intervals_data(high_rp_gene, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_gene", parameter))

high_p_down_interval <- mcmc_intervals_data(high_rp_down, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_down", parameter))

high_p_up_interval <- mcmc_intervals_data(high_rp_up, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_up", parameter))

high_p_intron_interval <- mcmc_intervals_data(high_rp_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_intron", parameter))

high_p_exon_interval <- mcmc_intervals_data(high_rp_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_load_m_exon", parameter))

brms_high_regions_interval <- rbind(high_t_promoter_interval, high_t_tss_interval, high_t_gene_interval,
                                    high_t_down_interval, high_t_up_interval, 
                                    high_t_intron_interval, high_t_exon_interval,
                                    high_r_promoter_interval, high_r_tss_interval, high_r_gene_interval,
                                    high_r_down_interval, high_r_up_interval, 
                                    high_r_intron_interval, high_r_exon_interval,
                                    high_p_promoter_interval, high_p_tss_interval, high_p_gene_interval,
                                    high_p_down_interval, high_p_up_interval, 
                                    high_p_intron_interval, high_p_exon_interval)


brms_high_regions_interval$region <- c(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"),
                                           times = 3))
#brms_high_regions_interval$loadtype <- c(rep(c("Total", "Expressed", "Masked"), each = 7))
brms_high_regions_interval$loadtype <- c(rep(c("Additive", "Homozygous", "Heterozygous"), each = 7))

#rearrange order for visualization
brms_high_regions_interval$region  <- factor(as.factor(brms_high_regions_interval$region),
                                             levels= c( "Intron", "Exon","Upstream", "Downstream", "Gene body",
                                                        "TSS", "Promoter"))

# brms_high_regions_interval$loadtype  <- factor(as.factor(brms_high_regions_interval$loadtype),
#                                                 levels= c( "Total", "Expressed","Masked" ))
brms_high_regions_interval$loadtype  <- factor(as.factor(brms_high_regions_interval$loadtype),
                                               levels= c( "Additive", "Homozygous","Heterozygous" ))

write_tsv(brms_high_regions_interval, file="results/bayes_models/tables/high_tem_regions_lms_intervals.tsv")

#areas

#total
high_t_promoter_area <- mcmc_areas_data(high_t_promoter) %>%
  subset(grepl("b_scalehigh_load_t_add_promoter", parameter))

high_t_tss_area <- mcmc_areas_data(high_t_tss) %>%
  subset(grepl("b_scalehigh_load_t_add_tss", parameter))

high_t_gene_area <- mcmc_areas_data(high_t_gene) %>%
  subset(grepl("b_scalehigh_load_t_add_gene", parameter))

high_t_down_area <- mcmc_areas_data(high_t_down) %>%
  subset(grepl("b_scalehigh_load_t_add_down", parameter))

high_t_up_area <- mcmc_areas_data(high_t_up) %>%
  subset(grepl("b_scalehigh_load_t_add_up", parameter))

high_t_intron_area <- mcmc_areas_data(high_t_intron) %>%
  subset(grepl("b_scalehigh_load_t_add_intron", parameter))

high_t_exon_area <- mcmc_areas_data(high_t_exon) %>%
  subset(grepl("b_scalehigh_load_t_add_exon", parameter))

#realised
high_r_promoter_area <- mcmc_areas_data(high_rp_promoter) %>%
  subset(grepl("b_scalehigh_load_e_promoter", parameter))

high_r_tss_area <- mcmc_areas_data(high_rp_tss) %>%
  subset(grepl("b_scalehigh_load_e_tss", parameter))

high_r_gene_area <- mcmc_areas_data(high_rp_gene) %>%
  subset(grepl("b_scalehigh_load_e_gene", parameter))

high_r_down_area <- mcmc_areas_data(high_rp_down) %>%
  subset(grepl("b_scalehigh_load_e_down", parameter))

high_r_up_area <- mcmc_areas_data(high_rp_up) %>%
  subset(grepl("b_scalehigh_load_e_up", parameter))

high_r_intron_area <- mcmc_areas_data(high_rp_intron) %>%
  subset(grepl("b_scalehigh_load_e_intron", parameter))

high_r_exon_area <- mcmc_areas_data(high_rp_exon) %>%
  subset(grepl("b_scalehigh_load_e_exon", parameter))


#potential
#realised
high_p_promoter_area <- mcmc_areas_data(high_rp_promoter) %>%
  subset(grepl("b_scalehigh_load_m_promoter", parameter))

high_p_tss_area <- mcmc_areas_data(high_rp_tss) %>%
  subset(grepl("b_scalehigh_load_m_tss", parameter))

high_p_gene_area <- mcmc_areas_data(high_rp_gene) %>%
  subset(grepl("b_scalehigh_load_m_gene", parameter))

high_p_down_area <- mcmc_areas_data(high_rp_down) %>%
  subset(grepl("b_scalehigh_load_m_down", parameter))

high_p_up_area <- mcmc_areas_data(high_rp_up) %>%
  subset(grepl("b_scalehigh_load_m_up", parameter))

high_p_intron_area <- mcmc_areas_data(high_rp_intron) %>%
  subset(grepl("b_scalehigh_load_m_intron", parameter))

high_p_exon_area <- mcmc_areas_data(high_rp_exon) %>%
  subset(grepl("b_scalehigh_load_m_exon", parameter))

brms_high_regions_area <- rbind(high_t_promoter_area, high_t_tss_area, high_t_gene_area,
                                high_t_down_area, high_t_up_area, 
                                high_t_intron_area, high_t_exon_area,
                                high_r_promoter_area, high_r_tss_area, high_r_gene_area,
                                high_r_down_area, high_r_up_area, 
                                high_r_intron_area, high_r_exon_area,
                                high_p_promoter_area, high_p_tss_area, high_p_gene_area,
                                high_p_down_area, high_p_up_area, 
                                high_p_intron_area, high_p_exon_area)


obs <- nrow(high_t_promoter_area)
brms_high_regions_area$region <- c(rep(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"), each = obs), times = 3))
#brms_high_regions_area$loadtype <- c(rep(rep(c("Total", "Expressed", "Masked"), each =obs), each = 7))
brms_high_regions_area$loadtype <- c(rep(rep(c("Additive", "Homozygous", "Heterozygous"), each =obs), each = 7))

#rearrange order for visualization
brms_high_regions_area$region  <- factor(as.factor(brms_high_regions_area$region),
                                         levels= c("Intron", "Exon",  "Upstream", "Downstream", "Gene body",
                                                   "TSS", "Promoter"))

# brms_high_regions_area$loadtype  <- factor(as.factor(brms_high_regions_area$loadtype),
#                                             levels= c( "Total", "Expressed","Masked" ))

brms_high_regions_area$loadtype  <- factor(as.factor(brms_high_regions_area$loadtype),
                                           levels= c( "Additive", "Homozygous","Heterozygous" ))

### plot

# split by interval
brms_high_regions <- split(brms_high_regions_area, brms_high_regions_area$interval)

brms_high_regions$bottom <- brms_high_regions$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_high_regions$outer) +  
  aes(x = .data$x, y = .data$region) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = loadtype, col = loadtype))+
  geom_segment(data=brms_high_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_high_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_high_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Region", title = "High impact SnpEff")+
  facet_wrap(~loadtype, ncol = 3, scales="free_x") +
  scale_fill_manual(values =alpha(c(clr_high,"#437689", "#5B98AE"), 0.7)) +
  scale_color_manual(values =c(clr_high,"#437689", "#5B98AE")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> supp_posterior_highregions

ggsave(supp_posterior_highregions, file = "plots/final/supp_posteriors_high_regions.png", width=12, height=16)


### Save together ####

cowplot::plot_grid(supp_posterior_gerpregions, supp_posterior_highregions, 
                   ncol = 1, 
                   align = "hv", axis = "lb", 
                   labels = "auto", label_fontface = "plain", label_size = 22) -> plotb_region

ggsave(plotb_region, file="plots/final/figure_posteriors_per_region.png", height = 12, width = 14)

##### Annual mating success ####
load(file = "results/bayes_models/annual_ms_froh_year_zi_29scaf.RData")
load(file = "results/bayes_models/annual_ms_gerp5_total_add_year_zi_29scaf.RData")
load(file = "results/bayes_models/annual_ms_snpef_high_total_add_year_zi_29scaf.RData")


#extract intervals and areas
brms_froh_ams_interval <- mcmc_intervals_data(brm_froh_annual_MS, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalefroh1", parameter))

brms_gerp5_ams_interval_total <- mcmc_intervals_data(brm_gerp5_total_annual_MS, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_high_ams_interval_total <- mcmc_intervals_data(brm_high_total_annual_MS, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_ams_interval <- rbind(brms_froh_ams_interval,
                           brms_gerp5_ams_interval_total,
                           brms_high_ams_interval_total)


brms_ams_interval$model <- c(rep("Inbreeding", nrow(brms_froh_ams_interval)),
                             rep("Total GERP load", nrow(brms_gerp5_ams_interval_total)),
                             rep("Total high impact load", nrow(brms_high_ams_interval_total)))


brms_ams_interval$model  <- factor(as.factor(brms_ams_interval$model),
                                   levels= c("Total high impact load", "Total GERP load", "Inbreeding"))

#area
brms_froh_ams_areas <- mcmc_areas_data(brm_froh_annual_MS) %>%
  subset(grepl("b_scalefroh1", parameter))

brms_gerp5_ams_areas_total <- mcmc_areas_data(brm_gerp5_total_annual_MS) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_high_ams_areas_total <- mcmc_areas_data(brm_high_total_annual_MS) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_ams_areas <- rbind(brms_froh_ams_areas,
                        brms_gerp5_ams_areas_total, 
                        brms_high_ams_areas_total)

brms_ams_areas$model <- c(rep("Inbreeding", nrow(brms_froh_ams_areas)),
                          rep("Total GERP load", nrow(brms_gerp5_ams_areas_total)),
                          rep("Total high impact load", nrow(brms_high_ams_areas_total)))


brms_ams_areas$model  <- factor(as.factor(brms_ams_areas$model),
                                levels= c("Total high impact load", "Total GERP load", "Inbreeding"))


### plot

# split by interval
brms_ams <- split(brms_ams_areas, brms_ams_areas$interval)

brms_ams$bottom <- brms_ams$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_ams$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_ams_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_ams_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_ams_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Genetic quality parameter")+
  scale_fill_manual(values =c(alpha(clr_high, 0.7), alpha(clr_gerp, 0.7), alpha(clr_froh, 0.7))) +
  scale_color_manual(values =c(clr_high, clr_gerp, clr_froh)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") -> supp_posteriors_ams
supp_posteriors_ams

ggsave(supp_posteriors_ams, file = "plots/final/supp_posterior_AMS.png", width = 10, height = 10)

## frequentist
## data
#load data
load("data/phenotypes/phenotypes_wide_extra.RData")
load("data/phenotypes/phenotypes_long_formodels_withsurv.RData") #annual
load("results/all_loads_combined_da_nosex_30scaf.RData") #loads no sex chr only 30 scaf

pheno <- pheno %>% mutate(core = as.factor(case_when(is.na(strategy) ~ "no_core",
                                                     TRUE ~ "core")))
core <- subset(pheno, core == "core")

# combine

pheno_load_annual <- left_join(pheno_long_models_ly, loads, by = "id")
pheno_load_annual$year <- as.factor(pheno_load_annual$year)
pheno_load_annual <- pheno_load_annual %>% mutate(core = as.factor(case_when(id %in% core$id ~ "core", TRUE ~ "no core")))
pheno_load_annual$born <- as.numeric(as.character(pheno_load_annual$year)) - pheno_load_annual$age

summary(glmmTMB::glmmTMB(MS ~ scale(froh1) + age + sqrt(age) + (1|site/id), 
                         data = subset(pheno_load_annual, core == "core"),
                         family = "poisson",  ziformula = ~1))

summary(glmmTMB::glmmTMB(MS ~ scale(gerp_Lt_additive_cat5) + age + sqrt(age) + (1|site/id), 
                         data = subset(pheno_load_annual, core == "core"),
                         family = poisson(link="log"), ziformula = ~1))

summary(glmmTMB::glmmTMB(MS ~ scale(high_Lt_additive) + age + sqrt(age) + (1|site/id), 
                         data = subset(pheno_load_annual, core == "core"),
                         family = "poisson", ziformula = ~1))

#makes a huge difference whether age is fixed or random

ggplot(pheno_load_annual, aes(x = froh1, y = MS)) + geom_point()
ggplot(pheno_load_annual, aes(x = gerp_Lt_additive_cat5, y = MS)) + geom_point()
ggplot(pheno_load_annual, aes(x = high_Lt_additive, y = MS)) + geom_point() + 
  geom_smooth(method="lm")
ggplot(pheno_load, aes(x = high_Lt_additive, y = LMS)) + geom_point()

##### Dominant models #####

#also add in additive for comparison

#load model
load(file = "results/bayes_models/lms_gerp5_total_codom_zi_29scaf.RData")
bayes_R2(brm_load_t_codom_gerp5_lms)
load(file = "results/bayes_models/lms_snpef_high_total_codom_zi_29scaf.RData")
bayes_R2(brm_load_t_codom_high_lms)
load(file = "results/bayes_models/lms_gerp5_total_add_1mil_zi_29scaf.RData")
brm_load_t_add_gerp5 <- brm_load_t_add_high_lms
load(file = "results/bayes_models/lms_snpef_high_total_add_1mil_zi_29scaf.RData")

#extract intervals and areas
brms_gerp5_lms_interval_total_a <- mcmc_intervals_data(brm_load_t_add_gerp5, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_high_lms_interval_total_a <- mcmc_intervals_data(brm_load_t_add_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_gerp5_lms_interval_total_c <- mcmc_intervals_data(brm_load_t_codom_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_cat5", parameter))

brms_high_lms_interval_total_c <- mcmc_intervals_data(brm_load_t_codom_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lt", parameter))

brms_total_codom_interval <- rbind(brms_gerp5_lms_interval_total_a,
                                   brms_high_lms_interval_total_a,
                                   brms_gerp5_lms_interval_total_c,
                                   brms_high_lms_interval_total_c)

brms_total_codom_interval$model <- c("GERP", "SnpEff", "GERP", "SnpEff")

brms_total_codom_interval$model  <- factor(as.factor(brms_total_codom_interval$model),
                                           levels= c("GERP", "SnpEff"))

brms_total_codom_interval$Dominance <- c("Additive", "Additive", "Dominant", "Dominant")

brms_total_codom_interval$Dominance  <- factor(as.factor(brms_total_codom_interval$Dominance),
                                               levels= c("Dominant", "Additive"))

write_tsv(brms_total_codom_interval, file="results/bayes_models/tables/codom_total_load_intervals.tsv")

#area
brms_gerp5_lms_areas_total_a <- mcmc_areas_data(brm_load_t_add_gerp5) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_high_lms_areas_total_a <- mcmc_areas_data(brm_load_t_add_high_lms) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_gerp5_lms_areas_total_c <- mcmc_areas_data(brm_load_t_codom_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lt_cat5", parameter))

brms_high_lms_areas_total_c <- mcmc_areas_data(brm_load_t_codom_high_lms) %>%
  subset(grepl("b_scalehigh_Lt", parameter))


brms_total_codom_areas <- rbind(brms_gerp5_lms_areas_total_a,
                                brms_high_lms_areas_total_a,
                                brms_gerp5_lms_areas_total_c,
                                brms_high_lms_areas_total_c)

brms_total_codom_areas$model <- c(rep("GERP", nrow(brms_gerp5_lms_areas_total_a)),
                                  rep("SnpEff", nrow(brms_high_lms_areas_total_a)),
                                  rep("GERP", nrow(brms_gerp5_lms_areas_total_c)),
                                  rep("SnpEff", nrow(brms_high_lms_areas_total_c)))


brms_total_codom_areas$model  <- factor(as.factor(brms_total_codom_areas$model),
                                        levels= c("GERP", "SnpEff"))

brms_total_codom_areas$Dominance <- c(rep("Additive", nrow(brms_gerp5_lms_areas_total_a)),
                                      rep("Additive", nrow(brms_high_lms_areas_total_a)),
                                      rep("Dominant", nrow(brms_gerp5_lms_areas_total_c)),
                                      rep("Dominant", nrow(brms_high_lms_areas_total_c)))


brms_total_codom_areas$Dominance  <- factor(as.factor(brms_total_codom_areas$Dominance),
                                            levels= c("Dominant", "Additive"))

### plot

# split by interval
brms_total_codom <- split(brms_total_codom_areas, brms_total_codom_areas$interval)

brms_total_codom$bottom <- brms_total_codom$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_total_codom$outer) +  
  aes(x = .data$x, y = .data$Dominance) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_total_codom_interval, aes(x = l, xend = h, yend = Dominance), col = "black", linewidth=3)+
  geom_segment(data=brms_total_codom_interval, aes(x = ll, xend = hh, yend = Dominance), col = "black")+
  geom_point(data=brms_total_codom_interval, aes(x = m, y = Dominance), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Method")+
  facet_wrap(~model, ncol=2, scales="free_x")+
  scale_fill_manual(values =alpha(c(clr_gerp, clr_high), 0.7)) +
  scale_color_manual(values =c(clr_gerp, clr_high)) +
  scale_x_continuous(n.breaks=3) +
  theme(plot.subtitle = element_text(hjust = .5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing=unit(2, "lines")) -> codom_plot
codom_plot
ggsave(codom_plot, file = "plots/final/supp_posteriors_total_codominant.png", width=10, height=10)
