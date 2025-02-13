### Packages ####
pacman::p_load(tidyverse, brms, bayesplot, data.table, ggridges)

# theme
source("scripts/theme_ggplot.R")

#### Approach 1: combining GERP and high mutations in one total load parameter ####
# combined
load(file = "output/models/total_hom_het/lms_total_gerp45_plus_high.RData")
brm_load_gerp_high <- brm_load_t

# gerp
load(file = "output/models/total_hom_het/lms_total_gerp45.RData")
brm_gerp <- brm_load_t

# snpeff
load(file = "output/models/total_hom_het/lms_total_high.RData")
brm_high <- brm_load_t

# get intervals
brms_both_lms_interval <- mcmc_intervals_data(brm_load_gerp_high, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
brms_gerp5_lms_interval <- mcmc_intervals_data(brm_gerp, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
brms_high_lms_interval <- mcmc_intervals_data(brm_high, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")

brms_plota_interval <- rbind(brms_both_lms_interval,
                             brms_gerp5_lms_interval,
                             brms_high_lms_interval)

brms_plota_interval$model <- c("Total GERP + SnpEff load", "Total GERP load", "Total SnpEff load")

# get areas
brms_both_lms_area <- mcmc_areas_data(brm_load_gerp_high, pars = "b_scaletotal_load")
brms_gerp5_lms_area <- mcmc_areas_data(brm_gerp, pars = "b_scaletotal_load")
brms_high_lms_area <- mcmc_areas_data(brm_high, pars = "b_scaletotal_load")

brms_plota_areas <- rbind(brms_both_lms_area,
                          brms_gerp5_lms_area,
                          brms_high_lms_area)

brms_plota_areas$model <- c(rep("Total GERP + SnpEff load", nrow(brms_both_lms_area)),
                            rep("Total GERP load", nrow(brms_gerp5_lms_area)),
                            rep("Total SnpEff load", nrow(brms_high_lms_area)))


#rearrange order for visualization
brms_plota_interval$model  <- factor(as.factor(brms_plota_interval$model),
                                     levels= c("Total SnpEff load", "Total GERP load","Total GERP + SnpEff load"))

brms_plota_areas$model  <- factor(as.factor(brms_plota_areas$model),
                                  levels= c("Total SnpEff load", "Total GERP load","Total GERP + SnpEff load"))

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
  xlim(-0.35, 0.1)+
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp, clr_highlight), 0.7)) +
  scale_color_manual(values =c(clr_high, clr_gerp, clr_highlight)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plot_posteriors_combined

plot_posteriors_combined

ggsave(plot_posteriors_combined, file = "plots/sup/combined_high_gerp_load.png", width=10, height=12)


#### Approach 2: both loads in the same model ####
load(file = "output/models/total_hom_het/lms_total_gerp45_high_sep.RData")
brm_load_gerp_high_sep <- brm_load_t

# gerp
load(file = "output/models/total_hom_het/lms_total_gerp45.RData")
brm_gerp <- brm_load_t

# snpeff
load(file = "output/models/total_hom_het/lms_total_high.RData")
brm_high <- brm_load_t

# get intervals
brms_gerp_both_lms_interval <- mcmc_intervals_data(brm_load_gerp_high_sep, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load_gerp45")
brms_high_both_lms_interval <- mcmc_intervals_data(brm_load_gerp_high_sep, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load_high")
brms_gerp5_lms_interval <- mcmc_intervals_data(brm_gerp, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
brms_high_lms_interval <- mcmc_intervals_data(brm_high, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")

brms_plota_interval <- rbind(brms_gerp_both_lms_interval,
                             brms_high_both_lms_interval,
                             brms_gerp5_lms_interval,
                             brms_high_lms_interval)

brms_plota_interval$model <- c("Combined model - GERP load", "Combined model - SnpEff load", "Total GERP load", "Total SnpEff load")

# get areas
brms_gerp_both_lms_area <- mcmc_areas_data(brm_load_gerp_high_sep, pars = "b_scaletotal_load_gerp45")
brms_high_both_lms_area <- mcmc_areas_data(brm_load_gerp_high_sep, pars = "b_scaletotal_load_high")
brms_gerp5_lms_area <- mcmc_areas_data(brm_gerp, pars = "b_scaletotal_load")
brms_high_lms_area <- mcmc_areas_data(brm_high, pars = "b_scaletotal_load")

brms_plota_areas <- rbind(brms_gerp_both_lms_area,
                          brms_high_both_lms_area,
                          brms_gerp5_lms_area,
                          brms_high_lms_area)

brms_plota_areas$model <- c(rep("Combined model - GERP load", nrow(brms_gerp_both_lms_area)),
                            rep("Combined model - SnpEff load", nrow(brms_high_both_lms_area)),
                            rep("Total GERP load", nrow(brms_gerp5_lms_area)),
                            rep("Total SnpEff load", nrow(brms_high_lms_area)))


#rearrange order for visualization
brms_plota_interval$model  <- factor(as.factor(brms_plota_interval$model),
                                     levels= c("Total SnpEff load","Total GERP load","Combined model - SnpEff load", "Combined model - GERP load"))

brms_plota_areas$model  <- factor(as.factor(brms_plota_areas$model),
                                  levels= c("Total SnpEff load","Total GERP load","Combined model - SnpEff load", "Combined model - GERP load"))

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
  xlim(-0.35, 0.1)+
  scale_fill_manual(values =alpha(c(clr_gerp, clr_high, clr_gerp, clr_high), 0.7)) +
  scale_color_manual(values =c( clr_gerp, clr_high, clr_gerp, clr_high)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plot_posteriors_one_model

plot_posteriors_one_model

ggsave(plot_posteriors_one_model, file = "plots/sup/separate_one_model_high_gerp_load.png", width=10, height=12)

#### Combine in one figure ####
cowplot::plot_grid(plot_posteriors_one_model, plot_posteriors_combined,
                   ncol = 1, labels = c("a", "b"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> supp_different_modelling

ggsave(file = "plots/sup/sup_different_modelling.png", width=10, height=12)
