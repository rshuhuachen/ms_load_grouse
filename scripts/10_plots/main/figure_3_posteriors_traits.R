#### Plotting posterior distributions for plot with diagram and
#### grouse illustrations added using figma

#packages
pacman::p_load(brms, bayesplot, tidyverse)
source("scripts/theme_ggplot.R")

#### Traits on MS ####

# load data
load(file = "output/models/traits/ms_alltraits.RData")

#extract intervals and areas
brms_trait_ms_interval <- mcmc_intervals_data(fit_MS, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scale", parameter))

brms_trait_ms_interval$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_interval$parameter)
brms_trait_ms_interval$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_interval$parameter)
brms_trait_ms_interval$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_interval$parameter)
brms_trait_ms_interval$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_interval$parameter)
brms_trait_ms_interval$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_interval$parameter)
brms_trait_ms_interval$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_interval$parameter)

brms_trait_ms_interval$parameter <- factor(brms_trait_ms_interval$parameter, 
                                           levels = c("Centrality", "Fighting", "Attendance",
                                                                    "Blue chroma", "Eye comb", "Lyre size"))

#area
brms_trait_ms_areas <- mcmc_areas_data(fit_MS) %>%
  subset(grepl("b_scale", parameter))

brms_trait_ms_areas$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_areas$parameter)
brms_trait_ms_areas$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_areas$parameter)
brms_trait_ms_areas$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_areas$parameter)
brms_trait_ms_areas$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_areas$parameter)
brms_trait_ms_areas$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_areas$parameter)
brms_trait_ms_areas$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_areas$parameter)

brms_trait_ms_areas$parameter <- factor(brms_trait_ms_areas$parameter, 
                                           levels = c("Centrality", "Fighting", "Attendance",
                                                      "Blue chroma", "Eye comb", "Lyre size"))

### plot

# split by interval
brms_trait_ms <- split(brms_trait_ms_areas, brms_trait_ms_areas$interval)

brms_trait_ms$bottom <- brms_trait_ms$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_trait_ms$outer) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=brms_trait_ms_interval, aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=brms_trait_ms_interval, aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=brms_trait_ms_interval, aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 1.5)+
  labs(x = "Standardised beta coefficient", y = "Trait")+
  scale_fill_manual(values =alpha(c(clrs_related), 0.7)) +
  scale_color_manual(values =c(clrs_related)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )-> traits_ms_posterior


ggsave(traits_ms_posterior, file = "plots/Figure_3_traits_AMS", height = 12, width=8)

#### GERP load on traits ####

# load data
load(file = "output/models/traits/model_attend_gerp5_load_tadd_all.RData")
fit_gerp_attend <- fit
load(file = "output/models/traits/model_fight_gerp5_load_tadd_all.RData")
fit_gerp_fight <- fit
load(file = "output/models/traits/model_dist_gerp5_load_tadd_all.RData")
fit_gerp_dist <- fit
load(file = "output/models/traits/model_eyec_gerp5_load_tadd_all.RData")
fit_gerp_eyec <- fit
load(file = "output/models/traits/model_blue_gerp5_load_tadd_all.RData")
fit_gerp_blue <- fit
load(file = "output/models/traits/model_lyre_gerp5_load_tadd_all.RData")
fit_gerp_lyre <- fit
rm(fit)

#extract intervals and areas
attend_interval <- mcmc_intervals_data(fit_gerp_attend, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

fight_interval <- mcmc_intervals_data(fit_gerp_fight, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

dist_interval <- mcmc_intervals_data(fit_gerp_dist, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

eyec_interval <- mcmc_intervals_data(fit_gerp_eyec, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

blue_interval <- mcmc_intervals_data(fit_gerp_blue, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

lyre_interval <- mcmc_intervals_data(fit_gerp_lyre, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

gerptrait_interval <- rbind(attend_interval,
                            fight_interval,
                            dist_interval, 
                            eyec_interval, 
                            blue_interval, 
                            lyre_interval)


gerptrait_interval$trait <- c("Attendance", "Fighting", "Centrality", "Eye comb", "Blue chroma", "Lyre size")

gerptrait_interval$trait <- factor(gerptrait_interval$trait, 
                                           levels = c("Centrality", "Fighting", "Attendance",
                                                      "Blue chroma", "Eye comb", "Lyre size"))

#area
attend_area <- mcmc_areas_data(fit_gerp_attend) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

fight_area <- mcmc_areas_data(fit_gerp_fight) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

dist_area <- mcmc_areas_data(fit_gerp_dist) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

eyec_area <- mcmc_areas_data(fit_gerp_eyec) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

blue_area <- mcmc_areas_data(fit_gerp_blue) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

lyre_area <- mcmc_areas_data(fit_gerp_lyre) %>%
  subset(grepl("b_scalegerp5_load_tadd_all", parameter))

gerptrait_area <- rbind(attend_area,
                        fight_area, 
                        dist_area,
                        eyec_area,
                        blue_area,
                        lyre_area)

gerptrait_area$trait <- c(rep("Attendance", nrow(attend_area)),
                          rep("Fighting", nrow(fight_area)),
                          rep("Centrality", nrow(dist_area)),
                          rep("Eye comb", nrow(eyec_area)),
                          rep("Blue chroma", nrow(blue_area)),
                          rep("Lyre size", nrow(lyre_area)))

gerptrait_area$trait <- factor(gerptrait_area$trait, 
                                        levels = c("Centrality", "Fighting", "Attendance",
                                                   "Blue chroma", "Eye comb", "Lyre size"))

### plot

# split by interval
gerptrait <- split(gerptrait_area, gerptrait_area$interval)

gerptrait$bottom <- gerptrait$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = gerptrait$outer) +  
  aes(x = .data$x, y = .data$trait) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = trait, col = trait))+
  geom_segment(data=gerptrait_interval, aes(x = l, xend = h, yend = trait), col = "black", linewidth=3)+
  geom_segment(data=gerptrait_interval, aes(x = ll, xend = hh, yend = trait), col = "black")+
  geom_point(data=gerptrait_interval, aes(x = m, y = trait), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 1.5)+
  labs(x = "Standardised beta coefficient", y = "Trait")+
  scale_fill_manual(values =alpha(c(rep(clr_gerp, 6)), 0.7)) +
  scale_color_manual(values =c(rep(clr_gerp, 6))) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )-> traits_gerp_posterior


ggsave(traits_gerp_posterior, file = "plots/Figure_3_gerp_traits.png", height = 12, width=8)

### save intervals ####
write_tsv(gerptrait_interval, file="output/models/intervals/traits_gerp_intervals.tsv")
write_tsv(brms_trait_ms_interval, file="output/models/intervals/ms_traits_intervals.tsv")
