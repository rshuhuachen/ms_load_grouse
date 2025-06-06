#### Plotting posterior distributions for plot with diagram and
#### grouse illustrations added using figma

#packages
pacman::p_load(brms, bayesplot, tidyverse, ggridges, performance)
source("scripts/theme_ggplot.R")

#### GERP - main ####
##### Traits on MS #####

# load data
load(file = "output/models/annual/ams/model_trait_ams_gerp45.RData")
fit_ms_gerp <- fit

interval_full_gerp <- mcmc_intervals_data(fit_ms_gerp, prob =0.8, prob_outer = 0.95)
interval_full_gerp <- data.frame(parameter = interval_full_gerp$parameter,
                                        median = round(interval_full_gerp$m, 2),
                                        ci_95 = paste0(round(interval_full_gerp$ll, 2), ", ", round(interval_full_gerp$hh, 2)),
                                        ci_80 = paste0(round(interval_full_gerp$l, 2), ", ", round(interval_full_gerp$h, 2)))

write_tsv(interval_full_gerp, file = "output/models/intervals/ams_gerp45_full.tsv")


#extract intervals and areas
brms_trait_ms_gerp_interval <- mcmc_intervals_data(fit_ms_gerp, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scale", parameter))
brms_trait_ms_gerp_interval <- subset(brms_trait_ms_gerp_interval, parameter != "b_scaletotal_load")

brms_trait_ms_gerp_interval$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_gerp_interval$parameter)
brms_trait_ms_gerp_interval$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_gerp_interval$parameter)
brms_trait_ms_gerp_interval$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_gerp_interval$parameter)
brms_trait_ms_gerp_interval$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_gerp_interval$parameter)
brms_trait_ms_gerp_interval$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_gerp_interval$parameter)
brms_trait_ms_gerp_interval$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_gerp_interval$parameter)

brms_trait_ms_gerp_interval$parameter <- factor(brms_trait_ms_gerp_interval$parameter, 
                                        levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                   "Centrality", "Fighting", "Attendance"))


#area
brms_trait_ms_gerp_areas <- mcmc_areas_data(fit_ms_gerp) %>%
  subset(grepl("b_scale", parameter))
brms_trait_ms_gerp_areas <- subset(brms_trait_ms_gerp_areas, parameter != "b_scaletotal_load")

brms_trait_ms_gerp_areas$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_gerp_areas$parameter)
brms_trait_ms_gerp_areas$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_gerp_areas$parameter)
brms_trait_ms_gerp_areas$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_gerp_areas$parameter)
brms_trait_ms_gerp_areas$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_gerp_areas$parameter)
brms_trait_ms_gerp_areas$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_gerp_areas$parameter)
brms_trait_ms_gerp_areas$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_gerp_areas$parameter)

brms_trait_ms_gerp_areas$parameter <- factor(brms_trait_ms_gerp_areas$parameter, 
                                        levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                   "Centrality", "Fighting", "Attendance"))

### plot

# split by interval
brms_trait_ms_gerp <- split(brms_trait_ms_gerp_areas, brms_trait_ms_gerp_areas$interval)

brms_trait_ms_gerp$bottom <- brms_trait_ms_gerp$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_trait_ms_gerp$outer) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=brms_trait_ms_gerp_interval, aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=brms_trait_ms_gerp_interval, aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=brms_trait_ms_gerp_interval, aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 1.5)+
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clr_grey,clr_grey,
                                    clr_grey,clr_grey,
                                    clr_grey,clr_highlight), 0.5)) +
  scale_color_manual(values =c(clr_grey,clr_grey,
                               clr_grey,clr_grey,
                               clr_grey,clr_highlight)) +
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
  )-> traits_ms_gerp_posterior
traits_ms_gerp_posterior

png(file = "plots/main/fig_4_right_traits_ams.png", width=600, height=800, bg='transparent')
traits_ms_gerp_posterior
dev.off()

brms_trait_ms_gerp_interval <- brms_trait_ms_gerp_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_trait_ms_gerp_interval, file = "output/models/intervals/gerp_traits_ms.csv", quote=F, row.names = F)

##### GERP load on traits ####

# load data
load(file = "output/models/annual/traits/model_attend_gerp45.RData")
fit_gerp_attend <- fit
r2_bayes(fit_gerp_attend)
mcmc_intervals_data(fit_ms_gerp, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_attend.tsv")

load(file = "output/models/annual/traits/model_fight_gerp45.RData")
fit_gerp_fight <- fit
r2_bayes(fit_gerp_fight)
mcmc_intervals_data(fit_gerp_fight, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_fight.tsv")

load(file = "output/models/annual/traits/model_dist_gerp45.RData")
fit_gerp_dist <- fit
r2_bayes(fit_gerp_dist)
mcmc_intervals_data(fit_gerp_dist, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_dist.tsv")

load(file = "output/models/annual/traits/model_eyec_gerp45.RData")
fit_gerp_eyec <- fit
r2_bayes(fit_gerp_eyec)
mcmc_intervals_data(fit_gerp_eyec, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_eyec.tsv")

load(file = "output/models/annual/traits/model_blue_gerp45.RData")
fit_gerp_blue <- fit
r2_bayes(fit_gerp_blue)
mcmc_intervals_data(fit_gerp_blue, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_blue.tsv")

load(file = "output/models/annual/traits/model_lyre_gerp45.RData")
fit_gerp_lyre <- fit
r2_bayes(fit_gerp_lyre)
mcmc_intervals_data(fit_gerp_lyre, prob =0.8, prob_outer = 0.95) %>%
  write_tsv(file = "output/models/intervals/fit_gerp_lyre.tsv")

rm(fit)

#extract intervals and areas
attend_interval <- mcmc_intervals_data(fit_gerp_attend, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

fight_interval <- mcmc_intervals_data(fit_gerp_fight, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

dist_interval <- mcmc_intervals_data(fit_gerp_dist, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

eyec_interval <- mcmc_intervals_data(fit_gerp_eyec, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

blue_interval <- mcmc_intervals_data(fit_gerp_blue, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

lyre_interval <- mcmc_intervals_data(fit_gerp_lyre, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

gerptrait_interval <- rbind(attend_interval,
                            fight_interval,
                            dist_interval, 
                            eyec_interval, 
                            blue_interval, 
                            lyre_interval)


gerptrait_interval$trait <- c("Attendance", "Fighting", "Centrality", "Eye comb", "Blue chroma", "Lyre size")

gerptrait_interval$trait <- factor(gerptrait_interval$trait, 
                                             levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                        "Centrality", "Fighting", "Attendance"))

#area
attend_area <- mcmc_areas_data(fit_gerp_attend) %>%
  subset(grepl("scaletotal_load", parameter))

fight_area <- mcmc_areas_data(fit_gerp_fight) %>%
  subset(grepl("scaletotal_load", parameter))

dist_area <- mcmc_areas_data(fit_gerp_dist) %>%
  subset(grepl("scaletotal_load", parameter))

eyec_area <- mcmc_areas_data(fit_gerp_eyec) %>%
  subset(grepl("scaletotal_load", parameter))

blue_area <- mcmc_areas_data(fit_gerp_blue) %>%
  subset(grepl("scaletotal_load", parameter))

lyre_area <- mcmc_areas_data(fit_gerp_lyre) %>%
  subset(grepl("scaletotal_load", parameter))

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
                                       levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                  "Centrality", "Fighting", "Attendance"))

### plot

# split by interval
gerptrait <- split(gerptrait_area, gerptrait_area$interval)

gerptrait$bottom <- gerptrait$outer %>%
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
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clr_grey,clr_grey,
                                    clr_grey,clr_grey,
                                    clr_grey,clr_highlight), 0.5)) +
  scale_color_manual(values =c(clr_grey,clr_grey,
                               clr_grey,clr_grey,
                               clr_grey,clr_highlight)) +
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

traits_gerp_posterior


png(file = "plots/main/fig_4_left_load_traits.png", width=600, height=800, bg='transparent')
traits_gerp_posterior
dev.off()

write.csv(gerptrait_interval, file = "output/models/intervals/gerp_load_traits.csv", quote=F, row.names = F)


#### for NEE ###

small_font = 10
large_font = 12
in_box_font = 3.5

theme_set(theme_classic() + theme(title = element_text(small_font),
                                  plot.subtitle = element_text(size=large_font),
                                  axis.title = element_text(size = large_font, family = "Arial"),
                                  axis.text = element_text(size = small_font, family = "Arial"),
                                  text=element_text(size=small_font, family = "Arial"),
                                  legend.text =  element_text(size = small_font, family = "Arial"),
                                  legend.title = element_text(size = small_font, family = "Arial"),
                                  strip.text = element_text(size = small_font, family = "Arial"),
                                  axis.title.y = element_text(margin = margin(t = 0, r =5, b = 0, l = 0),
                                                              color = "black"),
                                  plot.margin = margin(0.3,0.3,0.3,0.3, "cm"),
                                  plot.title=element_text(margin=margin(0,0,5,0)),
                                  axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0),
                                                              color = "black"),
                                  panel.background = element_rect(fill = "white", colour = NA),
                                  plot.background = element_rect(fill = "white", colour = NA),))

ggplot(data = brms_trait_ms_gerp$outer) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=brms_trait_ms_gerp_interval, aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=brms_trait_ms_gerp_interval, aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=brms_trait_ms_gerp_interval, aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 3) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 0.5)+
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clr_grey,clr_grey,
                                    clr_grey,clr_grey,
                                    clr_grey,clr_highlight), 0.5)) +
  scale_color_manual(values =c(clr_grey,clr_grey,
                               clr_grey,clr_grey,
                               clr_grey,clr_highlight)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
 )-> traits_ms_gerp_posterior_nee

traits_ms_gerp_posterior_nee

ggsave(plot = traits_ms_gerp_posterior_nee, filename = 'plots/main/fig_4_right_traits_ams_nee.png', width = 90, height = 180,
       bg='transparent', units = 'mm')

ggplot(data = gerptrait$outer) +  
  aes(x = .data$x, y = .data$trait) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = trait, col = trait))+
  geom_segment(data=gerptrait_interval, aes(x = l, xend = h, yend = trait), col = "black", linewidth=3)+
  geom_segment(data=gerptrait_interval, aes(x = ll, xend = hh, yend = trait), col = "black")+
  geom_point(data=gerptrait_interval, aes(x = m, y = trait), fill="white",  col = "black", shape=21, size = 3) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 0.5)+
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clr_grey,clr_grey,
                                    clr_grey,clr_grey,
                                    clr_grey,clr_highlight), 0.5)) +
  scale_color_manual(values =c(clr_grey,clr_grey,
                               clr_grey,clr_grey,
                               clr_grey,clr_highlight)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )-> traits_gerp_posterior_nee

ggsave(plot = traits_ms_gerp_posterior_nee, filename = 'plots/main/fig_4_left_load_traits_nee.png', width = 90, height = 180,
       bg='transparent',
       units = 'mm')



