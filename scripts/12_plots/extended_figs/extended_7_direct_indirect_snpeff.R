#### Plotting posterior distributions for plot with diagram and
#### grouse illustrations added using figma

#packages
pacman::p_load(brms, bayesplot, tidyverse, ggridges, performance)
source("scripts/theme_ggplot.R")

#### SnpEff - supp ####
##### Traits on MS #####
# load data
load(file = "output/models/annual/ams/model_trait_ams_high.RData")
fit_ms_high <- fit
r2_bayes(fit_ms_high)

#extract intervals and areas
brms_trait_ms_high_interval <- mcmc_intervals_data(fit_ms_high, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scale", parameter))
brms_trait_ms_high_interval <- subset(brms_trait_ms_high_interval, parameter != "b_scaletotal_load")

brms_trait_ms_high_interval$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_high_interval$parameter)
brms_trait_ms_high_interval$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_high_interval$parameter)
brms_trait_ms_high_interval$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_high_interval$parameter)
brms_trait_ms_high_interval$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_high_interval$parameter)
brms_trait_ms_high_interval$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_high_interval$parameter)
brms_trait_ms_high_interval$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_high_interval$parameter)

brms_trait_ms_high_interval$parameter <- factor(brms_trait_ms_high_interval$parameter, 
                                                levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                           "Centrality", "Fighting", "Attendance"))


#area
brms_trait_ms_high_areas <- mcmc_areas_data(fit_ms_high) %>%
  subset(grepl("b_scale", parameter))
brms_trait_ms_high_areas <- subset(brms_trait_ms_high_areas, parameter != "b_scaletotal_load")

brms_trait_ms_high_areas$parameter <- gsub("b_scaleblue", "Blue chroma", brms_trait_ms_high_areas$parameter)
brms_trait_ms_high_areas$parameter <- gsub("b_scaleeyec", "Eye comb", brms_trait_ms_high_areas$parameter)
brms_trait_ms_high_areas$parameter <- gsub("b_scalelyre", "Lyre size", brms_trait_ms_high_areas$parameter)
brms_trait_ms_high_areas$parameter <- gsub("b_scaledist", "Centrality", brms_trait_ms_high_areas$parameter)
brms_trait_ms_high_areas$parameter <- gsub("b_scalefight", "Fighting", brms_trait_ms_high_areas$parameter)
brms_trait_ms_high_areas$parameter <- gsub("b_scaleattend", "Attendance", brms_trait_ms_high_areas$parameter)

brms_trait_ms_high_areas$parameter <- factor(brms_trait_ms_high_areas$parameter, 
                                             levels = c("Blue chroma", "Eye comb", "Lyre size",
                                                        "Centrality", "Fighting", "Attendance"))

### plot

# split by interval
brms_trait_ms_high <- split(brms_trait_ms_high_areas, brms_trait_ms_high_areas$interval)

brms_trait_ms_high$bottom <- brms_trait_ms_high$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_trait_ms_high$outer) +  
  aes(x = .data$x, y = .data$parameter) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = parameter, col = parameter))+
  geom_segment(data=brms_trait_ms_high_interval, aes(x = l, xend = h, yend = parameter), col = "black", linewidth=3)+
  geom_segment(data=brms_trait_ms_high_interval, aes(x = ll, xend = hh, yend = parameter), col = "black")+
  geom_point(data=brms_trait_ms_high_interval, aes(x = m, y = parameter), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 1.5)+
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1],clrs_hunting[1],
                                    clrs_hunting[1],clrs_hunting[1],
                                    clrs_hunting[1],clrs_hunting[1]), 0.5)) +
  scale_color_manual(values =c(clrs_hunting[1],clrs_hunting[1],
                               clrs_hunting[1],clrs_hunting[1],
                               clrs_hunting[1],clrs_hunting[1])) +
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
  )-> traits_ms_high_posterior
traits_ms_high_posterior

png(file = "plots/extended/extended_7_high_right_traits_ams.png", width=600, height=800, bg='transparent')
traits_ms_high_posterior
dev.off()

write.csv(brms_trait_ms_high_interval, file = "output/models/intervals/high_traits_ms.csv", quote=F, row.names = F)

##### SnpEff load on traits ####

# load data
load(file = "output/models/annual/traits/model_attend_high.RData")
fit_high_attend <- fit
r2_bayes(fit_high_attend)

load(file = "output/models/annual/traits/model_fight_high.RData")
fit_high_fight <- fit
r2_bayes(fit_high_fight)

load(file = "output/models/annual/traits/model_dist_high.RData")
fit_high_dist <- fit
r2_bayes(fit_high_dist)

load(file = "output/models/annual/traits/model_eyec_high.RData")
fit_high_eyec <- fit
r2_bayes(fit_high_eyec)

load(file = "output/models/annual/traits/model_blue_high.RData")
fit_high_blue <- fit
r2_bayes(fit_high_blue)

load(file = "output/models/annual/traits/model_lyre_high.RData")
fit_high_lyre <- fit
r2_bayes(fit_high_lyre)
rm(fit)

#extract intervals and areas
attend_interval <- mcmc_intervals_data(fit_high_attend, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

fight_interval <- mcmc_intervals_data(fit_high_fight, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

dist_interval <- mcmc_intervals_data(fit_high_dist, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

eyec_interval <- mcmc_intervals_data(fit_high_eyec, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

blue_interval <- mcmc_intervals_data(fit_high_blue, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

lyre_interval <- mcmc_intervals_data(fit_high_lyre, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

hightrait_interval <- rbind(attend_interval,
                            fight_interval,
                            dist_interval, 
                            eyec_interval, 
                            blue_interval, 
                            lyre_interval)


hightrait_interval$trait <- c("Attendance", "Fighting", "Centrality", "Eye comb", "Blue chroma", "Lyre size")

hightrait_interval$trait <- factor(hightrait_interval$trait, 
                                   levels = c("Blue chroma", "Eye comb", "Lyre size",
                                              "Centrality", "Fighting", "Attendance"))

#area
attend_area <- mcmc_areas_data(fit_high_attend) %>%
  subset(grepl("scaletotal_load", parameter))

fight_area <- mcmc_areas_data(fit_high_fight) %>%
  subset(grepl("scaletotal_load", parameter))

dist_area <- mcmc_areas_data(fit_high_dist) %>%
  subset(grepl("scaletotal_load", parameter))

eyec_area <- mcmc_areas_data(fit_high_eyec) %>%
  subset(grepl("scaletotal_load", parameter))

blue_area <- mcmc_areas_data(fit_high_blue) %>%
  subset(grepl("scaletotal_load", parameter))

lyre_area <- mcmc_areas_data(fit_high_lyre) %>%
  subset(grepl("scaletotal_load", parameter))

hightrait_area <- rbind(attend_area,
                        fight_area, 
                        dist_area,
                        eyec_area,
                        blue_area,
                        lyre_area)

hightrait_area$trait <- c(rep("Attendance", nrow(attend_area)),
                          rep("Fighting", nrow(fight_area)),
                          rep("Centrality", nrow(dist_area)),
                          rep("Eye comb", nrow(eyec_area)),
                          rep("Blue chroma", nrow(blue_area)),
                          rep("Lyre size", nrow(lyre_area)))

hightrait_area$trait <- factor(hightrait_area$trait, 
                               levels = c("Blue chroma", "Eye comb", "Lyre size",
                                          "Centrality", "Fighting", "Attendance"))

### plot

# split by interval
hightrait <- split(hightrait_area, hightrait_area$interval)

hightrait$bottom <- hightrait$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = hightrait$outer) +  
  aes(x = .data$x, y = .data$trait) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = trait, col = trait))+
  geom_segment(data=hightrait_interval, aes(x = l, xend = h, yend = trait), col = "black", linewidth=3)+
  geom_segment(data=hightrait_interval, aes(x = ll, xend = hh, yend = trait), col = "black")+
  geom_point(data=hightrait_interval, aes(x = m, y = trait), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", size = 1.5)+
  labs(x = expression("Standardised"~beta), y = "Trait")+
  scale_fill_manual(values =alpha(c(clrs_hunting[1],clrs_hunting[1],
                                    clrs_hunting[1],clrs_hunting[1],
                                    clrs_hunting[1],clrs_hunting[1]), 0.5)) +
  scale_color_manual(values =c(clrs_hunting[1],clrs_hunting[1],
                               clrs_hunting[1],clrs_hunting[1],
                               clrs_hunting[1],clrs_hunting[1])) +
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
  )-> traits_high_posterior

traits_high_posterior


png(file = "plots/extended/extended_7_high_left_load_traits.png", width=600, height=800, bg='transparent')
traits_high_posterior
dev.off()

write.csv(hightrait_interval, file = "output/models/intervals/high_load_traits.csv", quote=F, row.names = F)
