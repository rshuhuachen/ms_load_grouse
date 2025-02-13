###### This script contains the plot that contains the posterior distributions and raw data ####

### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, 
               extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

#### plot a - froh, total gerp, total snpeff ####

# froh
load(file = "output/models/from_lms.RData")
brm_froh_lms <- fit
r2_bayes(brm_froh_lms)

# load models

# gerp
load(file = "output/models/total_hom_het/lms_total_gerp45.RData")
brm_load_t_gerp5_lms <- brm_load_t
r2_bayes(brm_load_t_gerp5_lms)

interval_full_gerp <- mcmc_intervals_data(brm_load_t_gerp5_lms, prob =0.8, prob_outer = 0.95)
interval_full_gerp <- data.frame(parameter = interval_full_gerp$parameter,
                                 median = round(interval_full_gerp$m, 2),
                                 ci_95 = paste0(round(interval_full_gerp$ll, 2), ", ", round(interval_full_gerp$hh, 2)),
                                 ci_80 = paste0(round(interval_full_gerp$l, 2), ", ", round(interval_full_gerp$h, 2)))

write_tsv(interval_full_gerp, file = "output/models/intervals/total_gerp45_full.tsv")

# snpeff
load(file = "output/models/total_hom_het/lms_total_high.RData")
brm_load_t_add_high_lms <- brm_load_t
r2_bayes(brm_load_t_add_high_lms)

interval_full_snpeff <- mcmc_intervals_data(brm_load_t_add_high_lms, prob =0.8, prob_outer = 0.95)
interval_full_snpeff <- data.frame(parameter = interval_full_snpeff$parameter,
                                 median = round(interval_full_snpeff$m, 2),
                                 ci_95 = paste0(round(interval_full_snpeff$ll, 2), ", ", round(interval_full_snpeff$hh, 2)),
                                 ci_80 = paste0(round(interval_full_snpeff$l, 2), ", ", round(interval_full_snpeff$h, 2)))

write_tsv(interval_full_snpeff, file = "output/models/intervals/total_snpeff_full.tsv")

# get intervals
brms_froh_lms_interval <- mcmc_intervals_data(brm_froh_lms, prob =0.8, prob_outer = 0.95, pars = "b_scalefroh")
brms_gerp5_lms_interval <- mcmc_intervals_data(brm_load_t_gerp5_lms, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
brms_high_lms_interval <- mcmc_intervals_data(brm_load_t_add_high_lms, prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")

brms_plota_interval <- rbind(brms_froh_lms_interval,
                             brms_gerp5_lms_interval,
                             brms_high_lms_interval)
brms_plota_interval$model <- c("FROH", "Total GERP load", "Total SnpEff load")

# export full intervals for gerp load

# get areas
brms_froh_lms_area <- mcmc_areas_data(brm_froh_lms, pars = "b_scalefroh")
brms_gerp5_lms_area <- mcmc_areas_data(brm_load_t_gerp5_lms, pars = "b_scaletotal_load")
brms_high_lms_area <- mcmc_areas_data(brm_load_t_add_high_lms, pars = "b_scaletotal_load")

brms_plota_areas <- rbind(brms_froh_lms_area,
                          brms_gerp5_lms_area,
                          brms_high_lms_area)

brms_plota_areas$model <- c(rep("FROH", nrow(brms_froh_lms_area)),
                            rep("Total GERP load", nrow(brms_gerp5_lms_area)),
                            rep("Total SnpEff load", nrow(brms_high_lms_area)))


#rearrange order for visualization
brms_plota_interval$model  <- factor(as.factor(brms_plota_interval$model),
                                     levels= c("Total SnpEff load", "Total GERP load","FROH"))

brms_plota_areas$model  <- factor(as.factor(brms_plota_areas$model),
                                  levels= c("Total SnpEff load", "Total GERP load","FROH"))

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
  scale_fill_manual(values =alpha(c(clr_high, clr_gerp, clr_froh), 0.7)) +
  scale_y_discrete(labels = c("Total SnpEff load", "Total GERP load", 
                              expression(italic(F)[ROH])))+
  scale_color_manual(values =c(clr_high, clr_gerp, clr_froh)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> totals

totals

brms_plota_interval <- brms_plota_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_plota_interval, file = "output/models/intervals/total_froh_gerp45_high.csv", quote=F, row.names = F)

png(file = "plots/main/fig_2a.png", width=600, height=600)
totals
dev.off()

##### plot b - hom and het  #####

load(file = "output/models/total_hom_het/lms_het_hom_gerp45.RData")
brm_load_rp_gerp5_lms <- brm_load_het_hom
r2_bayes(brm_load_rp_gerp5_lms)

interval_full_homhet_gerp <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95)
interval_full_homhet_gerp <- data.frame(parameter = interval_full_homhet_gerp$parameter,
                                        median = round(interval_full_homhet_gerp$m, 2),
                                        ci_95 = paste0(round(interval_full_homhet_gerp$ll, 2), ", ", round(interval_full_homhet_gerp$hh, 2)),
                                        ci_80 = paste0(round(interval_full_homhet_gerp$l, 2), ", ", round(interval_full_homhet_gerp$h, 2)))

write_tsv(interval_full_homhet_gerp, file = "output/models/intervals/homhet_gerp45_full.tsv")


load(file = "output/models/total_hom_het/lms_het_hom_high.RData")
brm_load_rp_high_lms <- brm_load_het_hom
r2_bayes(brm_load_rp_high_lms)

#extract intervals and areas
brms_gerp5_lms_interval_hom <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehom_load", parameter))

brms_gerp5_lms_interval_het <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehet_load", parameter))

brms_high_lms_interval_hom <- mcmc_intervals_data(brm_load_rp_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehom_load", parameter))

brms_high_lms_interval_het <- mcmc_intervals_data(brm_load_rp_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehet_load", parameter))

brms_hom_het_lms_interval <- rbind(brms_gerp5_lms_interval_hom,
                                   brms_gerp5_lms_interval_het,
                                   brms_high_lms_interval_hom,
                                   brms_high_lms_interval_het)

brms_hom_het_lms_interval$model <- c("Hom", "Het","Hom", "Het")

brms_hom_het_lms_interval$model  <- factor(as.factor(brms_hom_het_lms_interval$model),
                                           levels= c("Het", "Hom"))

brms_hom_het_lms_interval$loadtype <- c("GERP", "GERP","SnpEff", "SnpEff")

brms_hom_het_lms_interval$loadtype  <- factor(as.factor(brms_hom_het_lms_interval$loadtype),
                                              levels= c("GERP", "SnpEff"))

#area
brms_gerp5_lms_areas_hom <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalehom_load", parameter))

brms_gerp5_lms_areas_het <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalehet_load", parameter))

brms_high_lms_areas_hom <- mcmc_areas_data(brm_load_rp_high_lms) %>%
  subset(grepl("b_scalehom_load", parameter))

brms_high_lms_areas_het <- mcmc_areas_data(brm_load_rp_high_lms) %>%
  subset(grepl("b_scalehet_load", parameter))


brms_hom_het_lms_areas <- rbind(brms_gerp5_lms_areas_hom, 
                                brms_gerp5_lms_areas_het,
                                brms_high_lms_areas_hom,
                                brms_high_lms_areas_het)

brms_hom_het_lms_areas$model <- c(rep("Hom", nrow(brms_gerp5_lms_areas_hom)),
                                  rep("Het", nrow(brms_gerp5_lms_areas_het)),
                                  rep("Hom", nrow(brms_high_lms_areas_hom)),
                                  rep("Het", nrow(brms_high_lms_areas_het)))

brms_hom_het_lms_areas$model  <- factor(as.factor(brms_hom_het_lms_areas$model),
                                        levels= c("Het", "Hom"))

brms_hom_het_lms_areas$loadtype <- c(rep("GERP", nrow(brms_gerp5_lms_areas_hom)),
                                     rep("GERP", nrow(brms_gerp5_lms_areas_het)),
                                     rep("SnpEff", nrow(brms_high_lms_areas_hom)),
                                     rep("SnpEff", nrow(brms_high_lms_areas_het)))

brms_hom_het_lms_areas$loadtype  <- factor(as.factor(brms_hom_het_lms_areas$loadtype),
                                           levels= c("GERP", "SnpEff"))

### plot

# split by interval
brms_hom_het_lms <- split(brms_hom_het_lms_areas, brms_hom_het_lms_areas$interval)

brms_hom_het_lms$bottom <- brms_hom_het_lms$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_hom_het_lms$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = loadtype, col = loadtype))+
  geom_segment(data=brms_hom_het_lms_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_hom_het_lms_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_hom_het_lms_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "Zygosity")+
  facet_grid(~loadtype, scales="free_x")+
  scale_x_continuous(breaks=c(-0.9, -0.6, -0.3, 0))+
  scale_fill_manual(values =alpha(c(clr_gerp, clr_high), 0.7)) + #"#E9A69B"
  scale_color_manual(values =c(clr_gerp, clr_high)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(3,"lines"),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> hom_het

hom_het

brms_hom_het_lms_interval <- brms_hom_het_lms_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_hom_het_lms_interval, file = "output/models/intervals/hom_het_gerp45_high.csv", quote=F, row.names = F)

png(file = "plots/main/fig_2b.png", width=800, height=600)
hom_het
dev.off()

#### plot c - GERP per region ####
load(file = "output/models/per_gene_region/lms_total_gerp45_promoter.RData")
gerp_promoter <- brm_load_t
r2_bayes(gerp_promoter)

load(file = "output/models/per_gene_region/lms_total_gerp45_tss.RData")
gerp_tss <- brm_load_t
r2_bayes(gerp_tss)

load(file = "output/models/per_gene_region/lms_total_gerp45_exon.RData")
gerp_exon <- brm_load_t
r2_bayes(gerp_exon)

load(file = "output/models/per_gene_region/lms_total_gerp45_intron.RData")
gerp_intron <- brm_load_t
r2_bayes(gerp_intron)

rm(brm_load_t)

#extract intervals and areas
#total
gerp_promoter_interval <- mcmc_intervals_data(gerp_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_tss_interval <- mcmc_intervals_data(gerp_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_intron_interval <- mcmc_intervals_data(gerp_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_exon_interval <- mcmc_intervals_data(gerp_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

brms_gerps_regions_interval <- rbind(gerp_promoter_interval, gerp_tss_interval, 
                                     gerp_intron_interval, gerp_exon_interval)

brms_gerps_regions_interval$region <- c("Promoter", "TSS", "Intron", "Exon")

#rearrange order for visualization
brms_gerps_regions_interval$region  <- factor(as.factor(brms_gerps_regions_interval$region),
                                              levels= c("Exon", "Intron", 
                                                        "TSS", "Promoter"))

#areas

#total
gerp_promoter_area <- mcmc_areas_data(gerp_promoter) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_tss_area <- mcmc_areas_data(gerp_tss) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_intron_area <- mcmc_areas_data(gerp_intron) %>%
  subset(grepl("scaletotal_load", parameter))

gerp_exon_area <- mcmc_areas_data(gerp_exon) %>%
  subset(grepl("scaletotal_load", parameter))

brms_gerps_regions_area <- rbind(gerp_promoter_area, gerp_tss_area, 
                                 gerp_intron_area, gerp_exon_area)

obs <- nrow(gerp_promoter_area)
brms_gerps_regions_area$region <- c(rep(c("Promoter", "TSS", "Intron", "Exon"), each = obs))

#rearrange order for visualization
brms_gerps_regions_area$region  <- factor(as.factor(brms_gerps_regions_area$region),
                                          levels= c( "Exon","Intron", 
                                                     "TSS", "Promoter"))

### plot

# split by interval
brms_gerps_regions <- split(brms_gerps_regions_area, brms_gerps_regions_area$interval)

brms_gerps_regions$bottom <- brms_gerps_regions$outer %>%
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
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = region, col = region))+
  geom_segment(data=brms_gerps_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_gerps_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_gerps_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  scale_fill_manual(values =alpha(c(clr_gerp, clr_gerp, clr_gerp, clr_gerp), 0.7)) + #
  scale_color_manual(values =c(clr_gerp, clr_gerp, clr_gerp, clr_gerp)) +
  labs(x = expression("Standardised"~beta), y = "Region", title = "GERP")+
  annotate("text", label = "Regulatory", y = 3.2, x = 0.5, size = 5.5, angle = -90, col = "grey30") +
  geom_segment(x = 0.4, y = 2, yend = 4.5,
               col = "grey30")+
  annotate("text", label = "Coding", y = 1.17, x = 0.5, size = 5.5, angle = -90, col = "grey30") +
  geom_segment(x = 0.4, y = 0.5, yend = 1.8,
               col = "grey30")+
  scale_x_continuous(labels = c("-0.25", "0.00", "0.25",  ""), breaks = c(-0.25, 0, 0.25, 0.5))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> posterior_gerpregions

posterior_gerpregions

brms_gerps_regions_interval <- brms_gerps_regions_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_gerps_regions_interval, file = "output/models/intervals/regions_gerp45.csv", quote=F, row.names = F)

png(file = "plots/main/fig_2c.png", width=600, height=800)
posterior_gerpregions
dev.off()

#### plot d - high per region ####
load(file = "output/models/per_gene_region/lms_total_high_promoter.RData")
high_promoter <- brm_load_t
r2_bayes(high_promoter)

load(file = "output/models/per_gene_region/lms_total_high_tss.RData")
high_tss <- brm_load_t
r2_bayes(high_tss)

load(file = "output/models/per_gene_region/lms_total_high_exon.RData")
high_exon <- brm_load_t
r2_bayes(high_exon)

load(file = "output/models/per_gene_region/lms_total_high_intron.RData")
high_intron <- brm_load_t
r2_bayes(high_intron)

rm(brm_load_t)

#extract intervals and areas
#total
high_promoter_interval <- mcmc_intervals_data(high_promoter, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

high_tss_interval <- mcmc_intervals_data(high_tss, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

high_intron_interval <- mcmc_intervals_data(high_intron, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

high_exon_interval <- mcmc_intervals_data(high_exon, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("scaletotal_load", parameter))

brms_high_regions_interval <- rbind(high_promoter_interval, high_tss_interval, 
                                    high_intron_interval, high_exon_interval)

brms_high_regions_interval$region <- c("Promoter", "TSS", "Intron", "Exon")

#rearrange order for visualization
brms_high_regions_interval$region  <- factor(as.factor(brms_high_regions_interval$region),
                                             levels= c("Exon", "Intron", 
                                                       "TSS", "Promoter"))

#areas

#total
high_promoter_area <- mcmc_areas_data(high_promoter) %>%
  subset(grepl("scaletotal_load", parameter))

high_tss_area <- mcmc_areas_data(high_tss) %>%
  subset(grepl("scaletotal_load", parameter))

high_intron_area <- mcmc_areas_data(high_intron) %>%
  subset(grepl("scaletotal_load", parameter))

high_exon_area <- mcmc_areas_data(high_exon) %>%
  subset(grepl("scaletotal_load", parameter))

brms_high_regions_area <- rbind(high_promoter_area, high_tss_area, 
                                high_intron_area, high_exon_area)

obs <- nrow(high_promoter_area)
brms_high_regions_area$region <- c(rep(c("Promoter", "TSS", "Intron", "Exon"), each = obs))

#rearrange order for visualization
brms_high_regions_area$region  <- factor(as.factor(brms_high_regions_area$region),
                                         levels= c( "Exon","Intron", 
                                                    "TSS", "Promoter"))

### plot

# split by interval
brms_high_regions <- split(brms_high_regions_area, brms_high_regions_area$interval)

brms_high_regions$bottom <- brms_high_regions$outer %>%
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
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = region, col = region))+
  geom_segment(data=brms_high_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_high_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_high_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  scale_fill_manual(values =alpha(c(clr_high, clr_high, clr_high, clr_high), 0.7)) + #
  scale_color_manual(values =c(clr_high, clr_high, clr_high, clr_high)) +
  labs(x = expression("Standardised"~beta), y = "Region", title = "SnpEff")+
  annotate("text", label = "Regulatory", y = 3.2, x = 0.2, size = 5.5, angle = -90, col = "grey30") +
  geom_segment(x = 0.15, y = 2, yend = 4.5,
               col = "grey30")+
  annotate("text", label = "Coding", y = 1.17, x = 0.2, size = 5.5, angle = -90, col = "grey30") +
  geom_segment(x = 0.15, y = 0.5, yend = 1.8,
               col = "grey30")+
  scale_x_continuous(labels = c("-0.25", "0.00", ""), breaks = c(-0.25, 0,0.75))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> posterior_highregions

posterior_highregions

brms_high_regions_interval <- brms_high_regions_interval %>% mutate(across(where(is.numeric), ~round(., 2)))
write.csv(brms_high_regions_interval, file = "output/models/intervals/regions_high.csv", quote=F, row.names = F)

png(file = "plots/main/fig_2d.png", width=600, height=800)
posterior_highregions
dev.off()

#### combine ####

cowplot::plot_grid(totals,  hom_het, posterior_gerpregions,  posterior_highregions, 
                   ncol = 2, align = "hv", axis = "lb",rel_heights = c(0.5, 1),
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig

png("plots/main/fig_4.png", height = 1000, width = 1000)
fig
dev.off()
