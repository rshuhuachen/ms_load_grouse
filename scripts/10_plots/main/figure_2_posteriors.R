###### This script contains the plot that contains the posterior distributions and raw data ####

### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

alpha_total = 0.8
alpha_expressed = 0.6
alpha_masked = 0.4

#### plot a: inbreeding depression ####

#load model
load(file = "output/models/genetic_quality/lms_froh_zi_29scaf.RData")

#extract intervals and areas
brms_froh_lms_interval <- mcmc_intervals_data(brm_froh_lms, prob =0.8, prob_outer = 0.95)

brms_froh_lms_interval$model <- c(rep("Inbreeding", nrow(mcmc_intervals_data(brm_froh_lms))))

brms_froh_lms_interval <- subset(brms_froh_lms_interval, grepl("b_scalefroh1", parameter))
brms_froh_lms_interval$model  <- as.factor(brms_froh_lms_interval$model )
#area
brms_froh_lms_area <- mcmc_areas_data(brm_froh_lms)

brms_froh_lms_area$model <- c(rep("Inbreeding", nrow(mcmc_areas_data(brm_froh_lms))))

brms_froh_lms_area <- subset(brms_froh_lms_area, grepl("b_scalefroh1", parameter))
brms_froh_lms_area$model <- as.factor(brms_froh_lms_area$model)
### plot

# split by interval
brms_froh_lms <- split(brms_froh_lms_area, brms_froh_lms_area$interval)

brms_froh_lms$bottom <- brms_froh_lms$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

ggplot(data = brms_froh_lms$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_froh_lms_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_froh_lms_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_froh_lms_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  scale_fill_manual(values =alpha(clr_froh, 0.7)) +
  scale_color_manual(values =clr_froh) +
  labs(x = "Standardised beta coefficient")+
  xlim(-1.3, 0.5)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")  -> fig4a_posteriors_froh

write_tsv(brms_froh_lms_interval, file="output/models/intervals/froh_lms_intervals.tsv")

##### plot c - GERP models, total and realized/potential #####

#load model
load(file = "output/models/genetic_quality/lms_gerp5_total_add_zi_29scaf.RData")
brm_load_t_gerp5_lms <- brm_load_t_add_high_lms
load(file = "output/models/genetic_quality/lms_gerp5_realized_potential_zi_29scaf.RData")
brm_load_rp_gerp5_lms <- brm_load_rp_high_lms

#extract intervals and areas
brms_gerp5_lms_interval_total <- mcmc_intervals_data(brm_load_t_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_gerp5_lms_interval_expressed <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lr_cat5", parameter))

brms_gerp5_lms_interval_potential <- mcmc_intervals_data(brm_load_rp_gerp5_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalegerp_Lp_cat5", parameter))

brms_gerp5_lms_interval <- rbind(brms_gerp5_lms_interval_total,
                                 brms_gerp5_lms_interval_expressed,
                                 brms_gerp5_lms_interval_potential)
                                 

brms_gerp5_lms_interval$model <- c(rep("Total", nrow(brms_gerp5_lms_interval_total)),
                                   rep("Expressed", nrow(brms_gerp5_lms_interval_expressed)),
                                   rep("Masked", nrow(brms_gerp5_lms_interval_potential)))


brms_gerp5_lms_interval$model  <- factor(as.factor(brms_gerp5_lms_interval$model),
                                         levels= c("Masked", "Expressed", "Total"))
#area
brms_gerp5_lms_areas_total <- mcmc_areas_data(brm_load_t_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lt_additive_cat5", parameter))

brms_gerp5_lms_areas_expressed <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lr_cat5", parameter))

brms_gerp5_lms_areas_potential <- mcmc_areas_data(brm_load_rp_gerp5_lms) %>%
  subset(grepl("b_scalegerp_Lp_cat5", parameter))

brms_gerp5_lms_areas <- rbind(brms_gerp5_lms_areas_total,
                              brms_gerp5_lms_areas_expressed, 
                              brms_gerp5_lms_areas_potential)

brms_gerp5_lms_areas$model <- c(rep("Total", nrow(brms_gerp5_lms_areas_total)),
                                   rep("Expressed", nrow(brms_gerp5_lms_areas_expressed)),
                                   rep("Masked", nrow(brms_gerp5_lms_areas_potential)))


brms_gerp5_lms_areas$model  <- factor(as.factor(brms_gerp5_lms_areas$model),
                                         levels= c("Masked", "Expressed", "Total"))


### plot

# split by interval
brms_gerp5_lms <- split(brms_gerp5_lms_areas, brms_gerp5_lms_areas$interval)

brms_gerp5_lms$bottom <- brms_gerp5_lms$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_gerp5_lms$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_gerp5_lms_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_gerp5_lms_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_gerp5_lms_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Load type", title= "GERP score 4-5")+
  scale_fill_manual(values =alpha(c( "#E9A69B", "#DE7A68",clr_gerp), 0.7)) +
  scale_color_manual(values =c("#E9A69B", "#DE7A68",clr_gerp)) +
  xlim(-1.3, 0.5)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig4c_posteriors_gerp

write_tsv(brms_gerp5_lms_interval, file="output/models/intervals/gerp5_tem_lms_intervals.tsv")

fig4c_posteriors_gerp

##### plot e - SNPEFF models, total and realized/potential #####

#load model
load(file = "output/models/genetic_quality/lms_snpef_high_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_snpef_high_total_add_1mil_zi_29scaf.RData")

#extract intervals and areas
brms_high_lms_interval_total <- mcmc_intervals_data(brm_load_t_add_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_high_lms_interval_expressed <- mcmc_intervals_data(brm_load_rp_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lr", parameter))

brms_high_lms_interval_potential <- mcmc_intervals_data(brm_load_rp_high_lms, prob =0.8, prob_outer = 0.95) %>%
  subset(grepl("b_scalehigh_Lp", parameter))

brms_high_lms_interval <- rbind(brms_high_lms_interval_total,
                                 brms_high_lms_interval_expressed,
                                 brms_high_lms_interval_potential)


brms_high_lms_interval$model <- c(rep("Total", nrow(brms_high_lms_interval_total)),
                                   rep("Expressed", nrow(brms_high_lms_interval_expressed)),
                                   rep("Masked", nrow(brms_high_lms_interval_potential)))


brms_high_lms_interval$model  <- factor(as.factor(brms_high_lms_interval$model),
                                         levels= c("Masked", "Expressed", "Total"))
#area
brms_high_lms_areas_total <- mcmc_areas_data(brm_load_t_add_high_lms) %>%
  subset(grepl("b_scalehigh_Lt_additive", parameter))

brms_high_lms_areas_expressed <- mcmc_areas_data(brm_load_rp_high_lms) %>%
  subset(grepl("b_scalehigh_Lr", parameter))

brms_high_lms_areas_potential <- mcmc_areas_data(brm_load_rp_high_lms) %>%
  subset(grepl("b_scalehigh_Lp", parameter))

brms_high_lms_areas <- rbind(brms_high_lms_areas_total,
                              brms_high_lms_areas_expressed, 
                              brms_high_lms_areas_potential)

brms_high_lms_areas$model <- c(rep("Total", nrow(brms_high_lms_areas_total)),
                                rep("Expressed", nrow(brms_high_lms_areas_expressed)),
                                rep("Masked", nrow(brms_high_lms_areas_potential)))


brms_high_lms_areas$model  <- factor(as.factor(brms_high_lms_areas$model),
                                      levels= c("Masked", "Expressed", "Total"))


### plot

# split by interval
brms_high_lms <- split(brms_high_lms_areas, brms_high_lms_areas$interval)

brms_high_lms$bottom <- brms_high_lms$outer %>%
  group_by(!!! groups) %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot

### divide up between the two traits
ggplot(data = brms_high_lms$outer) +  
  aes(x = .data$x, y = .data$model) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = model, col = model))+
  geom_segment(data=brms_high_lms_interval, aes(x = l, xend = h, yend = model), col = "black", linewidth=3)+
  geom_segment(data=brms_high_lms_interval, aes(x = ll, xend = hh, yend = model), col = "black")+
  geom_point(data=brms_high_lms_interval, aes(x = m, y = model), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Load type", title = "High impact SnpEff")+
  scale_fill_manual(values =alpha(c( "#5B98AE","#437689", clr_high), 0.7)) +
  scale_color_manual(values =c("#5B98AE","#437689", clr_high)) +
  xlim(-1.3, 0.5)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig4e_posteriors_high

write_tsv(brms_high_lms_interval, file="output/models/intervals/high_tem_lms_intervals.tsv")

#### Add total load per region in same plot ####
#### GERP per different region ####
#combine total, expressed and realized in one
load(file = "output/models/per_gene_region/lms_total_gerp5_region_promoter_zi_29scaf.RData")
gerp5_t_promoter <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_tss_zi_29scaf.RData")
gerp5_t_tss <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_gene_zi_29scaf.RData")
gerp5_t_gene <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_down_zi_29scaf.RData")
gerp5_t_down <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_up_zi_29scaf.RData")
gerp5_t_up <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_exon_zi_29scaf.RData")
gerp5_t_exon <- brm_region
load(file = "output/models/per_gene_region/lms_total_gerp5_region_intron_zi_29scaf.RData")
gerp5_t_intron <- brm_region

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
brms_gerps_regions_interval <- rbind(gerp5_t_promoter_interval, gerp5_t_tss_interval, gerp5_t_gene_interval,
                                     gerp5_t_down_interval, gerp5_t_up_interval, 
                                     gerp5_t_intron_interval, gerp5_t_exon_interval)


brms_gerps_regions_interval$region <- c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon")
#rearrange order for visualization
brms_gerps_regions_interval$region  <- factor(as.factor(brms_gerps_regions_interval$region),
                                              levels= c("Exon", "Intron", "Upstream", "Downstream", "Gene body",
                                                        "TSS", "Promoter"))

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

brms_gerps_regions_area <- rbind(gerp5_t_promoter_area, gerp5_t_tss_area, gerp5_t_gene_area,
                                 gerp5_t_down_area, gerp5_t_up_area, 
                                 gerp5_t_intron_area, gerp5_t_exon_area)


obs <- nrow(gerp5_t_promoter_area)
brms_gerps_regions_area$region <- c(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"), each = obs))

#rearrange order for visualization
brms_gerps_regions_area$region  <- factor(as.factor(brms_gerps_regions_area$region),
                                          levels= c( "Intron", "Exon","Upstream", "Downstream", "Gene body",
                                                     "TSS", "Promoter"))

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
  geom_ridgeline(aes(scale = 0.4, height = scaled_density), fill = alpha(clr_gerp, 0.7), col = clr_gerp)+
  geom_segment(data=brms_gerps_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_gerps_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_gerps_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Region", title = "Total load GERP score 4-5")+
  xlim(-0.42, 0.42)+ 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> posterior_gerpregions_t
posterior_gerpregions_t


#### SNPeff per different region ####
#combine total, expressed and realized in one
load(file = "output/models/per_gene_region/lms_total_high_region_promoter_zi_29scaf.RData")
high_t_promoter <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_tss_zi_29scaf.RData")
high_t_tss <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_gene_zi_29scaf.RData")
high_t_gene <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_down_zi_29scaf.RData")
high_t_down <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_up_zi_29scaf.RData")
high_t_up <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_exon_zi_29scaf.RData")
high_t_exon <- brm_region
load(file = "output/models/per_gene_region/lms_total_high_region_intron_zi_29scaf.RData")
high_t_intron <- brm_region

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


brms_high_regions_interval <- rbind(high_t_promoter_interval, high_t_tss_interval, high_t_gene_interval,
                                    high_t_down_interval, high_t_up_interval, 
                                    high_t_intron_interval, high_t_exon_interval)


brms_high_regions_interval$region <- c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon")

#rearrange order for visualization
brms_high_regions_interval$region  <- factor(as.factor(brms_high_regions_interval$region),
                                             levels= c( "Intron", "Exon","Upstream", "Downstream", "Gene body",
                                                        "TSS", "Promoter"))


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

brms_high_regions_area <- rbind(high_t_promoter_area, high_t_tss_area, high_t_gene_area,
                                high_t_down_area, high_t_up_area, 
                                high_t_intron_area, high_t_exon_area)


obs <- nrow(high_t_promoter_area)
brms_high_regions_area$region <- c(rep(c("Promoter", "TSS", "Gene body", "Downstream", "Upstream", "Intron", "Exon"), each = obs))

#rearrange order for visualization
brms_high_regions_area$region  <- factor(as.factor(brms_high_regions_area$region),
                                         levels= c("Intron", "Exon",  "Upstream", "Downstream", "Gene body",
                                                   "TSS", "Promoter"))

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
  geom_ridgeline(aes(scale = 0.4, height = scaled_density), fill = alpha(clr_high, 0.7), col = clr_high)+
  geom_segment(data=brms_high_regions_interval, aes(x = l, xend = h, yend = region), col = "black", linewidth=3)+
  geom_segment(data=brms_high_regions_interval, aes(x = ll, xend = hh, yend = region), col = "black")+
  geom_point(data=brms_high_regions_interval, aes(x = m, y = region), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "Region", title = "Total load high impact SnpEff")+
  xlim(-0.42, 0.42)+ 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> posterior_high_regions_t
posterior_high_regions_t

#### also add gerp categories ####

#combine total, expressed and realized in one
load(file = "output/models/genetic_quality/lms_gerp1_total_add_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp1_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp2_total_add_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp2_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp3_total_add_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp3_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp4_total_add_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp4_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp5_total_add_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_gerp5_realized_potential_zi_29scaf.RData")
brm_load_t_add_gerp5_lms <- brm_load_t_add_high_lms
brm_load_rp_gerp5_lms <- brm_load_rp_high_lms
rm(brm_load_t_add_high_lms)
rm(brm_load_rp_high_lms)

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

write_tsv(brms_gerps_interval, file="output/models/intervals/allgerpcats_tem_lms_intervals.tsv")

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

ggplot(data = subset(brms_gerps$outer, cat != "0-1" & cat != "1-2" & 
                       loadtype == "Total")) +  
  aes(x = .data$x, y = .data$cat) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = loadtype, col = loadtype))+
  geom_segment(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"& loadtype == "Total"), aes(x = l, xend = h, yend = cat), col = "black", linewidth=3)+
  geom_segment(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"& loadtype == "Total"), aes(x = ll, xend = hh, yend = cat), col = "black")+
  geom_point(data=subset(brms_gerps_interval,cat != "0-1" & cat != "1-2"& loadtype == "Total"), aes(x = m, y = cat), fill="white",  col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = "Standardised beta coefficient", y = "GERP score category", title = "Total load")+
#  facet_wrap(~loadtype, ncol = 3, scales="free_x") +
  xlim(-1.3, 0.5)+
  scale_fill_manual(values =alpha(c(clr_gerp, "#DE7A68", "#E9A69B"), 0.7)) +
  scale_color_manual(values =c(clr_gerp, "#DE7A68", "#E9A69B")) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> supp_posterior_gerpcats

#### combine ####
cowplot::plot_grid(fig4a_posteriors_froh,  supp_posterior_gerpcats,
                   fig4c_posteriors_gerp, 
                   fig4e_posteriors_high,  ncol = 2, align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig4_v2a

cowplot::plot_grid(posterior_gerpregions_t,  posterior_high_regions_t, 
                   ncol = 2, align = "hv", axis = "lb",
                   labels = c("e", "f"), label_fontface = "plain", label_size = 22) -> fig4_v2b

cowplot::plot_grid(fig4_v2a,  fig4_v2b, 
                   ncol = 1) -> fig4_v2


ggsave(fig4_v2, file = "plots/Figure_2.png", height = 16, width = 14)

