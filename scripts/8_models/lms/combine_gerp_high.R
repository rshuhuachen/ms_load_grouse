### Packages ####
pacman::p_load(tidyverse, brms, bayesplot, data.table, ggridges)

# function to calculate load
source("scripts/7_calculate_load/function_calculate_load.R")

# theme
source("scripts/theme_ggplot.R")

### Merge the two mutation types together ####
load(file = "output/load/snpeff/snpeff_high.RData")
load(file = "output/load/gerp/gerp_over4.RData")

names(gerp)[1:12]
names(snpeff)[1:12]

snpeff_c <- data.frame(chr = snpeff$CHROM,
                       start = snpeff$POS -1,
                       pos = snpeff$POS,
                       neutral_rate_n = NA,
                       rs_score = NA,
                       ancestral = snpeff$REF,
                       derived = snpeff$ALT,
                       qual = snpeff$QUAL,
                       info = snpeff$INFO,
                       format = snpeff$FORMAT)

snpeff_c <- cbind(snpeff_c, snpeff[c(10:199)])

combined <- rbind(gerp, snpeff_c)

load_combined <- calculate_load_gerp(vcf = combined, output_vcf = FALSE, loadtype= "gerp45_plus_high")

# load pheno data lifetime
load(file = "output/load/pheno_loads_lifetime.RData")

# combine
pheno_wide_load <- left_join(pheno_wide, load_combined, by = "id")
pheno_wide_load <- subset(pheno_wide_load, !is.na(total_load)) #some ids without genotypes, excluded for wgr

#### Model ####
summary(glmmTMB::glmmTMB(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                         family = "poisson", ziformula = ~1))


# brms model
brm_load_t <- brm(LMS_min ~ scale(total_load) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_t, file = "output/models/total_hom_het/lms_total_gerp45_plus_high.RData")

brm_load_gerp_high <- brm_load_t
bayesplot::mcmc_intervals_data(brm_load_gerp_high, prob = 0.8, prob_outer = 0.95) %>%
  write.csv(file = "output/models/total_hom_het/lms_total_gerp45_plus_high_intervals.csv")

#### Plot with gerp and snpeff ####
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
  scale_fill_manual(values =alpha(c(clr_highlight, clr_gerp, clr_high), 0.7)) +
  scale_color_manual(values =c(clr_highlight, clr_gerp, clr_high)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plot_posteriors_all

ggsave(plot_posteriors_all, file = "plots/sup/combined_high_gerp_load.png", width=10, height=12)



#### Approach 2: both in the same model ####
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

#### model ####
brm_load_t <- brm(LMS_min ~ scale(total_load_high) + scale(total_load_gerp45) + core + (1|site), data = pheno_wide_load,
                  family = "zero_inflated_poisson",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = 1000000, thin = 1000, warmup = 500000, seed = 1908)

save(brm_load_t, file = "output/models/total_hom_het/lms_total_gerp45_high_sep.RData")

brm_load_gerp_high_sep <- brm_load_t
bayesplot::mcmc_intervals_data(brm_load_gerp_high_sep, prob = 0.8, prob_outer = 0.95) %>%
  write.csv(file = "output/models/total_hom_het/lms_total_gerp45_high_sep_intervals.csv")

#### Plot with gerp and snpeff ####
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
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> plot_posteriors_all

ggsave(plot_posteriors_all, file = "plots/sup/separate_one_model_high_gerp_load.png", width=10, height=12)


