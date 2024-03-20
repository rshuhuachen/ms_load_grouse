
### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

alpha_total = 0.8
alpha_expressed = 0.6
alpha_masked = 0.4

## load pheno data and load
load("data/phenotypes/phenotypes_lifetime.RData") #LMS
load("output/load/all_loads_combined_da_nosex_29scaf.RData") #loads no sex chr only 30 scaf

pheno_load <- left_join(pheno_wide, loads, by = "id")

#### plot a: inbreeding depression ####

#load model
load(file = "output/models/genetic_quality/lms_froh_zi_29scaf.RData")

effect_froh <- conditional_effects(brm_froh_lms)

ggplot() + 
  geom_point(data=pheno_load, aes(froh1, LMS_min), size = 2.5) + 
  scale_y_log10(labels= c("1","3","10","30"), breaks = c(1,3,10,30)) + 
  geom_line(data=effect_froh$froh1, aes(x=froh1, y=estimate__), color=clr_froh, linewidth = 1) +
  geom_ribbon(data=effect_froh$froh1, aes(x=froh1, y=estimate__, ymin=lower__, ymax=upper__), alpha= 0.3, fill=clr_froh) +
  labs(x = expression(F[ROH]), y = "Lifetime mating success") +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig4b_fit_froh

##### plot b - GERP models, total and realized/potential #####

#load model
load(file = "output/models/genetic_quality/lms_gerp5_total_add_zi_29scaf.RData")
brm_load_t_gerp5_lms <- brm_load_t_add_high_lms
load(file = "output/models/genetic_quality/lms_gerp5_realized_potential_zi_29scaf.RData")
brm_load_rp_gerp5_lms <- brm_load_rp_high_lms


##### plot d - raw data GERP models, total and realized/potential #####
effect_gerpt<- conditional_effects(brm_load_t_gerp5_lms)
effect_gerpt <- effect_gerpt$gerp_Lt_additive_cat5
effect_gerprp<- conditional_effects(brm_load_rp_gerp5_lms)
effect_gerpr <- effect_gerprp$gerp_Lr_cat5
effect_gerpp <- effect_gerprp$gerp_Lp_cat5

## combine all in one df

effect_gerpt$var <- "Total"
effect_gerpr$var <- "Expressed"
effect_gerpp$var <- "Masked"

effect_gerpt<- effect_gerpt %>% dplyr::select(c(gerp_Lt_additive_cat5, var, estimate__, lower__, upper__))
effect_gerpr<- effect_gerpr %>% dplyr::select(c(gerp_Lr_cat5, var, estimate__, lower__, upper__))
effect_gerpp<- effect_gerpp %>% dplyr::select(c(gerp_Lp_cat5, var, estimate__, lower__, upper__))

names(effect_gerpt)[1] <- "load"
names(effect_gerpr)[1] <- "load"
names(effect_gerpp)[1] <- "load"

effect_gerp <- rbind(effect_gerpt, effect_gerpr, effect_gerpp)

## make a long df with the loads and pheno
pheno_load_long <- pheno_load %>% dplyr::select(id, LMS_min, gerp_Lt_additive_cat5, gerp_Lr_cat5, gerp_Lp_cat5)
pheno_load_long <- gather(pheno_load_long, var, load, gerp_Lt_additive_cat5:gerp_Lp_cat5, factor_key=T)
pheno_load_long$var <- gsub("gerp_Lt_additive_cat5", "Total", pheno_load_long$var)
pheno_load_long$var <- gsub("gerp_Lr_cat5", "Expressed", pheno_load_long$var)
pheno_load_long$var <- gsub("gerp_Lp_cat5", "Masked", pheno_load_long$var)

pheno_load_long$var <- factor(pheno_load_long$var , levels = c("Total", "Expressed", "Masked"))
effect_gerp$var <- factor(effect_gerp$var , levels = c("Total", "Expressed", "Masked"))
pacman::p_load(scales)

#divide up between total and then masked/expressed

ggplot() + 
  geom_point(data=subset(pheno_load_long, var == "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  geom_line(data=subset(effect_gerp, var == "Total"), aes(x=load, y=estimate__), col = clr_gerp,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_gerp, var == "Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_gerp, 
              alpha= 0.3) +
  labs(x = "Total GERP load", y = "Lifetime mating success") + 
  coord_cartesian(ylim = c(1, 30))+
  scale_y_log10(breaks = c(1, 3, 10, 30))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_gerp_total

#realized/masked

ggplot() + 
  geom_point(data=subset(pheno_load_long, var != "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  geom_line(data=subset(effect_gerp, var != "Total"), aes(x=load, y=estimate__), col = clr_gerp,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_gerp, var != "Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_gerp, 
              alpha= 0.3) +
  facet_wrap(~var, scales = "free", switch = "y", ncol = 1, nrow = 3) + 
  labs(x = "GERP load", y = "Lifetime mating success") + 
  coord_cartesian(ylim = c(1, 30))+
  scale_y_log10(breaks = c(1, 3, 10, 30))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_gerp_express_mask


##### plot c - SNPEFF models, total and realized/potential #####

#load model
load(file = "output/models/genetic_quality/lms_snpef_high_realized_potential_zi_29scaf.RData")
load(file = "output/models/genetic_quality/lms_snpef_high_total_add_zi_29scaf.RData")

##### plot f - raw data SNPeff models, total and realized/potential #####

effect_hight<- conditional_effects(brm_load_t_add_high_lms)
effect_highrp<- conditional_effects(brm_load_rp_high_lms)

effect_hight<- conditional_effects(brm_load_t_add_high_lms)
effect_hight <- effect_hight$high_Lt_additive
effect_highrp<- conditional_effects(brm_load_rp_high_lms)
effect_highr <- effect_highrp$high_Lr
effect_highp <- effect_highrp$high_Lp

## combine all in one df

effect_hight$var <- "Total"
effect_highr$var <- "Expressed"
effect_highp$var <- "Masked"

effect_hight<- effect_hight %>% dplyr::select(c(high_Lt_additive, var, estimate__, lower__, upper__))
effect_highr<- effect_highr %>% dplyr::select(c(high_Lr, var, estimate__, lower__, upper__))
effect_highp<- effect_highp %>% dplyr::select(c(high_Lp, var, estimate__, lower__, upper__))

names(effect_hight)[1] <- "load"
names(effect_highr)[1] <- "load"
names(effect_highp)[1] <- "load"

effect_high <- rbind(effect_hight, effect_highr, effect_highp)

## make a long df with the loads and pheno
pheno_load_long_high <- pheno_load %>% dplyr::select(id, LMS_min, high_Lt_additive, high_Lr, high_Lp)
pheno_load_long_high <- gather(pheno_load_long_high, var, load, high_Lt_additive:high_Lp, factor_key=T)
pheno_load_long_high$var <- gsub("high_Lt_additive", "Total", pheno_load_long_high$var)
pheno_load_long_high$var <- gsub("high_Lr", "Expressed", pheno_load_long_high$var)
pheno_load_long_high$var <- gsub("high_Lp", "Masked", pheno_load_long_high$var)

pheno_load_long_high$var <- factor(pheno_load_long_high$var , levels = c("Total", "Expressed", "Masked"))
effect_high$var <- factor(effect_high$var , levels = c("Total", "Expressed", "Masked"))

ggplot() + 
  geom_point(data=subset(pheno_load_long_high, var == "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  scale_y_log10()+
  geom_line(data=subset(effect_high, var=="Total"), aes(x=load, y=estimate__), col = clr_high,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_high, var=="Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_high,
              alpha= 0.3) +
  labs(x = "Total high impact load", y = "Lifetime mating success") + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_high_total

ggplot() + 
  geom_point(data=subset(pheno_load_long_high, var != "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  scale_y_log10()+
  geom_line(data=subset(effect_high, var!="Total"), aes(x=load, y=estimate__), col = clr_high,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_high, var!="Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_high,
              alpha= 0.3) +
  facet_wrap(~var, scales = "free", switch = "y", ncol = 1, nrow = 3) + 
  labs(x = "High impact load", y = "Lifetime mating success") + 
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_high_express_mask

#### Combine #### 

cowplot::plot_grid(fig4b_fit_froh,  fig_raw_gerp_total, fig_raw_high_total,
                   fig_raw_gerp_express_mask, fig_raw_high_express_mask, ncol = 2, 
                   align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig_rawdata

png(file = "plots/final/Supp_4_rawdata.png", height=1000, width=800)
fig_rawdata
dev.off()
