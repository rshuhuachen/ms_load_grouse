
### load packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

## load pheno data and load
load("data/phenotypes/phenotypes_lifetime.RData") #LMS
load("output/load/all_loads_combined_da_nosex_29scaf.RData") #loads no sex chr only 30 scaf

pheno_load <- left_join(pheno_wide, load, by = "id")


##### plots a+b - GERP models, total and hom/het #####

#load model
load(file = "output/models/total_hom_het/lms_total_gerp45.RData")
brm_load_t_gerp5_lms <- brm_load_t
load(file = "output/models/total_hom_het/lms_het_hom_gerp45.RData")
brm_load_rp_gerp5_lms <- brm_load_het_hom

##### plot b - raw data GERP models, total and realized/potential #####
effect_gerpt<- conditional_effects(brm_load_t_gerp5_lms)
effect_gerpt <- effect_gerpt$total_load
effect_gerprp<- conditional_effects(brm_load_rp_gerp5_lms)
effect_gerpr <- effect_gerprp$hom_load
effect_gerpp <- effect_gerprp$het_load

## combine all in one df

effect_gerpt$var <- "Total"
effect_gerpr$var <- "Homozygous"
effect_gerpp$var <- "Heterozygous"

effect_gerpt<- effect_gerpt %>% dplyr::select(c(total_load, var, estimate__, lower__, upper__))
effect_gerpr<- effect_gerpr %>% dplyr::select(c(hom_load, var, estimate__, lower__, upper__))
effect_gerpp<- effect_gerpp %>% dplyr::select(c(het_load, var, estimate__, lower__, upper__))

names(effect_gerpt)[1] <- "load"
names(effect_gerpr)[1] <- "load"
names(effect_gerpp)[1] <- "load"

effect_gerp <- rbind(effect_gerpt, effect_gerpr, effect_gerpp)

## make a long df with the loads and pheno
pheno_load_long <- pheno_load %>% filter(loadtype == "gerp45") %>% dplyr::select(id, LMS_min, total_load, hom_load, het_load)
pheno_load_long <- gather(pheno_load_long, var, load, total_load:het_load, factor_key=T)
pheno_load_long$var <- gsub("total_load", "Total", pheno_load_long$var)
pheno_load_long$var <- gsub("hom_load", "Homozygous", pheno_load_long$var)
pheno_load_long$var <- gsub("het_load", "Heterozygous", pheno_load_long$var)

pheno_load_long$var <- factor(pheno_load_long$var , levels = c("Total", "Homozygous", "Heterozygous"))
effect_gerp$var <- factor(effect_gerp$var , levels = c("Total", "Homozygous", "Heterozygous"))
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

#hom/het

ggplot() + 
  geom_point(data=subset(pheno_load_long, var != "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  geom_line(data=subset(effect_gerp, var != "Total"), aes(x=load, y=estimate__), col = clr_gerp,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_gerp, var != "Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_gerp, 
              alpha= 0.3) +
  facet_wrap(~var, scales = "free", switch = "y", ncol = 1, nrow = 3) + 
  labs(x = "GERP", y = "Lifetime mating success") + 
  coord_cartesian(ylim = c(1, 30))+
  scale_y_log10(breaks = c(1, 3, 10, 30))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_gerp_hom_het

##### plot c+d - SNPEFF models, total and realized/potential #####

#load model
load(file = "output/models/total_hom_het/lms_total_high.RData")
brm_load_t_high5_lms <- brm_load_t
load(file = "output/models/total_hom_het/lms_het_hom_high.RData")
brm_load_rp_high5_lms <- brm_load_het_hom

##### plot d - raw data high models, total and realized/potential #####
effect_hight<- conditional_effects(brm_load_t_high5_lms)
effect_hight <- effect_hight$total_load
effect_highrp<- conditional_effects(brm_load_rp_high5_lms)
effect_highr <- effect_highrp$hom_load
effect_highp <- effect_highrp$het_load

## combine all in one df

effect_hight$var <- "Total"
effect_highr$var <- "Homozygous"
effect_highp$var <- "Heterozygous"

effect_hight<- effect_hight %>% dplyr::select(c(total_load, var, estimate__, lower__, upper__))
effect_highr<- effect_highr %>% dplyr::select(c(hom_load, var, estimate__, lower__, upper__))
effect_highp<- effect_highp %>% dplyr::select(c(het_load, var, estimate__, lower__, upper__))

names(effect_hight)[1] <- "load"
names(effect_highr)[1] <- "load"
names(effect_highp)[1] <- "load"

effect_high <- rbind(effect_hight, effect_highr, effect_highp)

## make a long df with the loads and pheno
pheno_load_long <- pheno_load %>% filter(loadtype == "high") %>% dplyr::select(id, LMS_min, total_load, hom_load, het_load)
pheno_load_long <- gather(pheno_load_long, var, load, total_load:het_load, factor_key=T)
pheno_load_long$var <- gsub("total_load", "Total", pheno_load_long$var)
pheno_load_long$var <- gsub("hom_load", "Homozygous", pheno_load_long$var)
pheno_load_long$var <- gsub("het_load", "Heterozygous", pheno_load_long$var)

pheno_load_long$var <- factor(pheno_load_long$var , levels = c("Total", "Homozygous", "Heterozygous"))
effect_high$var <- factor(effect_high$var , levels = c("Total", "Homozygous", "Heterozygous"))
pacman::p_load(scales)

#divide up between total and then masked/expressed

ggplot() + 
  geom_point(data=subset(pheno_load_long, var == "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  geom_line(data=subset(effect_high, var == "Total"), aes(x=load, y=estimate__), col = clr_high,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_high, var == "Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_high, 
              alpha= 0.3) +
  labs(x = "Total SnpEff load", y = "Lifetime mating success") + 
  coord_cartesian(ylim = c(1, 30))+
  scale_y_log10(breaks = c(1, 3, 10, 30))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_high_total

#hom/het

ggplot() + 
  geom_point(data=subset(pheno_load_long, var != "Total"), aes(x=load, y=LMS_min), size = 2.5) + 
  geom_line(data=subset(effect_high, var != "Total"), aes(x=load, y=estimate__), col = clr_high,
            linewidth = 1) +
  geom_ribbon(data=subset(effect_high, var != "Total"), aes(x=load, y=estimate__, ymin=lower__, ymax=upper__), fill = clr_high, 
              alpha= 0.3) +
  facet_wrap(~var, scales = "free", switch = "y", ncol = 1, nrow = 3) + 
  labs(x = "SnpEff", y = "Lifetime mating success") + 
  coord_cartesian(ylim = c(1, 30))+
  scale_y_log10(breaks = c(1, 3, 10, 30))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") -> fig_raw_high_hom_het

#### Combine #### 

cowplot::plot_grid(fig_raw_gerp_total, 
                   fig_raw_high_total, 
                   fig_raw_gerp_hom_het,
                   fig_raw_high_hom_het, ncol = 2, 
                   rel_heights = c(0.6,1),
                 #  align = "hv", axis = "lb",
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig_rawdata

png(file = "plots/sup/rawdata_total_homhet_gerp_snp.png", height=1000, width=800)
fig_rawdata
dev.off()
