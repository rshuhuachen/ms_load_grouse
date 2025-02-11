#### packages ####
pacman::p_load(brms, bayesplot, tidyverse, ggridges, prismatic, ggpubr, 
               extrafont, cowplot, performance, data.table, effects)

source("scripts/theme_ggplot.R")

#### gerp ####

scafs_posteriors_gerp <- list.files(path = "output/models/per_chr", pattern = "lms_total_gerp45*", full.names = T)

scafs_gerp <- list()
for (i in 1:length(scafs_posteriors_gerp)){
  load(file = scafs_posteriors_gerp[i])
  name <- gsub("output/models/per_chr/lms_total_gerp45_", "", scafs_posteriors_gerp[i])
  name <- gsub(".RData", "", name)
  scafs_gerp[[i]] <- brm_load_t
}

scafs_posteriors_names_gerp <- gsub("output/models/per_chr/lms_total_gerp45_", "", scafs_posteriors_gerp)
scafs_posteriors_names_gerp <- gsub(".RData", "", scafs_posteriors_names_gerp)
names(scafs_gerp) <- scafs_posteriors_names_gerp

## nr of snps per chr ##
#system('cat output/ancestral/ltet_filtered_ann_aa.vcf | grep -v "^#" | cut -f 1 | sort | uniq -c > output/n_snps_per_chr.txt')
load("data/genomic/raw/metadata/scaffold_names_dovetail.RData")

n_snp <- read.table("output/n_snps_per_chr.txt")
names(n_snp) <- c("n_snp", "contig")
genome <- left_join(genome[,c("contig", "scaf_nr")], n_snp, by = "contig")
genome$scaf_nr <- as.character(genome$scaf_nr)

## get all intervals and areas 

all_intervals_gerp <- data.frame()
all_areas_gerp <- data.frame()
for (i in 1:length(scafs_gerp)){
  # get interval
  interval <- mcmc_intervals_data(scafs_gerp[[i]], prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
  interval$scaf_nr <- gsub("scaf_", "", names(scafs_gerp[i]))
  interval <- left_join(interval, genome[,c("scaf_nr", "n_snp")], by = "scaf_nr")
  
  # merge
  all_intervals_gerp <- rbind(all_intervals_gerp, interval)
  
  # get area
  area <- mcmc_areas_data(scafs_gerp[[i]], pars = "b_scaletotal_load")
  area$scaf_nr <- gsub("scaf_", "", names(scafs_gerp[i]))
  area <- left_join(area, genome[,c("scaf_nr", "n_snp")], by = "scaf_nr")
  
  
  # merge
  all_areas_gerp <- rbind(all_areas_gerp, area)
}


# split by interval
brms_plota_gerp <- split(all_areas_gerp, all_areas_gerp$interval)

brms_plota_gerp$bottom <- brms_plota_gerp$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot
brms_plota_gerp$outer$scaf_nr <- as.numeric(brms_plota_gerp$outer$scaf_nr)
all_intervals_gerp$scaf_nr <- as.numeric(all_intervals_gerp$scaf_nr)

lm <- lm(m ~ scaf_nr, data = all_intervals_gerp)

all_intervals_gerp <- all_intervals_gerp %>% mutate(sig = case_when(
  m < 0 & l < 0 & h < 0 ~ "neg",
  m > 0 & l > 0 & h > 0 ~ "pos"
))

ggplot(data = all_intervals_gerp, aes(x = m, y = scaf_nr)) +  
  geom_segment(aes(x = l, xend = h, y = scaf_nr), col = "black", linewidth=3)+
  geom_segment(aes(x = ll, xend = hh, y = scaf_nr), col = "black")+
  geom_point(aes(fill = sig), col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "Scaffold", title = "GERP")+
 # geom_smooth(method='lm')+
  scale_fill_manual(values =alpha(c(clr_highlight, clr_grey), 0.7)) +
  scale_color_manual(values =c(clr_highlight, clr_grey)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  coord_flip() -> gerp_load_per_chr

#### snpeff ####
#### gerp ####

scafs_posteriors <- list.files(path = "output/models/per_chr", pattern = "lms_total_high*", full.names = T)

scafs <- list()
for (i in 1:length(scafs_posteriors)){
  load(file = scafs_posteriors[i])
  name <- gsub("output/models/per_chr/lms_total_high", "", scafs_posteriors[i])
  name <- gsub(".RData", "", name)
  scafs[[i]] <- brm_load_t
}

scafs_posteriors_names <- gsub("output/models/per_chr/lms_total_high_", "", scafs_posteriors)
scafs_posteriors_names <- gsub(".RData", "", scafs_posteriors_names)
names(scafs) <- scafs_posteriors_names

## get all intervals and areas 

all_intervals <- data.frame()
all_areas <- data.frame()
for (i in 1:length(scafs)){
  # get interval
  interval <- mcmc_intervals_data(scafs[[i]], prob =0.8, prob_outer = 0.95, pars = "b_scaletotal_load")
  interval$scaf_nr <- gsub("scaf_", "", names(scafs[i]))
  interval <- left_join(interval, genome[,c("scaf_nr", "n_snp")], by = "scaf_nr")
  
  # merge
  all_intervals <- rbind(all_intervals, interval)
  
  # get area
  area <- mcmc_areas_data(scafs[[i]], pars = "b_scaletotal_load")
  area$scaf_nr <- gsub("scaf_", "", names(scafs[i]))
  area <- left_join(area, genome[,c("scaf_nr", "n_snp")], by = "scaf_nr")
  
  
  # merge
  all_areas <- rbind(all_areas, area)
}


# split by interval
brms_plota <- split(all_areas, all_areas$interval)

brms_plota$bottom <- brms_plota$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

## plot
brms_plota$outer$scaf_nr <- as.numeric(brms_plota$outer$scaf_nr)
all_intervals$scaf_nr <- as.numeric(all_intervals$scaf_nr)

brms_plota$outer$n_snp <- as.numeric(brms_plota$outer$n_snp)
all_intervals$n_snp <- as.numeric(all_intervals$n_snp)

lm <- lm(m ~ scaf_nr, data = all_intervals)

all_intervals <- all_intervals %>% mutate(sig = case_when(
  m < 0 & l < 0 & h < 0 ~ "neg",
  m > 0 & l > 0 & h > 0 ~ "pos"
))

ggplot(data = all_intervals, aes(x = m, y = scaf_nr)) +  
  geom_segment(aes(x = l, xend = h, y = scaf_nr), col = "black", linewidth=3)+
  geom_segment(aes(x = ll, xend = hh, y = scaf_nr), col = "black")+
  geom_point(aes(fill = sig), col = "black", shape=21, size = 6) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = expression("Standardised"~beta), y = "Scaffold", title = "SnpEff")+
  # geom_smooth(method='lm')+
  scale_fill_manual(values =alpha(c(clr_highlight, clr_grey), 0.7)) +
  scale_color_manual(values =c(clr_highlight, clr_grey)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  coord_flip() -> high_load_per_chr

### save plots together ####
cowplot::plot_grid(gerp_load_per_chr, high_load_per_chr,
                   ncol = 1, labels = "auto", label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb")  -> loads_per_chr

ggsave(loads_per_chr, file = "plots/sup/extended_loads_per_chr.png", width=16, height = 12)
