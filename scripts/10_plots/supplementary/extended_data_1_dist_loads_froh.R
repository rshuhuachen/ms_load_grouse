extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra)
source("scripts/theme_ggplot.R")

### Load data ####
# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load <- load_per_region %>% filter(loadtype == "gerp45"|loadtype=="high")

# rename
load$loadtype <- gsub("gerp45", "GERP ≥ 4", load$loadtype)
load$loadtype <- gsub("high", "High impact SnpEff", load$loadtype)

# load inbreeding
load("output/inbreeding/froh.RData")

#### Fig 3a and b - distribution of loads  ####
long_load <- gather(load, zygosity, load, total_load:het_load, factor_key=T)
long_load$zygosity <- gsub("total_load", "Total", long_load$zygosity)
long_load$zygosity <- gsub("hom_load", "Homozygous", long_load$zygosity)
long_load$zygosity <- gsub("het_load", "Heterozygous", long_load$zygosity)

long_load$zygosity  <- factor(long_load$zygosity, levels=c("Total", "Homozygous", "Heterozygous"))

ggplot(subset(long_load, loadtype == "GERP ≥ 4"), aes(x = load, fill = zygosity, col = zygosity)) + 
  geom_histogram(position="identity", alpha=0.8, bins = 40) + 
  facet_wrap( ~zygosity, scales="free")+
  scale_fill_manual(values = alpha(c(clrs_hunting[1:3]), 0.8)) +
  scale_color_manual(values = c(clrs_hunting[1:3])) +
  theme(strip.background = element_blank(),
        legend.position="none")+
  guides(col = "none")+
  labs(x = "Load", y = "Count", fill = "Load type", title = "GERP ≥ 4") -> hist_loads_gerp

hist_loads_gerp

png(file = "plots/sup/extended_data_1a.png", width=800, height=600)
hist_loads_gerp
dev.off()

ggplot(subset(long_load, loadtype == "High impact SnpEff"), aes(x = load, fill = zygosity, col = zygosity)) + 
  geom_histogram(position="identity", alpha=0.8, bins = 40) + 
  facet_wrap(~zygosity, scales="free")+ 
  scale_fill_manual(values = alpha(c(clrs_hunting[1:3]), 0.8)) +
  scale_color_manual(values = c(clrs_hunting[1:3])) +
  guides(col = "none")+
  theme(strip.background = element_blank(),
        legend.position="none")+
  labs(x = "Load", y = "Count", fill = "Load type", title = "High impact SnpEff") -> hist_loads_high

hist_loads_high

png(file = "plots/sup/extended_data_1b.png", width=800, height=600)
hist_loads_high
dev.off()

### Fig 3c - relationship between inbreeding and loads ####

# merge
load_froh <- left_join(load, froh, by = "id")

# make long
load_froh_long <- gather(load_froh, zygosity, load, het_load:total_load, factor_key=T)

load_froh_long$zygosity <- gsub("hom_load", "Homozygous", load_froh_long$zygosity)
load_froh_long$zygosity <- gsub("het_load", "Heterozygous", load_froh_long$zygosity)
load_froh_long$zygosity <- gsub("total_load", "Total", load_froh_long$zygosity)

load_froh_long$zygosity  <- factor(load_froh_long$zygosity, levels=c("Total", "Homozygous", "Heterozygous"))

ggplot(load_froh_long, aes(x = froh, y = load, col = zygosity, fill = zygosity)) + 
  geom_point(size=4, shape=21) + geom_smooth(method='lm', fill = clr_grey, show.legend = F)+
  facet_wrap(~loadtype, scales="free") +
  theme(strip.background = element_blank(),
        legend.position="bottom")+
  labs(x = expression(italic(F)[ROH]), y= "Load", fill = "Load type") +
  scale_fill_manual(values = alpha(c(clrs_hunting[1:3]), 0.6)) +
  guides(col = "none")+
  scale_color_manual(values = c(clrs_hunting[1:3]))  -> froh_load 

froh_load

png(file = "plots/sup/extended_data_1c.png", width=800, height=600)
froh_load
dev.off()

### Cor tests ####

# hom and froh
cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Homozygous")],
         load_froh_long$load[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Homozygous")])

cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Homozygous")],
         load_froh_long$load[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Homozygous")])

# het and froh
cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Heterozygous")],
         load_froh_long$load[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Heterozygous")])

cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Heterozygous")],
         load_froh_long$load[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Heterozygous")])

# total and froh
cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Total")],
         load_froh_long$load[which(load_froh_long$loadtype == "GERP ≥ 4" & load_froh_long$zygosity == "Total")])

cor.test(load_froh_long$froh[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Total")],
         load_froh_long$load[which(load_froh_long$loadtype == "High impact SnpEff" & load_froh_long$zygosity == "Total")])

### Combine in one figure ####
cowplot::plot_grid(hist_loads_gerp, hist_loads_high, froh_load,
                   ncol = 1, labels = c("a", "b", "c"), label_fontface = "plain", label_size = 22,
                   rel_heights = c(0.8, 0.8, 1),
                   align = "hv", axis = "lb") -> fig_loads


png(file = "plots/sup/extended_data_1.png", width=800, height=1000)
fig_loads
dev.off()

