### load packages ####
pacman::p_load(tidyverse)

source("scripts/theme_ggplot.R")

#### Number of SNPs per region ####
load(file = "output/load/snps_by_region_plot.RData")

ggplot(subset(snps_by_region, method == "GERP scores >= 4" & region != "Non-coding"), aes(y = region, x = n_snps)) + 
  geom_col(position="dodge", fill = clr_gerp) +
  labs(x = "Number of SNPs", y = "Region", title = "GERP score 4-5")+
  geom_text(aes(label = paste0("N = ", prettyNum(round(n_snps, 2), big.mark=",")), y = region, hjust=-0.3),
            position=position_dodge(width=0.9), size = 6) + 
  scale_x_continuous(labels = c("0", "50k", "100k", "150k", "200k"), 
                     breaks = c(0,50000,100000,150000, 200000), limits = c(0,250000))+
  theme(legend.position = "bottom",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> snps_per_region_plot_gerp

ggplot(subset(snps_by_region, method == "High impact SnpEff"& region != "Non-coding"), aes(y = region, x = n_snps)) + 
  geom_col(position="dodge", fill = clr_high) +
  labs(x = "Number of SNPs", y = "Region", title = "High impact SnpEff")+
  geom_text(aes(label = paste0("N = ", prettyNum(round(n_snps, 2), big.mark=",")), 
                y = region, hjust=-0.3),
            position=position_dodge(width=0.9), size = 6) + 
  scale_x_continuous(labels = c("0", "2,500", "5,000", "7,500", "10,000"), 
                     breaks = c(0,2500,5000,7500,10000), limits = c(0,10000))+
  theme(legend.position = "bottom",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) -> snps_per_region_plot_high

cowplot::plot_grid(snps_per_region_plot_gerp, snps_per_region_plot_high,
                   ncol = 1, 
                   align = "hv", axis = "lb", 
                   labels = "auto", label_fontface = "plain", label_size = 22) -> plota_region

ggsave(plota_region, file="plots/Supp_6_snps_per_region.png", height = 12, width = 16)

