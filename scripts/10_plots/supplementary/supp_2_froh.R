#### Comparing FROH plink vs bcftools #####
pacman::p_load(tidyverse)

source("scripts/theme_ggplot.R")

load(file="output/inbreeding/froh_combined.RData")

cor <-cor.test(combined$froh_bcftools, combined$froh_plink)


ggplot(combined, aes(froh_bcftools, froh_plink)) + 
  geom_point(col = alpha(clr_froh, 0.7), size=3) + 
  geom_smooth(method="lm", col = "black", linewidth=1) +
  annotate("text", label = paste0("r = ", round(cor$estimate, 2),
                                  ", p < 0.001"), 
           x = 0.20, y = 0.10, family = "Arial", size = 6)+
  labs(x = expression(F[ROH]~bcftools), y = expression(F[ROH]~plink)) -> froh

ggsave(froh, file="plots/Supp_2_froh_plink_bcf.png", width=8, height=8)
