### Effect of site on FROH and/or load? ###
## packages
pacman::p_load(tidyverse, data.table, lme4, lmerTest, rcartocolor, prismatic)

## data
load(file = "output/load/pheno_loads_lifetime.RData")


#### Site vs FROH1 ####

site_froh <- lm(scale(froh_bcftools) ~ site , data = pheno_wide_load)
ggplot(pheno_wide_load, aes(x = site, y = froh_bcftools)) + geom_point(alpha = 0.8, position = "jitter") + geom_boxplot(alpha = 0.5) +
  labs(subtitle = paste0("p-values glmer froh ~ site:
LEH = ", round(summary(site_froh)$coef[2,4], 3), ", NYR = ", round(summary(site_froh)$coef[3,4], 3),
                         ", SAA = ", round(summary(site_froh)$coef[4,4], 3), ", TEE = ", round(summary(site_froh)$coef[5,4], 3)),
       title = "Site vs FROH", x = "Site", y = "FROH") 

write.csv(as.data.frame(summary(site_froh)$coefficients), file = "output/inbreeding/site_vs_froh.csv", quote=F, col.names=T)

#### Site vs high impact load ####

site_high <- lm(scale(high_Lt_additive) ~ site , data = pheno_wide_load)
ggplot(pheno_wide_load, aes(x = site, y = high_Lt_additive)) + geom_point(alpha = 0.8, position = "jitter") + geom_boxplot(alpha = 0.5) +
  labs(subtitle = paste0("p-values glmer load ~ site :
LEH = ", round(summary(site_high)$coef[2,4], 3), ", NYR = ", round(summary(site_high)$coef[3,4], 3),
       ", SAA = ", round(summary(site_high)$coef[4,4], 3), ", TEE = ", round(summary(site_high)$coef[5,4], 3)),
       title = "Site vs total high impact load", x = "Site", y = "Total load") 

write.csv(as.data.frame(summary(site_high)$coefficients), file = "output/load/snpeff/site_vs_high.csv", quote=F, col.names=T)


#### Site vs gerp load ####

site_gerp <- lm(scale(gerp_Lt_additive_cat5) ~ site , data = pheno_wide_load)
ggplot(pheno_wide_load, aes(x = site, y = high_Lt_additive)) + geom_point(alpha = 0.8, position = "jitter") + geom_boxplot(alpha = 0.5) +
  labs(subtitle = paste0("p-values glmer load ~ site :
LEH = ", round(summary(site_gerp)$coef[2,4], 3), ", NYR = ", round(summary(site_gerp)$coef[3,4], 3),
       ", SAA = ", round(summary(site_gerp)$coef[4,4], 3), ", TEE = ", round(summary(site_gerp)$coef[5,4], 3)),
       title = "Site vs total high impact load", x = "Site", y = "Total load") 

write.csv(as.data.frame(summary(site_gerp)$coefficients), file = "output/load/gerp/site_vs_gerp5.csv", quote=F, col.names=T)
