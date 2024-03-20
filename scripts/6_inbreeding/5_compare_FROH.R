#### Compare FROH from different runs, sMLH and phenotypes ####

#load libraries
pacman::p_load(tidyverse,data.table,reshape2, lme4, performance)

#### Compare FROH to WGR sMLH #####

#load different froh runs
froh_bcf <- read.csv("output/inbreeding/bcftools_mcmc_froh.csv")
froh_plink <- read.csv("output/inbreeding/plink_froh.csv")
smlh <- read.csv("output/inbreeding/sMLH.csv")

#merge with smlh
combined <- left_join(froh_bcf, froh_plink, by = "id") %>%
  left_join(smlh, by = "id")

names(combined) <- c("id", 
                     "sum_bcftools", "froh_bcftools", 
                     "sum_plink", "froh_plink", 
                     "smlh")
combined[is.na(combined)] <- 0 #when FROH is NA, FROH = 0

save(combined, file = "output/inbreeding/froh_combined.RData")

#correlation tests
cor_bcf_plink <- cor.test(combined$froh_bcftools, combined$froh_plink)
cor_bcf_smlh <- cor.test(combined$froh_bcftools, combined$smlh)
cor_plink_smlh <- cor.test(combined$froh_plink, combined$smlh)

#plot relationship
ggplot(combined, aes(x = froh_bcftools, y = froh_plink))+ geom_point() + theme_classic() +
  geom_smooth(method = "lm") + labs(x = "FROH bcftools", y = "FROH plink", 
                                    subtitle=paste0("Pearson correlation coefficient = ", 
                                                    round(cor_bcf_plink$estimate,2), "
with p = ", round(cor_bcf_plink$p.value,2),"
N with ROH = ", length(combined$id[which(!is.na(combined$froh_bcftools))])))

ggplot(combined, aes(x = froh_bcftools, y = smlh))+ geom_point() + theme_classic() +
  geom_smooth(method = "lm") + labs(x = "FROH bcftools", y = "sMLH",
                                    subtitle=paste0("Pearson correlation coefficient = ", 
                                                    round(cor_bcf_smlh$estimate,2), "
with p = ", round(cor_bcf_smlh$p.value,2),"
N with ROH = ", length(combined$id[which(!is.na(combined$froh_bcftools))])))

ggplot(combined, aes(x = froh_plink, y = smlh))+ geom_point() + theme_classic() +
  geom_smooth(method = "lm") + labs(x = "FROH plink", y = "sMLH", 
                                    subtitle=paste0("Pearson correlation coefficient = ", 
                                                    round(cor_plink_smlh$estimate,2), "
with p = ", round(cor_plink_smlh$p.value,2),"
N with ROH = ", length(combined$id[which(!is.na(combined$froh_plink))])))
