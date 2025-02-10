pacman::p_load(tidyverse, data.table)

##### load in all gerps #####

gerp_scafs <- list.files(path = "output/gerp", pattern = "count_mutations*", full.names = T)

gerp_raw <- data.frame()
for (i in 1:length(gerp_scafs)){
  scaf <- fread(gerp_scafs[i])
  gerp_raw <- rbind(gerp_raw, scaf)
}

gerp <- gerp_raw %>% select(-c(scafno)) %>% group_by(id, gerp_cat) %>%
  summarise(across(n_total:hom_data, sum))

gerp <- gerp %>% mutate(het_load = het_data / n_genotyped,
                        hom_load = hom_data / n_genotyped,
                        total_load = (het_data * 0.5 + hom_data) / n_genotyped)

### explore relationship between different gerp classes ####
source("scripts/theme_ggplot.R")
wide <- spread(gerp[,c("id", "gerp_cat", "total_load")], gerp_cat, total_load)

ggplot(wide, aes(x = `4-5`, y = `< 0`)) + geom_point() + geom_smooth(method='lm') + 
  labs(x = "GERP load RS ≥ 4", y = "GERP load RS < 0") -> plot_gerp0_gerp45

ggsave(plot_gerp0_gerp45, file = "plots/sup/rel_gerp0_gerp45_load.png", width=8, height=8)

ggplot(wide, aes(x = `4-5`, y = `0-1`)) + geom_point()
ggplot(wide, aes(x = `4-5`, y = `1-2`)) + geom_point()
ggplot(wide, aes(x = `4-5`, y = `2-3`)) + geom_point()
ggplot(wide, aes(x = `4-5`, y = `3-4`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `3-4`, y = `2-3`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `1-2`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `2-3`, y = `1-2`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `2-3`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `2-3`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `2-3`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `1-2`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `1-2`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

load("data/phenotypes/phenotypes_lifetime.RData")

wide <- left_join(wide, pheno_wide[,c("id", "LMS_min", "site", "core")], by = "id")
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`4-5`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`3-4`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`2-3`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`1-2`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`0-1`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`< 0`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(`< 0`) + scale(`0-1`) + scale(`1-2`)+
                           scale(`2-3`) + scale(`3-4`)+ scale(`4-5`) + (1|site), data = wide, family = "poisson", ziformula = ~1))

## hom load
wide <- spread(gerp[,c("id", "gerp_cat", "hom_load")], gerp_cat, hom_load)

ggplot(wide, aes(x = `4-5`, y = `< 0`)) + geom_point() + geom_smooth(method='lm') + 
  labs(x = "GERP load RS ≥ 4", y = "GERP load RS < 0") -> plot_gerp0_gerp45

ggsave(plot_gerp0_gerp45, file = "plots/sup/rel_gerp0_gerp45_hom_load.png", width=8, height=8)

ggplot(wide, aes(x = `4-5`, y = `0-1`)) + geom_point()+ geom_smooth(method='lm')
ggplot(wide, aes(x = `4-5`, y = `1-2`)) + geom_point()+ geom_smooth(method='lm')
ggplot(wide, aes(x = `4-5`, y = `2-3`)) + geom_point()+ geom_smooth(method='lm')
ggplot(wide, aes(x = `4-5`, y = `3-4`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `3-4`, y = `2-3`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `1-2`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `3-4`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `2-3`, y = `1-2`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `2-3`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `2-3`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `2-3`, y = `0-1`)) + geom_point() + geom_smooth(method='lm')
ggplot(wide, aes(x = `1-2`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

ggplot(wide, aes(x = `1-2`, y = `< 0`)) + geom_point() + geom_smooth(method='lm')

load("data/phenotypes/phenotypes_lifetime.RData")

wide <- left_join(wide, pheno_wide[,c("id", "LMS_min", "site", "core")], by = "id")
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`4-5`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`3-4`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`2-3`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`1-2`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`0-1`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))
summary(glmmTMB::glmmTMB(LMS_min ~ scale(`< 0`) + core + (1|site), data = wide, family = "poisson", ziformula = ~1))

summary(glmmTMB::glmmTMB(LMS_min ~ scale(`< 0`) + scale(`0-1`) + scale(`1-2`)+
                           scale(`2-3`) + scale(`3-4`)+ scale(`4-5`) + (1|site), data = wide, family = "poisson", ziformula = ~1))

