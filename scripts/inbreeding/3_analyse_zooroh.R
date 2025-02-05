#### load packages ####
pacman::p_load(tidyverse, data.table, RZooRoH, lme4, glue)

### load theme
source("scripts/theme_ggplot.R")

### load rzooroh output ###
load(file = "output/inbreeding/zooroh_out_full2.RData")

### load headers
names <- fread("output/inbreeding/chicks_adults_samples.txt", header = F)
names(names) <- "file_id"

### load metadata
load("data/metadata/metadata_adult_chick.RData")

### load id info
#add real id
ids <- read.csv("/vol/cluster-data/rchen/git/inbreeding-sexual-traits/data/metadata/file_list_all_bgi_clean.csv")
ids$file_id <- gsub(".sorted.bam", "", ids$file)

names <- left_join(names, ids[,c("file_id", "id")], by = c("file_id"))

names <- names %>% mutate(id = case_when(
  is.na(id) ~ file_id,
  TRUE ~ id
))

names <- names %>% mutate(age = as.factor(case_when(
  grepl("C", names$id) ~ "chick",
  grepl("D", names$id) ~ "adult"
)))

names$n_id <- rownames(names)

save(names, file = "output/inbreeding/names_id.RData")

### divide out rzooroh ###
hbd <- realized(out)
hbd$file_id <- out@sampleids
hbd <- left_join(hbd, names[,c("n_id", "id", "age")], by =c("file_id" = "n_id"))
hbd <- hbd %>% relocate(id, .before=R_2)
hbd$file_id <- NULL

save(hbd, file = "output/inbreeding/hbd_froh.RData")

ggplot(hbd, aes(x = age, y = R_2))+ geom_boxplot(outlier.shape = NA) + geom_point(position="jitter") + 
  scale_y_log10()
summary(lm(R_2 ~ age, data = hbd)) #ns
ggplot(hbd, aes(x = age, y = R_20))+ geom_boxplot(outlier.shape = NA) + geom_point(position="jitter") + 
  scale_y_log10()
summary(lm(R_20 ~ age, data = hbd)) #na
ggplot(hbd, aes(x = age, y = R_200))+ geom_boxplot(outlier.shape = NA) + geom_point(position="jitter") + 
  scale_y_log10()
summary(lm(R_200 ~ age, data = hbd)) #ns
ggplot(hbd, aes(x = age, y = R_2000))+ geom_boxplot(outlier.shape = NA) + geom_point(position="jitter") + 
  scale_y_log10()
summary(lm(R_2000 ~ age, data = hbd)) #chick slightly lower (0.04)
ggplot(hbd, aes(x = age, y = 1-NonHBD))+ geom_boxplot(outlier.shape = NA) + geom_point(position="jitter") 
summary(lm(1-NonHBD ~ age, data = hbd)) #ns

hbd <- left_join(hbd, meta, by = "id")

summary(lmerTest::lmer(1-NonHBD ~ age + (1|site) + (1|year), data = hbd)) #ns

#compare to bcftools
load(file = "/vol/cluster-data/rchen/wgr/chicks_genomes/results/froh_29scaf.RData")

frohs <- left_join(hbd[,c("id", "NonHBD", "age")], froh_29scaf[,c("id", "froh")],by = "id")
ggplot(frohs, aes(1-NonHBD, froh)) + geom_point(aes(col = age), alpha = 0.6, size = 3) + geom_smooth(method='lm', col = "black") + 
  labs(y = "BCFtools", x = "ZooRoH") -> compare_tools
compare_tools
cor.test(1-frohs$NonHBD, frohs$froh) #0.93, p < 0.001

png("plots/inbreeding/compare_bcftools_rzooroh.png", width=600, height=600)
compare_tools
dev.off()

## ID in LMS?
load("data/metadata/phenotypes_wide_extra.RData")
pheno <- inner_join(pheno[,c("id", "core", "LMS_min")], hbd, by = "id")
pheno$froh <- 1-pheno$NonHBD

summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + core + (1|site), data = pheno, family = "poisson", ziformula = ~1)) # still sig negative!

summary(glmmTMB::glmmTMB(LMS_min ~ R_2 + R_20 + R_200 + R_2000 + core + (1|site), data = pheno, family = "poisson", ziformula = ~1)) # still sig negative!
summary(glmmTMB::glmmTMB(LMS_min ~ R_2 + R_4 + R_8 + R_16 + R_32 + R_64 + R_128 + R_256 + R_512 + core + (1|site), data = pheno, family = "poisson", ziformula = ~1)) # still sig negative!

## compare to derived homozygosity
load("output/inbreeding/homozygosity_derived_per_id.RData")

ggplot(load, aes(x = hom_anc_load)) + geom_histogram()
ggplot(load, aes(x = hom_der_load)) + geom_histogram()

# merge with age
load <- left_join(load, names[,c("id", "age")], by = "id")
load <- left_join(load, meta, by = "id")

ggplot(load, aes(x = age, y = hom_der_load)) + geom_boxplot() + geom_point()
summary(lm(hom_anc_load ~ age, data = load)) # chick have a higher ancestral homozygosity (p = 0.048)
summary(lm(hom_der_load ~ age, data = load)) # tendency that chicks have higher derived load (p = 0.07)
summary(lmerTest::lmer(hom_der_load ~ age + (1|site) + (1|year), data = load)) # tendency that chick have a higher derived homozygosity (p = 0.09)

#### plotting ####

### froh plots ####

## bar plot
long_hbd <- gather(hbd, HBDclass, F, R_2:NonHBD, factor_key=T)

long_hbd %>% subset(age == "chick" & HBDclass != "NonHBD") %>%
  ggplot() + 
  geom_col(aes(x = id, y = F, fill = gsub("R_", "", HBDclass)), col = "black",position="stack", width = 1) +
  theme(axis.text.x = element_blank()) + 
  labs(x = "ID", y = "F", fill = "HBDclass") + 
  scale_fill_manual(values = clrs_hunting[2:6]) -> froh_col_chick

froh_col_chick

png(file = "plots/inbreeding/froh_col_chick.png", width=800, height=600)
froh_col_chick
dev.off()

long_hbd %>% subset(age == "adult" & HBDclass != "NonHBD") %>%
  ggplot() + 
  geom_col(aes(x = id, y = F, fill = gsub("R_", "", HBDclass)), col = "black",position="stack", width = 1) +
  theme(axis.text.x = element_blank()) + 
  labs(x = "ID", y = "F", fill = "HBDclass") + 
  scale_fill_manual(values = clrs_hunting[2:5]) -> froh_col_adult 

froh_col_adult


png(file = "plots/inbreeding/froh_col_adult.png", width=800, height=600)
froh_col_adult
dev.off()

## line plot cum roh

long_hbd %>% subset(HBDclass != "NonHBD") %>% 
  ggplot(aes(x = log2(as.numeric(as.character(gsub("R_", "", HBDclass)))), y = F, group = id)) +
  geom_line(alpha = 0.3) +
  geom_point(alpha = 0.6, size = 3) +
  labs(x = "HBD class") -> cum_roh

cum_roh

png(file = "plots/inbreeding/cum_roh.png", width=600, height=400)
cum_roh
dev.off()

## boxplot prop in each class
long_hbd %>% subset(HBDclass != "NonHBD") %>% 
  ggplot(aes(x = gsub("R_", "", HBDclass), y = F, fill = HBDclass)) + 
  geom_point(position = "jitter", aes(col = HBDclass), size = 2)+
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  labs(x = "HBD class") +
  scale_fill_manual(values = clrs_hunting[2:5]) + 
  scale_color_manual(values = clrs_hunting[2:5]) + 
  theme(legend.position = "none") -> prop_hbd_class

prop_hbd_class

png(file = "plots/inbreeding/prop_hbd_class.png", width=600, height=400)
prop_hbd_class
dev.off()

### plot ROH segments
rohs <- out@hbdseg

# merge correct ids
hbd$id_nr <- out@ids
rohs <- left_join(rohs, hbd[,c("id_nr", "id")], by = c("id" = "id_nr"))  
rohs <- rohs %>% rename(id_nr = id, id = id.y)

# merge correct chrom
#data <- zoodata(genofile = "test.gen", zformat = "gp", chrcol = 1, poscol = 3, supcol = 5)
load("output/inbreeding/zooroh_chrnames.RData")

chrnames$chrom <- as.numeric(rownames(chrnames))
names(chrnames) <- c("scaf", "chrom")
rohs <- left_join(rohs, chrnames, by = "chrom")

rohs <- rohs %>% select(c(id, scaf, start_snp, end_snp, start_pos, end_pos, number_snp, length, HBDclass))
rohs <- rohs %>% mutate(hbdclass = as.factor(case_when(HBDclass == 1 ~ "2",
                                             HBDclass == 2 ~ "20",
                                             HBDclass == 3 ~ "200",
                                             HBDclass == 4 ~ "2000")))
chrms <- sort(unique(rohs$scaf))

# merge with age
rohs <- left_join(rohs, hbd[,c("id", "age")], by = 'id')
save(rohs, file = "output/inbreeding/rohs_zooroh.RData")

rohs %>%
  filter(scaf == "ScEsiA3_18278__HRSCAF_21663" & age == "chick")  %>% # filter(id == "ES3351")
  ggplot(aes(y = id)) +
  geom_linerange(aes(xmin = -Inf, xmax = Inf),
                 linewidth = .2) +
  geom_linerange(aes(xmin = start_pos, xmax = end_pos, color = hbdclass),
                 linewidth = 4) +
  labs(y = "ID", color = "HBD class") +
  scale_color_manual(values=clrs_hunting[2:5]) +
  scale_x_continuous(glue("Position on {chrms[[1]]}"),
                     labels = \(x){sprintf("%.0f Mb", x*1e-6)}) -> rohs_chicks
rohs_chicks

png(file = "plots/inbreeding/rohs_chicks.png", width=800, height=600)
rohs_chicks
dev.off()

rohs %>% filter(age == "adult") %>% select(id) %>% unique() %>% head(30) -> adult_30


rohs %>%
  filter(scaf == "ScEsiA3_18278__HRSCAF_21663" & id %in% adult_30$id)  %>% # filter(id == "ES3351")
  ggplot(aes(y = id)) +
  geom_linerange(aes(xmin = -Inf, xmax = Inf),
                 linewidth = .2) +
  geom_linerange(aes(xmin = start_pos, xmax = end_pos, color = hbdclass),
                 linewidth = 4) +
  labs(y = "ID", color = "HBD class") +
  scale_color_manual(values=clrs_hunting[2:5]) +
  scale_x_continuous(glue("Position on {chrms[[1]]}"),
                     labels = \(x){sprintf("%.0f Mb", x*1e-6)}) -> rohs_adult30
rohs_adult30

png(file = "plots/inbreeding/rohs_adult30.png", width=800, height=600)
rohs_adult30
dev.off()

### do chicks have a larger prop of their genome in shorter rohs? ####

short_rohs <- subset(rohs, length < 100000)
short_rohs_froh <- short_rohs %>% group_by(id) %>% summarise(froh = sum(length))
short_rohs_froh_n <- short_rohs %>% group_by(id) %>% count()
short_rohs_froh <- left_join(short_rohs_froh, unique(rohs[,c("id", "age")]), by = "id")
short_rohs_froh_n <- left_join(short_rohs_froh_n, unique(rohs[,c("id", "age")]), by = "id")
summary(lm(froh ~ age, data = short_rohs_froh)) #chicks have smaller froh of short rohs
summary(lm(n ~ age, data = short_rohs_froh_n)) #chicks have less short rohs

### per snp### per snpfrohs

#id 1
local_c06 <- as.data.frame(t(out@hbdp[[6]]))
names(local_c06) <- c("R_2", "R_20", "R_200", "R_2000", "nonHBD")
local_c06 <- local_c06 %>% mutate(pos = row_number(),
                                  id = hbd$id[2])

local_c06_long <- gather(local_c06, HBD, prob, R_2:nonHBD, factor_key = T) 

local_c06_long %>% subset(pos < 100000) %>%
  ggplot(aes(x = pos, y = prob, color = HBD)) + geom_line()+
  facet_grid(str_remove(HBD, "^R_") ~.) +
  scale_color_manual(values=clrs_hunting[2:6]) +
  theme(legend.position = "none") +
  ylim(0,1) + 
  labs(y = "Local HBD probabiliy", title = "ID C06") -> local_c06_plot

png(file = "plots/inbreeding/local_c06.png", width=800, height=600)
local_c06_plot
dev.off()


#id 2
local_d229 <- as.data.frame(t(out@hbdp[[100]]))
names(local_d229) <- c("R_2", "R_20", "R_200", "R_2000", "nonHBD")
local_d229 <- local_d229 %>% mutate(pos = row_number(),
                                  id = names[100]$id)

local_d229_long <- gather(local_d229, HBD, prob, R_2:nonHBD, factor_key = T) 

local_d229_long %>% subset(pos < 100000 ) %>%
  ggplot(aes(x = pos, y = prob, color = HBD)) + geom_line()+
  facet_grid(str_remove(HBD, "^R_") ~.) +
  scale_color_manual(values=clrs_hunting[2:6]) +
  theme(legend.position = "none") +
  ylim(0,1) + 
  labs(y = "Local HBD probabiliy", title = "ID D229032") -> local_d229_plot


png(file = "plots/inbreeding/local_d229.png", width=600, height=400)
local_d229_plot
dev.off()

### combine in cowplot ###
cowplot::plot_grid(cum_roh, prop_hbd_class, 
                   froh_col_chick, froh_col_adult, 
                   rohs_chicks, rohs_adult30, 
                   local_c06_plot, local_d229_plot, ncol = 2, labels = "auto",
                   rel_heights=c(1,1,2,2), align = "hv", axis = "lb") -> zooroh_plots

png(file = "plots/inbreeding/zooroh_out.png", width=1000, height=1600)
zooroh_plots
dev.off()

### Follow stoffel: froh contributions per size class to fitness ####
load(file = "output/inbreeding/hbd_froh.RData")
load(file = "output/inbreeding/rohs_zooroh.RData")

# mean roh per id
roh_length_per_id <- rohs %>% group_by(id) %>% summarise(mean_roh_length = mean(length, na.rm=T))
summary(roh_length_per_id$mean_roh_length)

#merge
roh_summary <- left_join(hbd, roh_length_per_id, by = "id")
#add pheno data
load("data/metadata/phenotypes_wide_extra.RData")
roh_summary <- inner_join(pheno[,c("id", "core", "LMS_min", "site")], roh_summary, by = "id")
roh_summary$froh <- 1-roh_summary$NonHBD

summary(glmmTMB::glmmTMB(LMS_min ~ scale(froh) + mean_roh_length + core + (1|site), data = roh_summary, family = "poisson", ziformula = ~1)) # still sig negative!
# only froh sig, mean roh length is not

summary(glmmTMB::glmmTMB(LMS_min ~ R_2 + R_20 + R_200+ R_2000 + core + (1|site), data = roh_summary, family = "poisson", ziformula = ~1)) # still sig negative!
#shorter rohs contribute more: R200 and R2000 most
