#### Plots for MS #####
extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra)
source("scripts/theme_ggplot.R")

#### average gerp score for each of the four categories of SnpEff ####
#### load all gerp scores ####
gerp_snp_scafs <- list.files(path = "output/gerp/beds", pattern = "gerp_overlapSNP*", full.names = T)
gerp_snp_scafs <- gerp_snp_scafs[-22] #empty

gerp_snp <- data.frame()
for (i in 1:length(gerp_snp_scafs)){
  scaf <- read.table(gerp_snp_scafs[i])
  gerp_snp <- rbind(gerp_snp, scaf)
}

head(gerp_snp)
## metadata on filenames and ids
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(gerp_snp) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns


### separate in groups of snpeff
high <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
mod <- subset(gerp_snp, grepl("MODERATE", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
low <- subset(gerp_snp, grepl("LOW", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
modify <- subset(gerp_snp, grepl("MODIFIER", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))

sum_gerp_per_snpeff <- data.frame(cat = c("high", "mod", "low", "modify"),
                                  mean = c(mean(high$rs_score), mean(mod$rs_score), mean(low$rs_score), mean(modify$rs_score)),
                                  median = c(median(high$rs_score), median(mod$rs_score), median(low$rs_score), median(modify$rs_score)),
                                  quant_25 = c(as.vector(quantile(high$rs_score, probs = 0.25)),
                                               as.vector(quantile(mod$rs_score, probs = 0.25)),
                                               as.vector(quantile(low$rs_score, probs = 0.25)),
                                               as.vector(quantile(modify$rs_score, probs = 0.25))),
                                  quant_75 = c(as.vector(quantile(high$rs_score, probs = 0.75)),
                                               as.vector(quantile(mod$rs_score, probs = 0.75)),
                                               as.vector(quantile(low$rs_score, probs = 0.75)),
                                               as.vector(quantile(modify$rs_score, probs = 0.75))))

write.csv(sum_gerp_per_snpeff, file = "output/load/sum_stats_gerp_per_snpeff_cat.csv", quote=F, row.names = F)

summary(high$rs_score) # mean = -0.57
summary(mod$rs_score) # mean = -0.02
summary(low$rs_score) # mean = -1.59
summary(modify$rs_score) # mean = -1.47


### relationship between things
# load mutation load measures
load("output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData") #loads no sex chr only 30 scaf

# subset only the relevant method/loadtype
load <- load_per_region %>% filter(loadtype == "gerp45"|loadtype=="high")

# rename
load$loadtype <- gsub("gerp45", "GERP ≥ 4", load$loadtype)
load$loadtype <- gsub("high", "High impact", load$loadtype)


#### Figure: inbreeding froh histogram #### 
# load froh mcmc
rohs <- fread("output/inbreeding/bcftools_roh_rg_clean.txt")
load("output/inbreeding/froh.RData")

ggplot(roh_autoscaf_bp_perid, aes(x = froh_auto)) + 
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  labs(x = expression(Inbreeding~coefficient~italic(F)[ROH]), y = "Number of individuals") -> fig_hist_froh

#### Figure: inbreeding roh number histogram #### 
num_roh_per_ind <- roh_autoscaf %>% group_by(id) %>% tally() 

ggplot(num_roh_per_ind, aes(n)) +
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  ylab("Number of individuals") +
  xlab("Number of ROHs")  -> fig_hist_nroh

#### Figure: distribution ROHs most inbred ####

### Distribution ROH most inbred (N = 5) and least inbred individuals (N = 5)
### ROH for most and least inbred individuals ###
all_roh <- roh_autoscaf %>% 
  group_by(id) %>% 
  summarise(sum_roh = sum(length)) %>% 
  ungroup() %>% 
  arrange(desc(sum_roh))

longest_roh <- all_roh %>% 
  top_n(5)
shortest_roh <- all_roh %>% 
  top_n(-5)
num_ind <- 10

extreme_roh <- rbind(longest_roh, shortest_roh)

df <- roh_autoscaf %>%
  mutate(POS1 = start / 1e+6,
         POS2 = end / 1e+6,
         MB = length / 1000)

df <- df %>% filter(id %in% extreme_roh$id) %>% 
  mutate(id = factor(id, levels = extreme_roh$id))

yax <- data.frame(id = forcats::fct_inorder(levels(df$id))) %>%
  mutate(yax = seq(from = 2,
                   to = 2*length(unique(df$id)),
                   by = 2)) 

df <- left_join(df, yax, by = "id")

#change from CHR to scaffold number
load(file="/vol/cluster-data/rchen/git/scaffold_names_dovetail.RData")
df <- left_join(df, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

#df$CHR <- as.numeric(as.factor(df$chr))
shade <- df %>%
  filter(scaf_nr < 11)%>%
  group_by(scaf_nr) %>%
  summarise(min = min(POS1), max = max(POS2)) %>%
  mutate(min = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10  ~ 0,
                         TRUE ~ min)) %>%
  mutate(max = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10  ~ 0,
                         TRUE ~ max))

col <- c("#008080", "#b4c8a8")

#take FROH measures from these individuals
yax <- left_join(yax, roh_autoscaf_bp_perid[,c("id", "froh_auto")])

df %>% 
  filter(MB > 1 & scaf_nr < 11) %>% 
  ggplot() +
  geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
            alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
  geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
  geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                fill = as.factor(scaf_nr)),  col = "grey", size = 0, alpha = 1) + 
  scale_y_reverse(expand = c(0, 0)) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
  scale_fill_manual(values = rep(clr_2, 28)) + 
  scale_color_manual(values = rep(clr_2, 28)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0, "lines"),
    # plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
    axis.line.x = element_blank(),
    legend.position="none",
    axis.title.x = element_text(margin=margin(t=10)),
    axis.title.y = element_text(margin=margin(r=1)),
    axis.text.y = element_text(colour = "white"),
    axis.line.y = element_blank()) +
  coord_cartesian(clip = 'off') +
  xlab("Scaffold") +
  ylab("Individuals") -> fig_roh_dist


# cor tests

homhet_gerp <- cor.test(load$hom_load[which(load$loadtype=="GERP ≥ 4")], load$het_load[which(load$loadtype=="GERP ≥ 4")]) # sig -0.94
homtot_gerp <- cor.test(load$hom_load[which(load$loadtype=="GERP ≥ 4")], load$total_load[which(load$loadtype=="GERP ≥ 4")]) #NS
hettot_gerp <- cor.test(load$het_load[which(load$loadtype=="GERP ≥ 4")], load$total_load[which(load$loadtype=="GERP ≥ 4")]) # sig 0.22

homhet_high <- cor.test(load$hom_load[which(load$loadtype=="High impact")], load$het_load[which(load$loadtype=="High impact")]) # sig -0.48
homtot_high <- cor.test(load$hom_load[which(load$loadtype=="High impact")], load$total_load[which(load$loadtype=="High impact")]) # sig 0.59
hettot_high <- cor.test(load$het_load[which(load$loadtype=="High impact")], load$total_load[which(load$loadtype=="High impact")]) # sig 0.42

cors <- data.frame(loadtype = c(rep("GERP ≥ 4", times = 3), rep("High impact", times = 3)),
                   x = c(rep(c("hom_load", "het_load", "hom_load"), times = 2)),
                   y = c(rep(c("total_load", "total_load", "het_load"), times = 2)),
                   pearson = round(c(homhet_gerp$estimate,
                                     homtot_gerp$estimate,
                                     hettot_gerp$estimate,
                                     homhet_high$estimate,
                                     homtot_high$estimate,
                                     hettot_high$estimate), 2),
                   pval = round(c(homhet_gerp$p.value,
                                  homtot_gerp$p.value,
                                  hettot_gerp$p.value,
                                  homhet_high$p.value,
                                  homtot_high$p.value,
                                  hettot_high$p.value), 3))

# plot relationships
ggplot(load, aes(x = hom_load, y = total_load)) + geom_point(aes(col = loadtype)) + geom_smooth(method='lm', col = "black")+
  facet_wrap(~loadtype, scales="free") +
  labs(x = "Homozygous load", y= "Total load")+
  scale_color_manual(values = c(clr_gerp, clr_high)) +
  theme(legend.position="none") -> tot_hom
tot_hom

ggplot(load, aes(x = het_load, y = total_load)) + geom_point(aes(col = loadtype)) + geom_smooth(method='lm', col = "black")+
  facet_wrap(~loadtype, scales="free") +
  labs(x = "Heterozygous load", y= "Total load")+
  scale_color_manual(values = c(clr_gerp, clr_high)) +
  theme(legend.position="none") -> tot_het

ggplot(load, aes(x = hom_load, y = het_load)) + geom_point(aes(col = loadtype)) + geom_smooth(method='lm', col = "black")+
  facet_wrap(~loadtype, scales="free") +
  labs(x = "Homozygous load", y= "Heterozygous load")+
  scale_color_manual(values = c(clr_gerp, clr_high)) +
  theme(legend.position="none") -> hom_het

plot_grid(tot_hom, tot_het, hom_het, 
          labels = "auto", label_fontface = "plain", label_size = 22,
          ncol = 1,
          align = "hv", axis = "lb") -> rels_load

ggsave(rels_load, file = "plots/main/fig_revision_relationships_load.png", width=10, height=12)

## distribution of loads 
long_load <- gather(load, zygosity, load, total_load:het_load, factor_key=T)
long_load$zygosity <- gsub("total_load", "Total", long_load$zygosity)
long_load$zygosity <- gsub("hom_load", "Hom", long_load$zygosity)
long_load$zygosity <- gsub("het_load", "Het", long_load$zygosity)

long_load$zygosity  <- factor(long_load$zygosity, levels=c("Total", "Hom", "Het"))

ggplot(long_load, aes(x = load, fill = zygosity)) + geom_histogram(position="identity", alpha=0.8) + 
  facet_wrap(~loadtype, scales="free")+ scale_fill_manual(values = c(clrs_hunting[1:3])) +
  labs(x = "Load", y = "Count", fill = "Load type") -> hist_loads

ggsave(hist_loads, file = "plots/main/fig_revision_hist_loads.png", width=12, height=8)

