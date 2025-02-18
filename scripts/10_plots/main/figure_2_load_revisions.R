extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra, patchwork)
source("scripts/theme_ggplot.R")

### Figure 2a: Distribution of GERP load cat ####

load(file = "output/load/gerp/gerps_count_per_cat.RData")
gerp_count$gerp_cat <- gsub("4-5", "≥ 4", gerp_count$gerp_cat)
gerp_count$gerp_cat <- factor(gerp_count$gerp_cat, levels = c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"))

ggplot(gerp_count, aes(x = gerp_cat, y = n_total)) + 
  geom_col(fill = alpha(c(clr_grey, clr_grey, clr_grey, clr_grey, clr_gerp, clr_grey), 0.7), 
           color=c(clr_grey,clr_grey, clr_grey, clr_grey,  clr_gerp, clr_grey)) + 
  # scale_y_log10(labels = c("0", expression(10^1),
  #                               expression(10^3), expression(10^5), expression(10^7)),
  #                    breaks = c(0, 10, 1000, 100000, 10000000)) +
  scale_y_log10(limits=c(1,10000000), labels = c(expression(10^1), expression(10^3), expression(10^6)),
                breaks=c(1, 1000,1000000)) +
  labs(x = "Category", y= expression('Number of SNPs (log'[10]*')'), title = "GERP") +
  geom_text(aes(label = prettyNum(n_total, big.mark=","), y = n_total), 
            hjust=1.5, size = 6) +
  coord_flip() + theme(plot.margin = margin(0.75,0,0.75,0.75, "cm"))-> fig_countgerp

## how many high impact per gerp score cat
# load all gerp scores
# gerp_snp_scafs <- list.files(path = "output/gerp/beds", pattern = "gerp_overlapSNP*", full.names = T)
# gerp_snp_scafs <- gerp_snp_scafs[-22] #empty
# 
# gerp_snp <- data.frame()
# for (i in 1:length(gerp_snp_scafs)){
#   scaf <- read.table(gerp_snp_scafs[i])
#   gerp_snp <- rbind(gerp_snp, scaf)
# }
# 
# head(gerp_snp)
# 
# ## metadata on filenames and ids
# filenames <- fread("data/genomic/raw/metadata/idnames.fam")
# ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
# 
# #merge
# idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))
# 
# names(gerp_snp) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns
# 
# 
# ### number of high impact in each cat
# high_gerp45 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score >=4)
# high_gerp34 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score >= 3 & rs_score< 4)
# high_gerp23 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score >= 2 & rs_score< 3)
# high_gerp12 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score >= 1 & rs_score< 2)
# high_gerp01 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score >= 0 & rs_score< 1) 
# high_gerp0 <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info) & rs_score < 0)
# 
# sum_high_per_gerp <- data.frame(gerp_cat = c("≥ 4", "3-4", "2-3", "1-2", "0-1", "< 0"),
#                                 n_high = c(nrow(high_gerp45), nrow(high_gerp34),
#                                            nrow(high_gerp23), nrow(high_gerp12),
#                                            nrow(high_gerp01), nrow(high_gerp0)))

# write.csv(sum_high_per_gerp, file = "output/load/sum_stats_high_per_gerp_cat.csv", quote=F, row.names = F)

sum_high_per_gerp <- read.csv(file = "output/load/sum_stats_high_per_gerp_cat.csv")
sum_high_per_gerp$gerp_cat <- factor(sum_high_per_gerp$gerp_cat, levels = c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"))
sum_high_per_gerp <- left_join(sum_high_per_gerp, gerp_count,by = "gerp_cat")
sum_high_per_gerp$n_perc <- round(sum_high_per_gerp$n_high / sum_high_per_gerp$n_total * 100, 2)

ggplot(sum_high_per_gerp, aes(x = n_perc, y = gerp_cat)) +
  geom_col(fill = alpha(clr_high, 0.7), col = clr_high) + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  xlim(0 ,0.25)+
  theme(plot.margin = margin(0.75, 0.75, 0.75, -1,unit= "cm"))+
  geom_text(aes(label = paste0(n_perc, " %"), y = gerp_cat), 
            hjust=-0.2, size = 6) +
  labs(x = "Percentage of high impact 
SnpEff SNPs", title = "") -> fig_count_high_per_gerp

fig_count_high_per_gerp

fig_countgerp+ fig_count_high_per_gerp-> fig_count_gerp_plus_high

png(file = "plots/main/fig_2a.png", width=1000, height=800)
fig_count_gerp_plus_high
dev.off()

### Figure 2b: number of SNPs in each category for snpeff ####

high_counts <- read.csv(file = "output/load/snpeff/n_mutations_per_type")

## only per impact category
count_high <- high_counts$n_mutations[which(high_counts$type == "total_high")]
count_moderate <- high_counts$n_mutations[which(high_counts$type == "total_mod")]
count_low <- high_counts$n_mutations[which(high_counts$type == "total_low")]
count_modify <- high_counts$n_mutations[which(high_counts$type == "total_modifier")]

## number of mutations per type
n_mutations_per_impact <- data.frame("type" = c("High", "Moderate", "Low", "Modifier"),
                                     "n_mutations" = c(count_high,
                                                       count_moderate,
                                                       count_low,
                                                       count_modify))

n_mutations_per_impact$type <- forcats::fct_relevel(n_mutations_per_impact$type, "Modifier", "Low", "Moderate", "High")

ggplot(n_mutations_per_impact, aes(x = type, y = n_mutations)) + 
  geom_col(color=c(clr_high, clr_grey, clr_grey, clr_grey), 
           fill=alpha(c(clr_high, clr_grey, clr_grey, clr_grey), 0.7)) + 
  scale_y_log10(labels = c(expression(paste(~10^1)), expression(paste(~10^4)), 
                           expression(paste(~10^7))),
                breaks=c(10, 10000, 10000000),
                limits = c(1, 10000000))+
  labs(y= expression('Number of SNPs (log'[10]*')'), title = "SnpEff", x = "Category")+
  scale_fill_manual(values = c(clrs[5], clrs[5], clrs[7], clrs[9]))+
  geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), hjust = 1.5, size =6 ) +
  theme(legend.position="none",
        plot.margin = margin(0.75,0,0.75,0.75, "cm"),
        text = element_text(size = 18)) + coord_flip() -> fig_countsnpef

## add: average gerp score for each of the four categories of SnpEff

### separate in groups of snpeff
# high <- subset(gerp_snp, grepl("HIGH", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
# mod <- subset(gerp_snp, grepl("MODERATE", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
# low <- subset(gerp_snp, grepl("LOW", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
# modify <- subset(gerp_snp, grepl("MODIFIER", gerp_snp$info) & !grepl("WARNING", gerp_snp$info))
# 
# sum_gerp_per_snpeff <- data.frame(cat = c("high", "mod", "low", "modify"),
#                                   mean = c(mean(high$rs_score), mean(mod$rs_score), mean(low$rs_score), mean(modify$rs_score)),
#                                   median = c(median(high$rs_score), median(mod$rs_score), median(low$rs_score), median(modify$rs_score)),
#                                   quant_25 = c(as.vector(quantile(high$rs_score, probs = 0.25)),
#                                                as.vector(quantile(mod$rs_score, probs = 0.25)),
#                                                as.vector(quantile(low$rs_score, probs = 0.25)),
#                                                as.vector(quantile(modify$rs_score, probs = 0.25))),
#                                   quant_75 = c(as.vector(quantile(high$rs_score, probs = 0.75)),
#                                                as.vector(quantile(mod$rs_score, probs = 0.75)),
#                                                as.vector(quantile(low$rs_score, probs = 0.75)),
#                                                as.vector(quantile(modify$rs_score, probs = 0.75))))
# 
# 
# ## number of mutations per type
# 
# write.csv(sum_gerp_per_snpeff, file = "output/load/sum_stats_gerp_per_snpeff_cat.csv", row.names = F, quote=F)

sum_gerp_per_snpeff <- read.csv("output/load/sum_stats_gerp_per_snpeff_cat.csv")
sum_gerp_per_snpeff$cat <- gsub("high", "High", sum_gerp_per_snpeff$cat)
sum_gerp_per_snpeff$cat <- gsub("modify", "Modifier", sum_gerp_per_snpeff$cat)
sum_gerp_per_snpeff$cat <- gsub("mod", "Moderate", sum_gerp_per_snpeff$cat)
sum_gerp_per_snpeff$cat <- gsub("low", "Low", sum_gerp_per_snpeff$cat)

sum_gerp_per_snpeff$mean <- round(sum_gerp_per_snpeff$mean, 1)
sum_gerp_per_snpeff$quant_25 <- round(sum_gerp_per_snpeff$quant_25, 1)
sum_gerp_per_snpeff$quant_75 <- round(sum_gerp_per_snpeff$quant_75, 1)

n_mutations_per_impact <- left_join(n_mutations_per_impact, sum_gerp_per_snpeff, by = c("type" = "cat"))
n_mutations_per_impact$type <- factor(n_mutations_per_impact$type, 
                                      levels=c("Modifier", "Low", "Moderate", "High"))

# plot

ggplot(n_mutations_per_impact, aes(x = median, y = type)) + 
  geom_segment(aes(x = quant_25, xend = quant_75, y = type), col = alpha(clr_grey, 0.7), size = 1) + 
  geom_point(size = 4, fill = "black", col = "black", shape=21) + 
  labs(x = "GERP score", title = "") +
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted") + 
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        plot.margin = margin(0.75,0.75,0.75,0, "cm")) -> fig_gerp_per_snpeff

fig_gerp_per_snpeff

fig_countsnpef+fig_gerp_per_snpeff-> fig_count_snpef_plus_gerp

png(file = "plots/main/fig_2b.png", width=1000, height=800)
fig_count_snpef_plus_gerp
dev.off()

### Figure 2c: count of each SNPeff mutation variant ####

n_mutations_pertype <- subset(high_counts, !is.na(impact))
n_mutations_pertype$impact <- gsub("low", "Low", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("moderate", "Moderate", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("high", "High", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("modify", "Modifier", n_mutations_pertype$impact)
n_mutations_pertype$impact <- factor(n_mutations_pertype$impact, 
                                     levels = c("High", "Moderate", "Low", "Modifier"))

n_mutations_pertype$abb <- c("Downstream gene variant",
                             "5' UTR premature start codon",
                             "5' UTR variant",
                             "Initiator codon variant",
                             "Intron variant",
                             "Loss of Function",
                             "Missense variant",
                             "Nonsense mediated decay",
                             "Splice acceptor variant",
                             "Splice donor variant",
                             "Splice region variant",
                             "Start codon lost",
                             "Stop codon gained",
                             "Stop codon lost",
                             "Stop codon retained",
                             "Synonymous variant",
                             "3' UTR varaint",
                             "Upstream gene variant",
                             "Intergenic region")


ggplot(n_mutations_pertype, aes(x = reorder(abb, desc(n_mutations)), 
                                y = n_mutations)) + 
  geom_col(aes(fill = impact, col = impact)) + 
  scale_y_log10(labels = c(expression(paste(~10^1)), expression(paste(~10^3)), expression(paste(~10^5))), 
                breaks = c(10, 1000, 100000))+
  labs(y = expression('Number of SNPs (log'[10]*')'), fill = "Impact Class", x = "Mutation type")+
  scale_fill_manual(values = alpha(c(clr_high, "#8EA4CC","#703D57",  "#FFCD70"), 0.7))+
  scale_color_manual(values = c(clr_high, "#8EA4CC","#703D57",  "#FFCD70"))+
  coord_flip() +geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), 
                          hjust = 1.5, size = 6)+
  guides(col="none") +
  theme(legend.position = c(0.7,0.8)) -> fig_countsnpef_cat
fig_countsnpef_cat

png(file = "plots/main/fig_2c.png", width=600, height=800)
fig_countsnpef_cat
dev.off()

### Figure 2d and 2e: allele freq combining GERP and SnpEFf ####

# deleterious
load(file = "output/load/gerp/allelefreq_gerp5.derived.RData")
load(file = "output/load/snpeff/allelefreq_snpeff_high.derived.RData")

# neutral
load(file = "output/load/snpeff/allelefreq_snpeff_low.derived.RData") #allele_freq_low_long_derived
load(file = "output/load/gerp/allelefreq_gerp0.derived.RData") #allele_freq_gerp0_long_derived

allele_freq_gerp_long_derived$Method <- "GERP"
allele_freq_high_long_derived$Method <- "SnpEff"
allele_freq_gerp0_long_derived$Method <- "GERP"
allele_freq_low_long_derived$Method <- "SnpEff"

allele_freq_gerp_long_derived$Deleteriousness <- "GERP ≥ 4"
allele_freq_high_long_derived$Deleteriousness <- "High impact"
allele_freq_gerp0_long_derived$Deleteriousness <- "GERP < 0"
allele_freq_low_long_derived$Deleteriousness <- "Low impact"

names(allele_freq_gerp_long_derived)[7] <- "allele_freq"
names(allele_freq_high_long_derived)[7] <- "allele_freq"
names(allele_freq_gerp0_long_derived)[7] <- "allele_freq"
names(allele_freq_low_long_derived)[7] <- "allele_freq"
names(allele_freq_high_long_derived)[6] <- names(allele_freq_gerp_long_derived)[6]
names(allele_freq_low_long_derived)[6] <- names(allele_freq_gerp_long_derived)[6]

allelefreq <- rbind(allele_freq_gerp_long_derived, allele_freq_high_long_derived,
                    allele_freq_gerp0_long_derived, allele_freq_low_long_derived)

ggplot(subset(allelefreq, Method == "GERP"), aes(x = frequency, fill = Deleteriousness, col = Deleteriousness)) + 
  # scale_y_continuous(breaks=(c(50000, 100000, 150000)),
  #                    labels=c("50 k", "100 k", "150 k"))+
  geom_histogram(aes(y=after_stat(c(count[group==1]/sum(count[group==1]),
                                    count[group==2]/sum(count[group==2]))*100)),
                 binwidth=0.07, position="dodge") +
  scale_color_manual(values = c("grey", clr_gerp)) + 
  scale_fill_manual(values = alpha(c("grey", clr_gerp), 0.7))+
  labs(x = "Allele frequency", y = "Percentage of SNPs", title = "GERP") +
  theme(legend.position = c(0.8,0.9),
        legend.title = element_blank())-> fig_af_gerp
fig_af_gerp

png(file = "plots/main/fig_2d.png", width=600, height=600)
fig_af_gerp
dev.off()

snpeff_af <- subset(allelefreq, Method == "SnpEff")
snpeff_af$Deleteriousness <- factor(snpeff_af$Deleteriousness, levels = c("Low impact", "High impact"))

ggplot(snpeff_af, aes(x = frequency, fill = Deleteriousness, col = Deleteriousness)) + 
  geom_histogram(aes(y=after_stat(c(count[group==1]/sum(count[group==1]),
                                    count[group==2]/sum(count[group==2]))*100)),
                 binwidth=0.07, position="dodge") +
  scale_color_manual(values = c("grey", clr_high)) + 
  scale_fill_manual(values = alpha(c("grey", clr_high), 0.7))+
  labs(x = "Allele frequency", y = "Frequency", title = "SnpEff") +
  theme(legend.position = c(0.8,0.9),
        legend.title = element_blank())-> fig_af_high

png(file = "plots/main/fig_2e.png", width=600, height=600)
fig_af_high
dev.off()

### Figure 2f + 2g: histogram number of loci hom het ######

# high
high <- read.table("data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf.gz")
source("scripts/7_calculate_load/function_calculate_load.R")
high_load <- calculate_load_snpeff(high, output_vcf = TRUE, loadtype = "high")
high <- high_load$load
high <- high %>% mutate(n_total_id_Heterozygosity = het_load*n_genotyped)
high <- high %>% mutate(n_total_id_Homozygosity = hom_load*n_genotyped)
high <- high %>% mutate(n_total_mutations = n_total_id_Heterozygosity+(n_total_id_Homozygosity*2))
sd(high$n_total_mutations)
mean(high$n_total_mutations)

high_loci <- gather(high[,c("id", "n_total_id_Heterozygosity", "n_total_id_Homozygosity")], type, n_loci, c(n_total_id_Heterozygosity:n_total_id_Homozygosity), factor_key = T)
high_loci$type <- gsub("n_total_id_", "", high_loci$type)
high_loci$type <- gsub("gosity", "gous", high_loci$type)
high_loci$type <- factor(high_loci$type, levels = c("Homozygous", "Heterozygous"))

ggplot(high_loci, aes(x = n_loci)) +
  geom_histogram(bins = 30, position = "identity", aes(fill = type, col = type)) +
  labs(x = "Number of loci", fill = "", y = "Number of individuals", title = "High impact SnpEff") +
  theme(legend.position = "top") + guides(color="none")+
  scale_fill_manual(values = c(alpha(clr_high, 0.8), alpha(clr_high, 0.4)))+
  scale_color_manual(values = c(clr_high, alpha(clr_high, 0.6)))-> hist_n_high

hist_n_high
png(file = "plots/main/fig_2g.png", width=600, height=600)
hist_n_high
dev.off()

# gerp
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

gerp <- subset(gerp, gerp_cat == "4-5")

gerp_loci <- gather(gerp[,c("id", "het_data", "hom_data")], type, n_loci, c(het_data:hom_data), factor_key = T)
gerp_loci$type <- gsub("het_data", "Heterozygous", gerp_loci$type)
gerp_loci$type <- gsub("hom_data", "Homozygous", gerp_loci$type)
gerp_loci$type <- factor(gerp_loci$type, levels = c("Homozygous", "Heterozygous"))

gerp <- gerp %>% mutate(n_total_mutations = het_data+(hom_data*2))
sd(gerp$n_total_mutations)
mean(gerp$n_total_mutations)

ggplot(gerp_loci, aes(x = n_loci)) +
  geom_histogram(bins = 30, position = "identity", aes(fill = type, col = type)) +
  labs(x = "Number of loci", fill = "", y = "Number of individuals", title = "GERP ≥ 4") +
  theme(legend.position = "top") + guides(color="none")+
  scale_x_continuous(labels = c("30,000", "40,000", "50,000"), breaks=c(30000,40000,50000))+
  scale_fill_manual(values = c(alpha(clr_gerp, 0.8), alpha(clr_gerp, 0.4)))+
  scale_color_manual(values = c(clr_gerp, alpha(clr_gerp, 0.6))) -> hist_n_gerp

hist_n_gerp
png(file = "plots/main/fig_2f.png", width=600, height=600)
hist_n_gerp
dev.off()

### Combine in one figure ####
cowplot::plot_grid(fig_count_gerp_plus_high, fig_count_snpef_plus_gerp,
                   ncol = 1, labels = c("a", "b"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> fig_mut_1

cowplot::plot_grid(fig_mut_1, fig_countsnpef_cat, 
                   labels = c("", "c"), label_fontface = "plain", label_size = 22,
                   ncol = 2) -> fig_mut_2

cowplot::plot_grid(fig_af_gerp, 
                   fig_af_high,
                   hist_n_gerp, 
                   hist_n_high, 
                   labels = c("d", "e", "f", "g"), label_fontface = "plain", label_size = 22,
                   ncol = 4,
                   align = "hv", axis = "lb") -> fig_mut_3

cowplot::plot_grid(fig_mut_2, fig_mut_3, rel_heights = c(1, 0.5),
                   ncol = 1, 
                   align = "hv", axis = "lb") -> fig_mut


png(file = "plots/main/fig_2.png", width=1500, height=1200)
fig_mut
dev.off()

