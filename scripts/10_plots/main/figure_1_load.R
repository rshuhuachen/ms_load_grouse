#### Plots for MS #####
extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra)
source("scripts/theme_ggplot.R")

### Figure: Distribution of GERP load cat ####

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
 coord_flip() -> fig_countgerp

png(file = "plots/main/fig_1a.png", width=600, height=600)
fig_countgerp
dev.off()

#### Figure: number of SNPs in each category for snpeff ####

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
        text = element_text(size = 18)) + coord_flip() -> fig_countsnpef

fig_countsnpef
png(file = "plots/main/fig_1b.png", width=600, height=600)
fig_countsnpef
dev.off()

#### Figure: count of each SNPeff mutation variant ####

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
  scale_fill_manual(values = alpha(c(clr_high, "#8EA4CC","#8FAD88",  "#FFCD70"), 0.7))+
  scale_color_manual(values = c(clr_high, "#8EA4CC","#8FAD88",  "#FFCD70"))+
  coord_flip() +geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), 
                          hjust = 1.5, size = 6)+
  guides(col="none") +
  theme(legend.position = c(0.7,0.8)) -> fig_countsnpef_cat
fig_countsnpef_cat

png(file = "plots/main/fig_1c.png", width=600, height=800)
fig_countsnpef_cat
dev.off()

#### Figure: allele freq combining GERP and SnpEFf ####

# deleterious
load(file = "output/load/gerp/allelefreq_gerp5.derived.RData")
load(file = "output/load/snpeff/allelefreq_snpeff_high.derived.RData")
# 
# allele_freq_gerp_long_derived$Method <- "GERP score ≥ 4"
# allele_freq_high_long_derived$Method <- "High impact SnpEff"
# names(allele_freq_gerp_long_derived)[7] <- "allele_freq"
# names(allele_freq_high_long_derived)[7] <- "allele_freq"
# names(allele_freq_high_long_derived)[6] <- names(allele_freq_gerp_long_derived)[6]
# 
# allelefreq <- rbind(allele_freq_gerp_long_derived, allele_freq_high_long_derived)
# 
# ggplot(subset(allelefreq, Method == "GERP score ≥ 4"), aes(x = frequency, fill = Method, col = Method)) + 
#   scale_y_continuous(breaks=(c(50000, 100000, 150000)),
#                      labels=c("50 k", "100 k", "150 k"))+
#   geom_histogram(aes(y = after_stat(count)),
#                  binwidth=0.07, position="dodge") +
#   scale_color_manual(values = c(clr_gerp)) + 
#   scale_fill_manual(values = alpha(c(clr_gerp), 0.7))+
#   labs(x = "Allele frequency", y = "Frequency", title = "GERP ≥ 4")+
#   theme(legend.position = "none") -> fig_af_gerp
# 
# png(file = "plots/main/fig_1d.png", width=600, height=600)
# fig_af_gerp
# dev.off()
# 
# ggplot(subset(allelefreq, Method == "High impact SnpEff"), aes(x = frequency, fill = Method, col = Method)) + 
#   scale_y_continuous(breaks=(c(0, 1000, 2000)),
#                      labels=c("0", "1 k", "2 k"))+
#   geom_histogram(aes(y = after_stat(count)),
#                  binwidth=0.07, position="dodge") +
#   scale_color_manual(values = c(clr_high)) + 
#   scale_fill_manual(values = alpha(c(clr_high), 0.7))+
#   labs(x = "Allele frequency", y = "Frequency", title = "High impact SnpEff") +
#   theme(legend.position = "none")-> fig_af_high
# 
# png(file = "plots/main/fig_1e.png", width=600, height=600)
# fig_af_high
# dev.off()


### combine neutral and del in one fig ####
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
  theme(legend.position = c(0.8,0.9))-> fig_af_gerp
fig_af_gerp

png(file = "plots/main/fig_1d.png", width=600, height=600)
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
  theme(legend.position = c(0.8,0.9))-> fig_af_high

png(file = "plots/main/fig_1e.png", width=600, height=600)
fig_af_high
dev.off()

#### Figure: histogram number of loci hom het ######

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
png(file = "plots/main/fig_1g.png", width=600, height=600)
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
png(file = "plots/main/fig_1f.png", width=600, height=600)
hist_n_gerp
dev.off()

#### Combine in one figure ####
cowplot::plot_grid(fig_countgerp, fig_countsnpef,
                   ncol = 1, labels = c("a", "b"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> fig_load_1

cowplot::plot_grid(fig_load_1, fig_countsnpef_cat, rel_widths = c(0.4,0.6),
                   labels = c("", "c"), label_fontface = "plain", label_size = 22,
                   ncol = 2) -> fig_load_2

cowplot::plot_grid(fig_af_gerp, 
                   fig_af_high,
                   hist_n_gerp, 
                   hist_n_high, 
                   labels = c("d", "e", "f", "g"), label_fontface = "plain", label_size = 22,
                   ncol = 2,
                   align = "hv", axis = "lb") -> fig_load_3_b

cowplot::plot_grid(fig_load_2, fig_load_3_b, rel_heights = c(1, 0.8),
                   ncol = 1, 
                   align = "hv", axis = "lb") -> fig_load


png(file = "plots/main/fig_1_load.png", width=1000, height=1400)
fig_load
dev.off()


