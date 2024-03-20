#### Plots for MS #####

pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra)
source("scripts/theme_ggplot.R")

### Figure: Distribution of GERP load cat ####

load(file = "output/load/gerp/gerps_count_per_cat.RData")

ggplot(gerps_select_count, aes(x = gerp_cat, y = n)) + 
  geom_col(fill = alpha(c(clr_grey, clr_grey, clr_grey, clr_grey, clr_grey, clr_gerp), 0.7), 
           color=c(clr_grey, clr_grey, clr_grey, clr_grey, clr_grey, clr_gerp)) + 
  scale_y_continuous(labels = c("0", expression(paste("1 x"~10^6)), 
                                expression(paste("2 x"~10^6)), expression(paste("3 x"~10^6)), expression(paste("4 x"~10^6))), 
                     breaks = c(0, 1000000, 2000000, 3000000, 4000000), 
                     limits=c(0, 4500000)) +
  labs(x = "Category", y= "Number of SNPs", title = "GERP") +
  geom_text(aes(label = prettyNum(n, big.mark=","), y = n), 
            hjust=-0.2, size = 6) +
 coord_flip() -> fig_countgerp

#### Figure: allele freq combining GERP and SnpEFf ####
load(file = "output/load/gerp/allelefreq_gerp5.derived.RData")
load(file = "output/load/snpeff/allelefreq_snpeff_high.derived.RData")

allele_freq_gerp_long_derived$Method <- "GERP score 4-5"
allele_freq_high_long_derived$Method <- "High impact SnpEff"
names(allele_freq_gerp_long_derived)[6] <- "allle_freq"
names(allele_freq_high_long_derived)[6] <- "allle_freq"

allelefreq <- rbind(allele_freq_gerp_long_derived, allele_freq_high_long_derived)

ggplot(allelefreq, aes(x = frequency, fill = Method, col = Method)) + 
  geom_histogram(binwidth=0.07, position="dodge", mapping = aes(y = after_stat(ncount))) + 
  scale_color_manual(values = c(clr_gerp, clr_high)) + 
  scale_fill_manual(values = alpha(c(clr_gerp, clr_high), 0.7))+
  labs(x = "Allele frequency", y = "Proportion") +
  theme(legend.position = c(0.7, 0.8),
        legend.title = element_blank()) -> fig_hist_af


#### Figure: number of SNPs in each category for snpeff ####

## only per impact category
#grep -v "#" data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf |  wc -l
count_high <- 5289
#grep -v "#" data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf |  wc -l
count_moderate <- 137408
#grep -v "#" data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf |  wc -l
count_low <- 118930

## number of mutations per type
n_mutations_per_impact <- data.frame("type" = c("High", "Moderate", "Low"),
                                     "n_mutations" = c(count_high,
                                                       count_moderate,
                                                       count_low))

n_mutations_per_impact$type <- forcats::fct_relevel(n_mutations_per_impact$type, "Low", "Moderate", "High")

ggplot(n_mutations_per_impact, aes(x = type, y = n_mutations)) + 
  geom_col(color=c(clr_high, clr_grey, clr_grey), 
           fill=alpha(c(clr_high, clr_grey, clr_grey), 0.7)) + 
  scale_y_continuous(labels = c("0", "50,000", "100,000", "150,000"), breaks = c(0,50000, 100000, 150000), 
                     limits = c(0,180000))+
  labs(y= "Number of SNPs", title = "SnpEff", x = "Category")+
  scale_fill_manual(values = c(clrs[5], clrs[7], clrs[9]))+
  geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), hjust = -0.2, size =6 ) +
  theme(legend.position="none",
        text = element_text(size = 18)) + coord_flip() -> fig_countsnpef

#### Figure: distribution load z trans ####

#### data ###
load("output/load/all_loads_combined_da_nosex_29scaf.RData") #loads no sex chr only 29 scaf

select_load <- loads %>% select(c(id, gerp_Lt_additive_cat4, gerp_Lt_additive_cat5, high_Lt_additive))
names(select_load) <- c("id", "GERP score 3-4", "GERP score 4-5", "High impact SnpEff")

#z-transform
select_load_z <- select_load
select_load_z$`GERP score 4-5` <- scale(select_load_z$`GERP score 4-5`)
select_load_z$`High impact SnpEff` <- scale(select_load_z$`High impact SnpEff`)

select_load_long <- gather(select_load_z, load_type, load, "GERP score 4-5":"High impact SnpEff", factor_key=T )

ggplot(select_load_long, aes(x = load, fill = load_type, col = load_type)) + geom_histogram(position="identity") + 
  labs(x = "z-transformed total genetic load", y = "Number of individuals", fill = "Method") + 
  scale_fill_manual(values = alpha(c(clr_gerp, clr_high), 0.7))+
  scale_color_manual(values = c(clr_gerp, clr_high))+
  guides(col = "none")+
  theme(legend.position = "bottom",
        legend.title = element_blank()) -> fig_distribution

#### Figure: relationship loads ####

# correlation test
cor_gerp_high <- cor.test(loads$gerp_Lt_additive_cat5, loads$high_Lt_additive)

ggplot(loads, aes(x=gerp_Lt_additive_cat5, y=high_Lt_additive)) + 
  geom_point(size=2, color = "black", fill = clr_grey) + 
  geom_smooth(method="lm", col = clr_grey) + 
  ylim(0.145, 0.161)+
  labs(x =  expression("Total load GERP">="4"), y = "Total high impact load") + 
  annotate("text", label = paste0("r = ", round(cor_gerp_high$estimate, 2),
                           ", p = ", round(cor_gerp_high$p.value, 2)), 
            x = 0.151, y = 0.160, family = "Arial", size = 6) -> fig_correlation

#### Figure: overlap heatmap ####
load(file = "output/load/overlap_mutations.RData")

ggplot(overlaps, aes(x = gerp, y = snpef, fill = perc)) +
  geom_tile() +
  labs(x = "GERP score category", y ="SnpEff category", fill = "Percentage\noverlap") +
  scale_y_discrete(labels=scales::label_wrap(10))+
  scale_fill_gradient(low = clrs_hunting[5], high=clrs_hunting[4]) -> fig_overlap_heat

#### Figure: count of each SNPeff mutation variant ####

n_mutations_pertype <- read.csv("output/load/snpeff/n_mutations_per_type.csv")
n_mutations_pertype <- subset(n_mutations_pertype, !is.na(n_mutations))
n_mutations_pertype$impact <- gsub("low", "Low", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("moderate", "Moderate", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("high", "High", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("modify", "Modifier", n_mutations_pertype$impact)
n_mutations_pertype$impact <- factor(n_mutations_pertype$impact, 
                                     levels = c("High", "Moderate", "Low", "Modifier"))

n_mutations_pertype$abb <- c("Downstream gene variant",
                             "5' UTR premature start codon",
                             "Initiator codon variant",
                             "Loss of Function",
                             "Missense",
                             "Nonsense mediated decay",
                             "Splice region variant",
                             "Start codon lost",
                             "Stop codon gained",
                             "Stop codon lost",
                             "Stop codon retained",
                             "Synonymous variant")

ggplot(n_mutations_pertype, aes(x = reorder(abb, desc(n_mutations)), 
                                  y = n_mutations)) + 
  geom_col(aes(fill = impact, col = impact)) + 
  scale_y_log10(labels = c("10", "1000", "100,000"), breaks = c(10, 1000, 100000))+
  labs(y = "Number of SNPs", fill = "Impact Class", x = "Mutation type")+
  scale_fill_manual(values = alpha(c(clr_high, "#549F93", "#8EA4CC", clrs_related[1]), 0.7))+
  scale_color_manual(values = c(clr_high, "#549F93", "#8EA4CC", clrs_related[1]))+
  coord_flip() +geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), 
                          hjust = 1.5, size = 6)+
  guides(col="none") +
  theme(legend.position = c(0.8,0.8)) -> fig_countsnpef_cat

#### Figure: inbreeding froh histogram #### 
# load froh mcmc
rohs <- fread("output/inbreeding/mcmc_roh_rg_filtered.txt")
froh_mcmc <- fread("output/inbreeding/bcftools_mcmc_froh.csv")

ggplot(froh_mcmc, aes(x = froh_mcmc)) + 
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  labs(x = expression(Inbreeding~coefficient~F[ROH]), y = "Number of individuals") -> fig_hist_froh

#### Figure: inbreeding roh number histogram #### 
num_roh_per_ind <- rohs %>% group_by(id) %>% tally() 

ggplot(num_roh_per_ind, aes(n)) +
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  ylab("Number of individuals") +
  xlab("Number of ROHs")  -> fig_hist_nroh

#### Figure: distribution ROHs most inbred ####

### Distribution ROH most inbred (N = 5) and least inbred individuals (N = 5)
### ROH for most and least inbred individuals ###
all_roh <- rohs %>% 
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

df <- rohs %>%
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
load(file="data/genomic/raw/metadata/scaffold_names_dovetail.RData")
df <- left_join(df, genome[,c("contig", "scaf_nr")], by = c("chr" = "contig"))

#df$CHR <- as.numeric(as.factor(df$chr))
shade <- df %>%
  filter(scaf_nr < 10)%>%
  group_by(scaf_nr) %>%
  summarise(min = min(POS1), max = max(POS2)) %>%
  mutate(min = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10  ~ 0,
                         TRUE ~ min)) %>%
  mutate(max = case_when(scaf_nr == 2 | scaf_nr == 4 | scaf_nr == 6 | scaf_nr == 8 | scaf_nr == 10  ~ 0,
                         TRUE ~ max))

col <- c("#008080", "#b4c8a8")

#take FROH measures from these individuals
yax <- left_join(yax, froh_mcmc[,c("id", "froh_mcmc")])

df %>% 
  filter(MB > 1 & scaf_nr < 10) %>% 
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


#### Combine in one figure ####

cowplot::plot_grid(fig_hist_froh, fig_hist_nroh, 
              labels = c("b", "c"), label_fontface = "plain", label_size = 22,
              ncol = 2) -> fig_load_frohroh


cowplot::plot_grid(fig_roh_dist, fig_load_frohroh, 
                   labels = c("a", ""), label_fontface = "plain", label_size = 22,
                   ncol = 1) -> fig_load_inbreeding

cowplot::plot_grid(fig_countgerp, fig_countsnpef,
                   ncol = 1, labels = c("d", "e"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> fig_load_g2

cowplot::plot_grid(fig_load_g2, fig_countsnpef_cat, rel_widths = c(0.4, 0.6),
                   labels = c("", "f"), label_fontface = "plain", label_size = 22,
                   ncol = 2) -> fig_load_h2

cowplot::plot_grid(fig_load_inbreeding, fig_load_h2, 
                  ncol = 1) -> fig_load_i

cowplot::plot_grid(fig_hist_af, fig_overlap_heat,
                   fig_distribution, fig_correlation,
                   labels = c("g", "h", "i", "j"), label_fontface = "plain", label_size = 22,
                   ncol = 2,
                   align = "hv", axis = "lb") -> fig_load_j2

cowplot::plot_grid(fig_load_i,fig_load_j2, rel_heights = c(1, 0.8),
                   ncol = 1, 
                   align = "hv", axis = "lb") -> fig_load_v5

ggsave(fig_load_v5, file = "plots/Figure_1.png", width = 16, height = 24)
