extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra, patchwork)
source("scripts/theme_ggplot.R")

#load all gerp scores
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

## violin plot
high$impact <- "High"
mod$impact <- "Moderate"
low$impact <- "Low"
modify$impact <- "Modify"

snpeff <- rbind(high, mod, low, modify)
snpeff$impact <- factor(snpeff$impact, levels = c("Modify", "Low", "Moderate", "High"))

sum_gerp_per_snpeff <- data.frame(impact = c("High", "Moderate", "Low", "Modify"),
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


test <- snpeff %>% sample_n(1000)
ggplot(snpeff, aes(x = rs_score, y = impact)) + 
  geom_vline(xintercept = 0, col = "darkred", linetype = "dotted") + 
  geom_violin(fill = "grey") + 
  #  geom_segment(data = sum_gerp_per_snpeff, aes(x = quant_25, xend = quant_75, yend = impact), 
  #               col = "black", linewidth=0.3)+
  geom_point(data = sum_gerp_per_snpeff, aes(x = median, y = impact), size = 3)+
  labs(x = "GERP score", y = "SnpEff impact category") -> violin_gerp

violin_gerp

ggsave(violin_gerp, file = "plots/sup/violin_gerp_per_snpeff.png", width=10, height=12)
