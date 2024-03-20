##### Here we will make some figures to understand the relationship between gerp and snpeff better #####


#### packages ####
pacman::p_load(tidyverse, data.table, cowplot, VennDiagram)


#### data #####
load("output/load/all_loads_combined_da_nosex_29scaf.RData") #loads no sex chr only 30 scaf

source("scripts/theme_ggplot.R")

#### Relationship loads ####

ggplot(loads, aes(x = high_Lr, y = gerp_Lr_cat5)) + geom_point() +
  labs(x = "SNPeff high impact expressed load", y = "GERP category 4-5 expressed load") +
  geom_smooth(method="lm", fill = clrs[2]) -> relation_exp_load

ggplot(loads, aes(x = high_Lp, y = gerp_Lp_cat5)) + geom_point() +
  labs(x = "SNPeff high impact masked load", y = "GERP category 4-5 masked load") +
  geom_smooth(method="lm", fill = clrs[2]) -> relation_mask_load

ggplot(loads, aes(x = high_Lt, y = gerp_Lt_cat5)) + geom_point() +
  labs(x = "SNPeff high impact total load", y = "GERP category 4-5 total load") +
  geom_smooth(method="lm", fill = clrs[2]) -> relation_total_load

ggplot(loads, aes(x = high_Lt_add, y = gerp_Lt_additive_cat5)) + geom_point() +
  labs(x = "SNPeff high impact total additiveload", y = "GERP category 4-5 total additive load") +
  geom_smooth(method="lm", fill = clrs[2]) -> relation_total_load_additive

cowplot::plot_grid(relation_exp_load, relation_mask_load, relation_total_load, relation_total_load_additive,
                   ncol = 2, labs="auto", align="hv", axis="lt") %>% ggsave(file="plots/compare_load_gerp_snpef.png", width=12, height=12)


##### Overlap SNPs ######
#### packages ####
pacman::p_load(tidyverse, data.table, cowplot)

#### data #####
load("results/all_loads_combined_da_nosex_30scaf.RData") #loads no sex chr only 30 scaf

source("scripts/theme_ggplot.R")

##### Overlap SNPs ######

#load in all gerps
gerps_vcf <- data.frame()
for (i in c(1:3,5:30)){
  scaf <- read.table(paste0("output/gerp/beds/gerp_overlapSNP_scaf_",i, ".tsv.gz"))
  gerps_vcf <- rbind(gerps_vcf, scaf)
}

### write out as vcf for allele frequencies
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

# merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(gerps_vcf) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns

##### Number of SNPs in each GERP category #####

### put in GERP categories
gerps_vcf %>% select(c(chr, start, rs_score)) -> gerps_select
gerps_select <- gerps_select %>% mutate(gerp_cat = as.factor(case_when(
  rs_score <= 0 ~ "0",
  rs_score >= 0 & rs_score < 1 ~ "0-1",
  rs_score >= 1 & rs_score < 2 ~ "1-2",
  rs_score >= 2 & rs_score < 3 ~ "2-3",
  rs_score >= 3 & rs_score < 4 ~ "3-4",
  rs_score >= 4 ~ "4-5")))

gerps_select_count <- gerps_select %>% group_by(gerp_cat) %>% count()
save(gerps_select_count, file = "output/load/gerp/gerps_count_per_cat.RData")

#### Continue with calculating overlap 

# loop manually over numbers
high_05 <- nrow(subset(gerps_vcf, rs_score >= 0 & 
                         grepl("HIGH", gerps_vcf$info)))
high_15 <- nrow(subset(gerps_vcf, rs_score >= 1 & 
                         grepl("HIGH", gerps_vcf$info)))
high_25 <- nrow(subset(gerps_vcf, rs_score >= 2 & 
                         grepl("HIGH", gerps_vcf$info)))
high_35 <- nrow(subset(gerps_vcf, rs_score >= 3 & 
                         grepl("HIGH", gerps_vcf$info)))
high_45 <- nrow(subset(gerps_vcf, rs_score >= 4 & 
                         grepl("HIGH", gerps_vcf$info)))

moderate_high_05 <- nrow(subset(gerps_vcf, rs_score >= 0 & 
                         (grepl("MODERATE", gerps_vcf$info) | grepl("HIGH", gerps_vcf$info))))
moderate_high_15 <- nrow(subset(gerps_vcf, rs_score >= 1 & 
                             (grepl("MODERATE", gerps_vcf$info) | grepl("HIGH", gerps_vcf$info))))
moderate_high_25 <- nrow(subset(gerps_vcf, rs_score >= 2 & 
                             (grepl("MODERATE", gerps_vcf$info) | grepl("HIGH", gerps_vcf$info))))
moderate_high_35 <- nrow(subset(gerps_vcf, rs_score >= 3 & 
                             (grepl("MODERATE", gerps_vcf$info) | grepl("HIGH", gerps_vcf$info))))
moderate_high_45 <- nrow(subset(gerps_vcf, rs_score >= 4 & 
                             (grepl("MODERATE", gerps_vcf$info) | grepl("HIGH", gerps_vcf$info))))

low_05 <- nrow(subset(gerps_vcf, rs_score >= 0 & 
                         grepl("LOW", gerps_vcf$info)))
low_15 <- nrow(subset(gerps_vcf, rs_score >= 1 & 
                         grepl("LOW", gerps_vcf$info)))
low_25 <- nrow(subset(gerps_vcf, rs_score >= 2 & 
                         grepl("LOW", gerps_vcf$info)))
low_35 <- nrow(subset(gerps_vcf, rs_score >= 3 & 
                         grepl("LOW", gerps_vcf$info)))
low_45 <- nrow(subset(gerps_vcf, rs_score >= 4 & 
                         grepl("LOW", gerps_vcf$info)))

gerp05 <-  nrow(subset(gerps_vcf, rs_score >= 0  ))
gerp15 <-  nrow(subset(gerps_vcf, rs_score >= 1  ))
gerp25 <-  nrow(subset(gerps_vcf, rs_score >= 2  ))
gerp35 <-  nrow(subset(gerps_vcf, rs_score >= 3  ))
gerp45 <-  nrow(subset(gerps_vcf, rs_score >= 4  ))

high <- 5289 #grep -v "#" data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf | wc -l
moderate <- 137408 #grep -v "#" data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf | wc -l

overlaps <- data.frame("snpef" = c(rep(c("High", "High + moderate"), each = 5)),
                       "gerp" = c(rep(c("0-5", "1-5", "2-5", "3-5", "4-5"), times = 2)),
                       "snpef_n" = c(rep(c(high, high+moderate), each = 5)),
                       "gerp_n" = c(rep(c(gerp05, gerp15, gerp25, gerp35, gerp45), times = 2)),
                       "overlap" = c(high_05, high_15,high_25, high_35, high_45,
                                     moderate_high_05, moderate_high_15, moderate_high_25, moderate_high_35, moderate_high_45))

overlaps$perc <- overlaps$overlap / (overlaps$snpef_n + overlaps$gerp_n - overlaps$overlap) * 100
overlaps$snpef <- as.factor(overlaps$snpef)
overlaps$snpef <- factor(overlaps$snpef, levels = c("High + moderate", "High"))
overlaps$gerp <- as.factor(overlaps$gerp)

save(overlaps, file = "output/load/overlap_mutations.RData")
