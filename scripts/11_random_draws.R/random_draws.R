### load packages
pacman::p_load(tidyverse, data.table)

#load data mutations
load(file = "output/load/snpeff/snpeff_high.RData")
load(file = "output/load/gerp/gerp_over4.RData")

# load annotation data
load(file = "output/load/snpeff/snpeff_high_annotated_region.RData")
snpef_all$CHROM <- gsub("__", ";", snpef_all$CHROM)
snpef_all$CHROM <- gsub("HRSCAF_", "HRSCAF=", snpef_all$CHROM)
snpef_all$chr_pos <- paste0(snpef_all$CHROM, "_", snpef_all$POS)

load(file = "output/load/gerp/gerp_annotated_region.RData")
gerp_all$chr <- gsub("__", ";", gerp_all$chr)
gerp_all$chr <- gsub("HRSCAF_", "HRSCAF=", gerp_all$chr)
gerp_all$chr_pos <- paste0(gerp_all$chr, "_", gerp_all$pos)

### make subsets of mutations based on region
# subset snpef based on annotation
snpeff_exons <- subset(snpef_all, region_exon == 1)
snpeff_tss <- subset(snpef_all, region_tss == 1)
snpeff_introns <- subset(snpef_all, region_intron == 1)
snpeff_promoter <- subset(snpef_all, region_promoter == 1 & is.na(region_tss))

snpeff$chr_pos <- paste0(snpeff$CHROM, "_", snpeff$POS)

snpeff_exons_gt <- subset(snpeff, chr_pos %in% snpeff_exons$chr_pos)
snpeff_tss_gt <- subset(snpeff, chr_pos %in% snpeff_tss$chr_pos)
snpeff_introns_gt <- subset(snpeff, chr_pos %in% snpeff_introns$chr_pos)
snpeff_promoter_gt <- subset(snpeff, chr_pos %in% snpeff_promoter$chr_pos)

# subset gerp based on annotation
gerp_exons <- subset(gerp_all, region_exon == 1)
gerp_tss <- subset(gerp_all, region_tss == 1)
gerp_introns <- subset(gerp_all, region_intron == 1)
gerp_promoter <- subset(gerp_all, region_promoter == 1 & is.na(region_tss))

gerp$chr_pos <- paste0(gerp$chr, "_", gerp$pos)

gerp_exons_gt <- subset(gerp, chr_pos %in% gerp_exons$chr_pos)
gerp_tss_gt <- subset(gerp, chr_pos %in% gerp_tss$chr_pos)
gerp_introns_gt <- subset(gerp, chr_pos %in% gerp_introns$chr_pos)
gerp_promoter_gt <- subset(gerp, chr_pos %in% gerp_promoter$chr_pos)

## create function to subset X number of SnpEffs and model its effect on fitness

random_draws <- function(geno, n_draws, n_mutations, file, method, emperical_beta){
  source("scripts/theme_ggplot.R")
  all_draws <- data.frame()
  
  if(method == "GERP"){
    for (i in 1:n_draws){
      draw <- geno[sample(nrow(geno), n_mutations),] #randomly draw snps
      ## load functions
      source("scripts/7_calculate_load/function_calculate_load.R")
      load <- calculate_load_gerp(draw, output_vcf = F, loadtype = "random_draw") #calculate load
      model_out <- model_load(load, i)
      model_out$method <- method
      
      all_draws <- rbind(all_draws, model_out)
    }}
  
  if(method == "High impact"){
    for (i in 1:n_draws){
      draw <- geno[sample(nrow(geno), n_mutations),] #randomly draw snps
      ## load functions
      source("scripts/7_calculate_load/function_calculate_load.R")
      load <- calculate_load_snpeff(draw, output_vcf = F, loadtype = "random_draw") #calculate load
      model_out <- model_load(load, i)
      model_out$method <- method
      
      all_draws <- rbind(all_draws, model_out)
    }}
  
  ### conclusion
  all_draws <- all_draws %>% mutate(conclusion = as.factor(case_when(
    beta < 0 & pval < 0.05 ~ "Significantly negative",
    beta > 0 & pval < 0.05 ~ "Significantly positive",
    TRUE ~ "Insignificant"
  )))
  
  return(all_draws)
  save(all_draws, file = file)
}

# run functions

# all mutations
draws_gerp <- random_draws(geno = gerp, n_draws = 5000, n_mutations = 1000, file = "output/random_draws/all_gerp.RData", method="GERP", emperical_beta = -0.21)
save(draws_gerp, file = "output/random_draws/all_gerp.RData")
draws_high <- random_draws(geno = snpeff, n_draws = 5000, n_mutations = 1000, file = "output/random_draws/all_high.RData", method = "High impact", emperical_beta = -0.07)
save(draws_high, file = "output/random_draws/all_high.RData")

# per region
# high
draws_high_promoter <- random_draws(geno = snpeff_promoter_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/promoter_high.RData", method = "High impact", emperical_beta = 0)
save(draws_high_promoter, file = "output/random_draws/promoter_high.RData")
draws_high_tss <- random_draws(geno = snpeff_tss_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/tss_high.RData", method = "High impact", emperical_beta = 0)
save(draws_high_tss, file = "output/random_draws/tss_high.RData")
draws_high_intron <- random_draws(geno = snpeff_introns_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/intron_high.RData", method = "High impact", emperical_beta = 0)
save(draws_high_intron, file = "output/random_draws/intron_high.RData")
draws_high_exon <- random_draws(geno = snpeff_exons_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/exon_high.RData", method = "High impact", emperical_beta = 0)
save(draws_high_exon, file = "output/random_draws/exon_high.RData")

# gerp
draws_gerp_promoter <- random_draws(geno = gerp_promoter_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/promoter_gerp.RData", method = "GERP", emperical_beta = 0)
save(draws_gerp_promoter, file = "output/random_draws/promoter_gerp.RData")
draws_gerp_tss <- random_draws(geno = gerp_tss_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/tss_gerp.RData", method = "GERP", emperical_beta = 0)
save(draws_gerp_tss, file = "output/random_draws/tss_gerp.RData")
draws_gerp_intron <- random_draws(geno = gerp_introns_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/intron_gerp.RData", method = "GERP", emperical_beta = 0)
save(draws_gerp_intron, file = "output/random_draws/intron_gerp.RData")
draws_gerp_exon <- random_draws(geno = gerp_exons_gt, n_draws = 5000, n_mutations = 500, file = "results/random_draws/exon_gerp.RData", method = "GERP", emperical_beta = 0)
save(draws_gerp_exon, file = "output/random_draws/exon_gerp.RData")

#### load in manuallly for prettier plots 
pacman::p_load(tidyverse, data.table)
source("scripts/theme_ggplot.R")

#load(file="results/gerp/random_draws_gerp.RData")
#draws_gerp <- out
#load(file="results/gerp/random_draws_snpeff_high.RData")
#draws_high <- out

### make manual plots to adjust ####

##### gerp ######
### plot distribution conclusion
ggplot(draws_gerp, aes(x = conclusion, fill = conclusion)) + geom_histogram(stat="count") +
  scale_fill_manual(values=c(clrs[1], clrs[2], clrs[3])) + coord_flip() + theme(legend.position="none", axis.title.y=element_blank())+
  scale_y_continuous(expand=c(0,1))+
  labs(title = "GERP ≥ 4", y = "Count") + geom_label(stat = "count", aes(label = paste0("n = ", after_stat(count))), hjust = 1, size = 6, fill = "white") -> gerp_conclusion

### plot distribution of effect sizes
ggplot(draws_gerp, aes(x = beta)) + geom_histogram(fill = clr_gerp, alpha = 0.7, col = "black") + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=1.5)+
  xlim(-1.2,1.2)+
  geom_segment(aes(x = quantile(beta, probs=c(0.025)), 
                   xend = quantile(beta, probs=c(0.975)), 
                   y = 1300), col = "black", linewidth=1)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)))+
  geom_text(label = "★", aes(x = -0.21, y = 1300), col = clr_highlight, size = 8, family = "HiraKakuPro-W3")+
  scale_y_continuous(expand=c(0,0), limits = c(0, 1500))+
  labs(title = "GERP ≥ 4", x = expression(beta~"total load on LMS"), y = "Count") -> gerp_hist

##### snpeff ######
### plot distribution conclusion
ggplot(draws_high, aes(x = conclusion, fill = conclusion)) + geom_histogram(stat="count") +
  scale_fill_manual(values=c(clrs[1], clrs[2], clrs[3])) + coord_flip() + theme(legend.position="none", axis.title.y=element_blank())+
  scale_y_continuous(expand=c(0,1))+
  labs(title = "High impact SnpEff", y = "Count") + geom_label(stat = "count", aes(label = paste0("n = ", after_stat(count))), hjust = 1, size = 6, fill = "white") -> snpeff_conclusion

### plot distribution of effect sizes
ggplot(draws_high, aes(x = beta)) + geom_histogram(fill = clr_gerp, alpha = 0.7, col = "black") + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=1.5)+
  xlim(-1.2,1.2)+
  geom_segment(aes(x = quantile(beta, probs=c(0.025)), 
                   xend = quantile(beta, probs=c(0.975)), 
                   y = 1300), col = "black", linewidth=1)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)))+
  geom_text(label = "★", aes(x = -0.07, y = 1300), col = clr_highlight, size = 8, family = "HiraKakuPro-W3")+
  scale_y_continuous(expand=c(0,0), limits = c(0, 1500))+
  labs(title = "High impact SnpEff", x = expression(beta~"total load on LMS"), y = "Count") -> snpeff_hist

### Combine ####
cowplot::plot_grid(gerp_conclusion, gerp_hist, snpeff_conclusion, snpeff_hist,
                   align = "hv", axis = "lb",labels = "auto", label_fontface = "plain", label_size = 22) -> plots

png("plots/final/plots_random_draws_conclusions_gerp_snpef.png", width=1000, height=1000)
plots
dev.off()

### regions
ggplot(draws_high_tss, aes(x = beta)) + geom_histogram(fill = clr_gerp, alpha = 0.7, col = "black") + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=1.5)+
  xlim(-1.2,1.2)+
  geom_segment(aes(x = quantile(beta, probs=c(0.025)), 
                   xend = quantile(beta, probs=c(0.975)), 
                   y = 1300), col = "black", linewidth=1)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)))+
  geom_text(label = "★", aes(x = -0.21, y = 1300), col = clr_highlight, size = 8, family = "HiraKakuPro-W3")+
  scale_y_continuous(expand=c(0,0), limits = c(0, 1500))+
  labs(title = "GERP ≥ 4", x = expression(beta~"total load on LMS"), y = "Count") 

