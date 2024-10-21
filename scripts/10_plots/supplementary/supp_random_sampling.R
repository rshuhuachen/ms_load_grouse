### load packages ####
pacman::p_load(tidyverse, ggpubr, extrafont, cowplot, data.table)

source("scripts/theme_ggplot.R")

#### total GERP #####
load(file="output/random_draws/all_gerp.RData")
load(file="output/random_draws/all_high.RData")

sum_gerp <- draws_gerp %>%
  summarize(lower_95 = quantile(beta, probs=c(0.025)),
            upper_95 = quantile(beta, probs=c(0.975)),
            lower_80 = quantile(beta, probs=c(0.1)),
            upper_80 = quantile(beta, probs=c(0.9)),
            mean = mean(beta))

ggplot(draws_gerp, aes(x = beta)) + 
  xlim(-1.2,1.2)+
  ylim(-50, 600)+
  geom_histogram(aes(fill = beta < 0, col = beta < 0), linewidth=0.5, bins=40)+
  scale_fill_manual(values =alpha(c("grey60", clr_gerp), 0.7)) + #
  scale_color_manual(values =c("grey60", clr_gerp)) +
  geom_segment(data=sum_gerp, aes(x = lower_95, 
                             xend = upper_95, 
                             y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum_gerp, aes(x = lower_80, 
                             xend = upper_80, 
                             y = 0), col = "black", linewidth=3)+
  geom_point(data=sum_gerp,aes(x = mean, y = 0), fill="white",  col = "black", shape=21, size = 6)+
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(3,"lines"),
        strip.background = element_blank(),
        legend.position="none")+
  labs(x = expression("Standardised"~beta), y = "# draws", title = "GERP") -> total_gerp

total_gerp

#### total snpeff #####
load(file="output/random_draws/all_gerp.RData")

totals <- rbind(draws_gerp, draws_high)
totals$method <- factor(totals$method, levels = c("GERP", "High impact"))

sum_high <- draws_high %>%
  summarize(lower_95 = quantile(beta, probs=c(0.025)),
            upper_95 = quantile(beta, probs=c(0.975)),
            lower_80 = quantile(beta, probs=c(0.1)),
            upper_80 = quantile(beta, probs=c(0.9)),
            mean = mean(beta))

ggplot(draws_high, aes(x = beta)) + 
  xlim(-1.2,1.2)+
  ylim(-50, 600)+
  geom_histogram(aes(fill = beta < 0, col = beta < 0), linewidth=0.5, bins=40)+
  scale_fill_manual(values =alpha(c("grey60", clr_high), 0.7)) + #
  scale_color_manual(values =c("grey60", clr_high)) +
  geom_segment(data=sum_high, aes(x = lower_95, 
                                    xend = upper_95, 
                                    y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum_high, aes(x = lower_80, 
                                    xend = upper_80, 
                                    y = 0), col = "black", linewidth=3)+
  geom_point(data=sum_high,aes(x = mean, y = 0), fill="white",  col = "black", shape=21, size = 6)+
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(3,"lines"),
        strip.background = element_blank(),
        legend.position = "none")+
  labs(x = expression("Standardised"~beta), y = "# draws", title = "SnpEff") -> total_snpeff

total_snpeff

#### regions gerp #####
## load data
load(file="output/random_draws/exon_gerp.RData")
load(file="output/random_draws/promoter_gerp.RData")
load(file="output/random_draws/tss_gerp.RData")
load(file="output/random_draws/intron_gerp.RData")

draws_gerp_promoter$region <- "Promoter"
draws_gerp_tss$region <- "TSS"
draws_gerp_intron$region <- "Intron"
draws_gerp_exon$region <- "Exon"

gerp_regions <- rbind(draws_gerp_promoter, draws_gerp_tss, draws_gerp_intron, draws_gerp_exon)
gerp_regions$region <- factor(gerp_regions$region, levels = c("Promoter", "TSS", "Intron", "Exon"))

sum_gerp <- gerp_regions %>%
  group_by(region) %>%
  summarize(lower_95 = quantile(beta, probs=c(0.025)),
            upper_95 = quantile(beta, probs=c(0.975)),
            lower_80 = quantile(beta, probs=c(0.1)),
            upper_80 = quantile(beta, probs=c(0.9)),
            mean = mean(beta))

ggplot(gerp_regions, aes(x = beta)) + 
 facet_wrap(~region, ncol=1, strip.position = "top", scales="free_y")+
  geom_histogram(aes(fill = beta < 0, col = beta < 0), linewidth=0.5, bins=40)+
  scale_fill_manual(values =alpha(c("grey60", clr_gerp), 0.7)) + #
  scale_color_manual(values =c("grey60", clr_gerp)) +geom_segment(data=sum_gerp, aes(x = lower_95, 
                             xend = upper_95, 
                             y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum_gerp, aes(x = lower_80, 
                             xend = upper_80, 
                             y = 0), col = "black", linewidth=3)+
  geom_point(data=sum_gerp,aes(x = mean, y = 0), fill="white",  col = "black", shape=21, size = 6)+
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(3,"lines"),
        strip.text = element_text(margin = margin(0,0,10,0)),
        strip.background = element_blank(),
        legend.position = "none")+
  ylim(-150, 600)+
  labs(x = expression("Standardised"~beta), y = "# draws", title = "GERP") -> regions_gerp
regions_gerp

#### regions high #####
## load data
load(file="output/random_draws/exon_high.RData")
load(file="output/random_draws/promoter_high.RData")
load(file="output/random_draws/tss_high.RData")
load(file="output/random_draws/intron_high.RData")

draws_high_promoter$region <- "Promoter"
draws_high_tss$region <- "TSS"
draws_high_intron$region <- "Intron"
draws_high_exon$region <- "Exon"

high_regions <- rbind(draws_high_promoter, draws_high_tss, draws_high_intron, draws_high_exon)
high_regions$region <- factor(high_regions$region, levels = c("Promoter", "TSS", "Intron", "Exon"))

sum_high <- high_regions %>%
  group_by(region) %>%
  summarize(lower_95 = quantile(beta, probs=c(0.025)),
            upper_95 = quantile(beta, probs=c(0.975)),
            lower_80 = quantile(beta, probs=c(0.1)),
            upper_80 = quantile(beta, probs=c(0.9)),
            mean = mean(beta))

ggplot(high_regions, aes(x = beta)) + 
  facet_wrap(~region, ncol=1, strip.position = "top", scales="free_y")+
  geom_histogram(aes(fill = beta < 0, col = beta < 0), linewidth=0.5, bins=40)+
  scale_fill_manual(values =alpha(c("grey60", clr_high), 0.7)) + #
  scale_color_manual(values =c("grey60", clr_high)) +geom_segment(data=sum_high, aes(x = lower_95, 
                                                                                xend = upper_95, 
                                                                                y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum_high, aes(x = lower_95, 
                             xend = upper_95, 
                             y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum_high, aes(x = lower_80, 
                             xend = upper_80, 
                             y = 0), col = "black", linewidth=3)+
  geom_point(data=sum_high,aes(x = mean, y = 0), fill="white",  col = "black", shape=21, size = 6)+
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(margin = margin(0,0,10,0)),
        panel.spacing = unit(3,"lines"),
        strip.background = element_blank(),
        legend.position = "none")+
  ylim(-150, 600)+
  labs(x = expression("Standardised"~beta), y = "# draws", title = "SnpEff") -> regions_high
regions_high

#### Combine ####
cowplot::plot_grid(total_gerp,  total_snpeff, regions_gerp,  regions_high, 
                   ncol = 2, align = "hv", axis = "lb",rel_heights = c(0.3, 1),
                   labels = "auto", label_fontface = "plain", label_size = 22) -> fig

png("plots/sup/fig_random_draws.png", height = 1000, width = 1000)
fig
dev.off()
