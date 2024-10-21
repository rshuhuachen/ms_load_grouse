
##### subsampling #####

#### gerp vs snpeff ####
### all gerp
load(file="output/random_draws/all_gerp.RData")

### all snpeff 
load(file="output/random_draws/all_high.RData")


## load data
load(file="results/gerp/random_draws_exon.RData")
subset_exon <- out
load(file="results/gerp/random_draws_intron.RData")
subset_intron <- out
load(file="results/gerp/random_draws_tss.RData")
subset_tss <- out
load(file="results/gerp/random_draws_promoter.RData")
subset_promoter <- out

all_subset <- rbind(subset_exon$all_draws, subset_intron$all_draws, subset_tss$all_draws, subset_promoter$all_draws)
all_subset$region <- factor(all_subset$region, levels = c("Promoter", "TSS", "Intron", "Exon"))

sum <- all_subset %>%
  group_by(region) %>%
  summarize(lower_95 = quantile(beta, probs=c(0.025)),
            upper_95 = quantile(beta, probs=c(0.975)),
            lower_80 = quantile(beta, probs=c(0.1)),
            upper_80 = quantile(beta, probs=c(0.9)),
            mean = mean(beta))

ggplot(all_subset, aes(x = beta)) + 
  xlim(-1.2,1.2)+
  facet_wrap(~region, ncol=2, strip.position = "left")+
  geom_histogram(fill = alpha(clr_gerp, 0.5), col = clr_gerp, linewidth=0.5)+
  geom_segment(data=sum, aes(x = lower_95, 
                             xend = upper_95, 
                             y = 0), col = "black", linewidth=1)+
  geom_segment(data=sum, aes(x = lower_80, 
                             xend = upper_80, 
                             y = 0), col = "black", linewidth=3)+
  geom_point(data=sum,aes(x = mean, y = 0), fill="white",  col = "black", shape=21, size = 6)+
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash", linewidth=0.6)+
  theme(plot.title = element_text(margin=margin(0,0,30,0)))+
  ylim(-100, 800)+  labs(x = expression("Standardised"~beta), y = "Number of random draws") -> draws
