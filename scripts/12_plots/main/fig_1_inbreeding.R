extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, cowplot, extrafont, ggtext, VennDiagram, ggvenn, gridExtra)
source("scripts/theme_ggplot.R")

#### Figure 1a: inbreeding froh histogram #### 
# load froh mcmc
rohs <- fread("output/inbreeding/bcftools_roh_rg_clean.txt")
load("output/inbreeding/froh.RData")

ggplot(froh, aes(x = froh)) + 
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  labs(x = expression(Inbreeding~coefficient~italic(F)[ROH]), y = "Number of individuals") -> fig_hist_froh

png(file = "plots/main/fig_1a.png", width=600, height=400)
fig_hist_froh
dev.off()

#### Figure 1b: inbreeding roh number histogram #### 
num_roh_per_ind <- rohs %>% group_by(id) %>% tally() 

ggplot(num_roh_per_ind, aes(n)) +
  geom_histogram(bins = 30,  position = "identity", 
                 col = clr_froh, fill = alpha(clr_froh, 0.7)) + 
  ylab("Number of individuals") +
  xlab("Number of ROHs")  -> fig_hist_nroh

# alternative: cum ROHs

#rohs_class <- rohs %>% mutate(class = as.factor(case_when(
#  length < 100000 ~ "< 100kb",
#  length >= 100000 & length < 1000000 ~ "100kb - 1Mb",
#  length >= 1*10^6 & length < 5*10^6 ~ "1 - 5Mb", 
#  length >= 5*10^6 ~ "> 5Mb"
#)))

rohs_class <- rohs %>% mutate(class = as.factor(case_when(
  length < 1*10^6 ~ "< 1 Mb",
  length >= 1*10^6 & length < 2*10^6 ~ "1 - 2Mb",
  length >= 2*10^6 ~ "> 2Mb"
)))

#rohs_class <- rohs_class %>% group_by(id, class) %>% summarise(froh = sum(length) / 1004266063)
#rohs_class$class <- factor(rohs_class$class, levels = c("< 100kb", "100kb - 1Mb", "1 - 5Mb", "> 5Mb"))
rohs_class$class <- factor(rohs_class$class, levels = c("", "< 1 Mb", "1 - 2Mb", "> 2Mb"))

rohs_class %>% group_by(id, class) %>% summarise(froh = sum(length) / 1004266063) %>% ungroup() %>%
  group_by(id) %>% arrange(id, class) %>% mutate(cum_froh = cumsum(froh)) %>% ungroup() %>%
  ggplot(aes(x = class, y = cum_froh)) + 
  geom_line(aes(group=id), col = clr_2[2])+ 
  geom_point(col = "#7E7E8B", size = 2) + 
  labs(x = "ROH length", y = expression(Cumulative~italic(F)[ROH])) -> fig_cum_froh

fig_cum_froh

png(file = "plots/main/fig_1b.png", width=600, height=400)
fig_cum_froh
dev.off()

#### Figure 1c: distribution ROHs most inbred ####

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
         KB = length / 1000)

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
  mutate(min = case_when(scaf_nr == 2 | scaf_nr == 5 | scaf_nr == 7 | scaf_nr == 9  ~ 0,
                         TRUE ~ min)) %>%
  mutate(max = case_when(scaf_nr == 2 | scaf_nr == 5 | scaf_nr == 7 | scaf_nr == 9  ~ 0,
                         TRUE ~ max))

col <- c("#008080", "#b4c8a8")

#take FROH measures from these individuals
yax <- left_join(yax, froh[,c("id", "froh")])

df %>% 
  filter(KB > 50 & scaf_nr < 11) %>% 
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
    strip.background = element_blank(),
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

png(file = "plots/main/fig_1c.png", width=800, height=600)
fig_roh_dist
dev.off()

#### Combine in one figure ####
cowplot::plot_grid(fig_hist_froh, fig_cum_froh,
                   ncol = 2, labels = c("a", "b"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> fig_froh_1

cowplot::plot_grid(fig_froh_1, fig_roh_dist,
                   ncol = 1, labels = c("", "c"), label_fontface = "plain", label_size = 22) -> fig_froh


png(file = "plots/main/fig_1.png", width=1100, height=800)
fig_froh
dev.off()

#### Figures for nee in PDF #####

# a
fig_hist_froh_nee <- fig_hist_froh + theme(axis.text = element_text(size = 12),
                                           axis.title = element_text(size = 14))


# b
fig_cum_froh_nee <- rohs_class %>% group_by(id, class) %>% summarise(froh = sum(length) / 1004266063) %>% ungroup() %>%
  group_by(id) %>% arrange(id, class) %>% mutate(cum_froh = cumsum(froh)) %>% ungroup() %>%
  ggplot(aes(x = class, y = cum_froh)) + 
  geom_line(aes(group=id), col = clr_2[2], linewidth=0.2)+ 
  geom_point(col = "#7E7E8B", size = 1) + 
  labs(x = "ROH length", y = expression(Cumulative~italic(F)[ROH])) + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
fig_cum_froh_nee

# c

df %>% 
  filter(KB > 100 & scaf_nr < 11) %>% 
  ggplot() +
  geom_rect(data=shade, aes(xmin=min, xmax=max, ymin=0, ymax=num_ind*2 + 1), 
            alpha=0.5, fill = "#eceff4") + # "#f7f7f7" "#eceff4"
  geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
  geom_rect(aes(xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9, 
                fill = as.factor(scaf_nr),  col = as.factor(scaf_nr)), size = 0, alpha = 1) + 
  scale_y_reverse(expand = c(0, 0)) +
  facet_grid(~scaf_nr,scales = 'free_x', space = 'free_x', switch = 'x') +
  scale_fill_manual(values = rep(clr_2, 28)) + 
  scale_color_manual(values = rep(clr_2, 28)) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    # plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
    axis.line.x = element_blank(),
    legend.position="none",
    axis.title.x = element_text(margin=margin(t=10)),
    axis.title.y = element_text(margin=margin(r=1)),
    axis.text.y = element_text(colour = "white"),
    axis.line.y = element_blank()) +
  coord_cartesian(clip = 'off') +
  xlab("Scaffold") +
  ylab("Individuals") -> fig_roh_dist_nee

cowplot::plot_grid(fig_hist_froh_nee, fig_cum_froh_nee,
                   ncol = 2, labels = c("a", "b"), label_fontface = "plain", label_size = 16,
                   align = "hv", axis = "lb") -> fig_froh_1_nee

cowplot::plot_grid(fig_froh_1_nee, fig_roh_dist_nee,
                   ncol = 1, labels = c("", "c"), label_fontface = "plain", label_size = 16) -> fig_froh_nee

ggsave(plot = fig_froh_nee, filename = 'plots/main/fig_1.pdf',width = 180,height = 140,
       units = 'mm', device = cairo_pdf)


