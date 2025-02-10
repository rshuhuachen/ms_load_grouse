### load packages ####
pacman::p_load(tidyverse, data.table, ggforce)
source("scripts/theme_ggplot.R")

### load average gerp scores ####
beds_windows_files <- list.files(path = "output/gerp/windows", pattern = "average_gerp_50kb*", full.names = T)

beds_windows <- data.frame()
for (i in 1:length(beds_windows_files)){
  scaf <- fread(beds_windows_files[i])
  beds_windows <- rbind(beds_windows, scaf)
}

names(beds_windows) <- c("scaf", "start_window", "end_window", "mean_gerp")

### load scaf numbers ###
load("data/genomic/raw/metadata/scaffold_names_dovetail.RData")
beds_windows <- left_join(beds_windows, genome[,c("contig", "scaf_nr")], by = c("scaf" = "contig"))
beds_windows$mean_gerp <- as.numeric(beds_windows$mean_gerp)
beds_windows <- beds_windows %>% mutate(pos = case_when(mean_gerp >0 ~ "pos",
                                                        mean_gerp < 0 ~"neg"))
### plot gerp averages ####

ggplot(subset(beds_windows, scaf_nr == 1), aes(x = start_window/10000000, y= mean_gerp)) + 
  geom_link2(aes(colour = after_stat(ifelse(y > 0, "positve", "negative")))) +
  labs(x = "Window start in Mb", y = "Mean GERP score") +
  theme(legend.position="none")+
  scale_color_manual(values = c(clr_highlight, clr_grey))+
  geom_hline(yintercept = 0, col = "darkred", linetype = "dotted")

ggplot(beds_windows, aes(x = start_window/10000000, y= mean_gerp)) + 
  geom_link2(aes(colour = after_stat(ifelse(y > 0, "positve", "negative")))) +
  labs(x = "Window start in Mb", y = "Mean GERP score") + geom_hline(yintercept = 0, col = "darkred", linetype = "dotted") + 
  facet_grid(scaf_nr~., scales="free") +
  theme(legend.position="none")+
  scale_color_manual(values = c("grey40", clr_gerp)) -> mean_gerp_across_scaf

ggsave(mean_gerp_across_scaf, file = "plots/sup/mean_gerp_across_genome_windows.png", width=20, height=24)
