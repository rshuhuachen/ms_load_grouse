#### packages ####
pacman::p_load(tidyverse, brms, bayesplot)

#### theme ####
source("scripts/theme_ggplot.R")

#### load all GO posteriors ####
files <- list.files(path = "output/biological_pathways", pattern = "lms_per_region*", full.names = T)

list <- list()
for (i in 1:length(files)){
  load(file=files[[i]])
  list[[i]] <- fit
}

names(list) <- gsub("output/biological_pathways/|.RData", "", files)

# ##### diagnose ####
# 
# diagnose_summary <- list()
# for (i in 1:length(list)){
#   #load fit
#   fit <- list[[i]]
#   #get posteriors
#   posterior <- as.array(fit)
#   log_ps <- log_posterior(fit)
#   nuts <- nuts_params(fit) #divergence
#   #get only beta and sd
#   betas <- variables(fit)[grep("b_", variables(fit))]
#   sd <- variables(fit)[grep("sd_", variables(fit))]
#   
#   #global patterns in divergence
#   diverge_beta <- mcmc_parcoord(posterior, np = nuts, pars= betas)
#   #diverge_sd <- mcmc_parcoord(posterior, np = nuts, pars= sd)
#   
#   #identify collinearity between parameters
#   collin_beta <- mcmc_pairs(posterior, np = nuts, pars= betas)
#   #collin_sd <- mcmc_pairs(posterior, np = nuts, pars= sd)
#   
#   #traceplot
#   trace_beta <- mcmc_trace(posterior, pars = betas, np = nuts)
#   trace_sd <- mcmc_trace(posterior, pars = sd, np = nuts)
#   
#   #rhat
#   rhat <- mcmc_rhat(brms::rhat(fit))
#   
#   #effective sample size
#   neff <- mcmc_neff(neff_ratio(fit))
#   
#   #autocorrelation
#   autocor_beta <- mcmc_acf(posterior, pars = betas)
#   autocor_sd <- mcmc_acf(posterior, pars=sd)
#   
#   #quick glance results
#   areas <- mcmc_areas(fit, pars=betas)
#   
#   #combine in list
#   diagnosis <- list(diverge_beta = diverge_beta, 
#                     # diverge_sd = diverge_sd, 
#                     collin_beta = collin_beta, 
#                     # collin_sd = collin_sd, 
#                     trace_beta = trace_beta, 
#                     trace_sd = trace_sd, 
#                     rhat = rhat, 
#                     neff = neff, 
#                     autocor_beta = autocor_beta, 
#                     autocor_sd = autocor_sd,
#                     areas = areas)
#   
#   
#   modelname <- sub(".*/", "", names(list[i]))
#   modelname <- sub("brm_out_", "", modelname)
#   
#   # save output
#   save(diagnosis, file=paste0("output/biological_pathways/diagnosis/", modelname, ".RData"))
#   
#   # add to summary
#   diagnose_summary[[modelname]] <- diagnosis
# }

##### prepare plot ####
source("scripts/10_biological_processes/function_posteriors.R")

# add a column that specifies the method (gerp/snpeff) and the region (exon/promoter/intron)
total_out <- list()
for(i in 1:length(list)){
  out <- get_posterior_load(posteriors = list[[i]], name = gsub("brm_out_", "", names(list[i])))
  if(grepl("gerp", names(list[i]))){method = "GERP"}else
    if(grepl("high", names(list[i]))){method = "SnpEff"}
  if(grepl("exon", names(list[i]))){region = "Exon"}else 
    if(grepl("promo", names(list[i]))){region = "Promoter"}else
      if(grepl("intron", names(list[i]))){region = "Intron"}
  
  out$method = method
  out$region = region
  total_out[[i]] <- out
}

names(total_out) <- names(list)

## add methods and region to area and interval
for (i in 1:length(total_out)){
  total_out[[i]]$interval$method <- total_out[[i]]$method[1] 
  total_out[[i]]$interval$region <- total_out[[i]]$region[1]
  total_out[[i]]$area$method <- total_out[[i]]$method[1] 
  total_out[[i]]$area$region <- total_out[[i]]$region[1] 
}

### GERP ####
### to make a big summary plot of those relationships we want to test for, bind all areas and all intervals 
all_areas_gerp <- rbind(total_out$lms_per_region_exon_gerp_androgen$area,
                   total_out$lms_per_region_exon_gerp_cellularresp$area,
                   total_out$lms_per_region_exon_gerp_developmentgrowth$area,
                   total_out$lms_per_region_exon_gerp_immune$area,
                   total_out$lms_per_region_exon_gerp_muscle$area,
                   total_out$lms_per_region_exon_gerp_oxidativestress$area,
                   total_out$lms_per_region_intron_gerp_androgen$area,
                   total_out$lms_per_region_intron_gerp_cellularresp$area,
                   total_out$lms_per_region_intron_gerp_developmentgrowth$area,
                   total_out$lms_per_region_intron_gerp_immune$area,
                   total_out$lms_per_region_intron_gerp_muscle$area,
                   total_out$lms_per_region_intron_gerp_oxidativestress$area,
                   total_out$lms_per_region_promo_gerp_androgen$area,
                   total_out$lms_per_region_promo_gerp_cellularresp$area,
                   total_out$lms_per_region_promo_gerp_developmentgrowth$area,
                   total_out$lms_per_region_promo_gerp_immune$area,
                   total_out$lms_per_region_promo_gerp_muscle$area,
                   total_out$lms_per_region_promo_gerp_oxidativestress$area)

all_intervals_gerp <- rbind(total_out$lms_per_region_exon_gerp_androgen$interval,
                            total_out$lms_per_region_exon_gerp_cellularresp$interval,
                            total_out$lms_per_region_exon_gerp_developmentgrowth$interval,
                            total_out$lms_per_region_exon_gerp_immune$interval,
                            total_out$lms_per_region_exon_gerp_muscle$interval,
                            total_out$lms_per_region_exon_gerp_oxidativestress$interval,
                            total_out$lms_per_region_intron_gerp_androgen$interval,
                            total_out$lms_per_region_intron_gerp_cellularresp$interval,
                            total_out$lms_per_region_intron_gerp_developmentgrowth$interval,
                            total_out$lms_per_region_intron_gerp_immune$interval,
                            total_out$lms_per_region_intron_gerp_muscle$interval,
                            total_out$lms_per_region_intron_gerp_oxidativestress$interval,
                            total_out$lms_per_region_promo_gerp_androgen$interval,
                            total_out$lms_per_region_promo_gerp_cellularresp$interval,
                            total_out$lms_per_region_promo_gerp_developmentgrowth$interval,
                            total_out$lms_per_region_promo_gerp_immune$interval,
                            total_out$lms_per_region_promo_gerp_muscle$interval,
                            total_out$lms_per_region_promo_gerp_oxidativestress$interval)

# relevel
all_areas_gerp$region <- factor(all_areas_gerp$region, levels = c("Promoter", "Intron", "Exon"))
all_intervals_gerp$region <- factor(all_intervals_gerp$region, levels = c("Promoter", "Intron", "Exon"))

# conditional colour of the plot

all_intervals_gerp <- all_intervals_gerp %>% mutate(Significance = as.factor(case_when(
  ll < 0 & hh < 0 ~ "95% CI doesn't overlap 0", 
  ll > 0 & hh > 0 ~ "95% CI doesn't overlap 0",
  TRUE ~ "95% CI overlaps 0"
)))

all_intervals_gerp$Significance <- factor(all_intervals_gerp$Significance, 
                                     levels=c("95% CI doesn't overlap 0", 
                                              "95% CI overlaps 0"))

all_intervals_gerp$name <- paste0(all_intervals_gerp$trait, "_", all_intervals_gerp$go, "_", all_intervals_gerp$region)
all_areas_gerp$name <- paste0(all_areas_gerp$trait, "_", all_areas_gerp$go, "_", all_areas_gerp$region)

# rename the go terms
all_intervals_gerp$go <- gsub("androgen", "Androgen metabolism", all_intervals_gerp$go)
all_intervals_gerp$go <- gsub("cellularresp", "Cellular respiration", all_intervals_gerp$go)
all_intervals_gerp$go <- gsub("developmentgrowth", "Developmental growth", all_intervals_gerp$go)
all_intervals_gerp$go <- gsub("immune", "Immune response", all_intervals_gerp$go)
all_intervals_gerp$go <- gsub("muscle", "Muscle tissue development", all_intervals_gerp$go)
all_intervals_gerp$go <- gsub("oxidativestress", "Response to oxidative stress", all_intervals_gerp$go)


all_areas_gerp$go <- gsub("androgen", "Androgen metabolism", all_areas_gerp$go)
all_areas_gerp$go <- gsub("cellularresp", "Cellular respiration", all_areas_gerp$go)
all_areas_gerp$go <- gsub("developmentgrowth", "Developmental growth", all_areas_gerp$go)
all_areas_gerp$go <- gsub("immune", "Immune response", all_areas_gerp$go)
all_areas_gerp$go <- gsub("muscle", "Muscle tissue development", all_areas_gerp$go)
all_areas_gerp$go <- gsub("oxidativestress", "Response to oxidative stress", all_areas_gerp$go)

# relevel
all_areas_gerp$go <- factor(all_areas_gerp$go, levels = c("Response to oxidative stress","Muscle tissue development","Immune response",
                                                          "Developmental growth", "Cellular respiration",
                                                          "Androgen metabolism"))

all_intervals_gerp$go <- factor(all_intervals_gerp$go, levels = c("Response to oxidative stress","Muscle tissue development","Immune response",
                                                                  "Developmental growth", "Cellular respiration",
                                                                  "Androgen metabolism"))

# add Significance to area
all_areas_gerp <- left_join(all_areas_gerp, all_intervals_gerp[,c("name", "Significance")], by = "name")

### plotting
# split by interval
brms <- split(all_areas_gerp, all_areas_gerp$interval)

brms$bottom <- brms$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

brms$outer$y <- 2
all_intervals_gerp$y <- 2

##### plot ####
ggplot(data = brms$outer) +  
  aes(x = .data$x, y = .data$go) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = Significance, col = Significance))+
  geom_segment(data=all_intervals_gerp, aes(x = l, xend = h, yend = go), col = "black", linewidth=3)+
  geom_segment(data=all_intervals_gerp, aes(x = ll, xend = hh, yend = go), col = "black")+
  geom_point(data=all_intervals_gerp, aes(x = m, y = go), fill="white",  col = "black", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = Standardised~beta, y = "GO term")+
  scale_fill_manual(values =alpha(c(clr_gerp, clr_grey), 0.7)) +
  scale_color_manual(values =c(clr_gerp, clr_grey)) +
  facet_grid(~region, scales = "free", labeller = label_wrap_gen(15))+
  scale_y_discrete(labels=scales::label_wrap(15))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(margin=margin(0,0,15,0)),
        strip.background = element_blank(),
        panel.spacing=unit(1.5,"lines"),
        legend.position = "none",
        strip.text.y = element_text(angle=0),
        panel.background = element_rect(fill = "white", color="white")) -> plot_all_gerp

plot_all_gerp

png(file = "plots/extended/extended_5_go_gerp.png", width=1200, height=800)
plot_all_gerp
dev.off()

#### SnpEff ####

### to make a big summary plot of those relationships we want to test for, bind all areas and all intervals to do a facetwrap
all_areas_high <- rbind(total_out$lms_per_region_exon_high_cellularresp$area,
                        total_out$lms_per_region_exon_high_developmentgrowth$area,
                        total_out$lms_per_region_exon_high_immune$area,
                        total_out$lms_per_region_exon_high_muscle$area,
                        total_out$lms_per_region_exon_high_oxidativestress$area,
                        total_out$lms_per_region_intron_high_cellularresp$area,
                        total_out$lms_per_region_intron_high_developmentgrowth$area,
                        total_out$lms_per_region_intron_high_immune$area,
                        total_out$lms_per_region_intron_high_muscle$area,
                        total_out$lms_per_region_intron_high_oxidativestress$area,
                        total_out$lms_per_region_promo_high_cellularresp$area,
                        total_out$lms_per_region_promo_high_developmentgrowth$area,
                        total_out$lms_per_region_promo_high_immune$area,
                        total_out$lms_per_region_promo_high_muscle$area,
                        total_out$lms_per_region_promo_high_oxidativestress$area)

all_intervals_high <- rbind(total_out$lms_per_region_exon_high_cellularresp$interval,
                            total_out$lms_per_region_exon_high_developmentgrowth$interval,
                            total_out$lms_per_region_exon_high_immune$interval,
                            total_out$lms_per_region_exon_high_muscle$interval,
                            total_out$lms_per_region_exon_high_oxidativestress$interval,
                            total_out$lms_per_region_intron_high_cellularresp$interval,
                            total_out$lms_per_region_intron_high_developmentgrowth$interval,
                            total_out$lms_per_region_intron_high_immune$interval,
                            total_out$lms_per_region_intron_high_muscle$interval,
                            total_out$lms_per_region_intron_high_oxidativestress$interval,
                            total_out$lms_per_region_promo_high_cellularresp$interval,
                            total_out$lms_per_region_promo_high_developmentgrowth$interval,
                            total_out$lms_per_region_promo_high_immune$interval,
                            total_out$lms_per_region_promo_high_muscle$interval,
                            total_out$lms_per_region_promo_high_oxidativestress$interval)

# relevel
all_areas_high$region <- factor(all_areas_high$region, levels = c("Promoter", "Intron", "Exon"))
all_intervals_high$region <- factor(all_intervals_high$region, levels = c("Promoter", "Intron", "Exon"))

# conditional colour of the plot

all_intervals_high <- all_intervals_high %>% mutate(Significance = as.factor(case_when(
  ll < 0 & hh < 0 ~ "95% CI doesn't overlap 0", 
  ll > 0 & hh > 0 ~ "95% CI doesn't overlap 0",
  TRUE ~ "95% CI overlaps 0"
)))

all_intervals_high$Significance <- factor(all_intervals_high$Significance, 
                                          levels=c("95% CI doesn't overlap 0", 
                                                   "95% CI overlaps 0"))

# go terms
all_intervals_high$go <- gsub("cellularresp", "Cellular respiration", all_intervals_high$go)
all_intervals_high$go <- gsub("developmentgrowth", "Developmental growth", all_intervals_high$go)
all_intervals_high$go <- gsub("immune", "Immune response", all_intervals_high$go)
all_intervals_high$go <- gsub("muscle", "Muscle tissue development", all_intervals_high$go)
all_intervals_high$go <- gsub("oxidativestress", "Response to oxidative stress", all_intervals_high$go)

all_areas_high$go <- gsub("cellularresp", "Cellular respiration", all_areas_high$go)
all_areas_high$go <- gsub("developmentgrowth", "Developmental growth", all_areas_high$go)
all_areas_high$go <- gsub("muscle", "Muscle tissue development", all_areas_high$go)
all_areas_high$go <- gsub("oxidativestress", "Response to oxidative stress", all_areas_high$go)
all_areas_high$go <- gsub("immune", "Immune response", all_areas_high$go)

# relevel
all_areas_high$go <- factor(all_areas_high$go, levels = c("Response to oxidative stress","Muscle tissue development","Immune response",
                                                          "Developmental growth", "Cellular respiration"))

all_intervals_high$go <- factor(all_intervals_high$go, levels = c("Response to oxidative stress","Muscle tissue development","Immune response",
                                                                  "Developmental growth", "Cellular respiration"))

#add Significance to area
all_areas_high <- left_join(all_areas_high, all_intervals_high[,c("name", "Significance")], by = "name")

### plotting
# split by interval
brms_high <- split(all_areas_high, all_areas_high$interval)

brms_high$bottom <- brms_high$outer %>%
  summarise(
    ll = min(.data$x),
    hh = max(.data$x),
    .groups = "drop_last"
  ) %>%
  ungroup()

brms_high$outer$y <- 2
all_intervals_high$y <- 2

##### plot ####
ggplot(data = brms_high$outer) +  
  aes(x = .data$x, y = .data$go) + 
  geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = Significance, col = Significance))+
  geom_segment(data=all_intervals_high, aes(x = l, xend = h, yend = go), col = "black", linewidth=3)+
  geom_segment(data=all_intervals_high, aes(x = ll, xend = hh, yend = go), col = "black")+
  geom_point(data=all_intervals_high, aes(x = m, y = go), fill="white",  col = "black", shape=21, size = 4) + 
  geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
  labs(x = Standardised~beta, y = "Density")+
  scale_fill_manual(values =alpha(c(clr_high,  clr_grey, clr_grey), 0.7)) +
  scale_color_manual(values =c(clr_high,  clr_grey, clr_grey)) +
  facet_grid(~region, scales = "free", labeller = label_wrap_gen(15))+
  scale_y_discrete(labels=scales::label_wrap(15))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing=unit(1.5,"lines"),
        strip.background = element_blank(),
        plot.title=element_text(margin=margin(0,0,15,0)),
        legend.position = "none",
        strip.text.y = element_text(angle=0),
        panel.background = element_rect(fill = "white", color="white")) -> plot_all_high
plot_all_high

png(file = "plots/extended/extended_6_go_snpeff.png", width=1200, height=800)
plot_all_high
dev.off()

## export intervals ####
interval_simple_gerp <- data.frame(parameter = all_intervals_gerp$parameter,
                                   region = all_intervals_gerp$region,
                                   go = all_intervals_gerp$go,
                                 median = round(all_intervals_gerp$m, 2),
                                 ci_95 = paste0(round(all_intervals_gerp$ll, 2), ", ", round(all_intervals_gerp$hh, 2)),
                                 ci_80 = paste0(round(all_intervals_gerp$l, 2), ", ", round(all_intervals_gerp$h, 2)))

write_tsv(interval_simple_gerp, file = "output/biological_pathways/intervals_gerp.tsv")

interval_simple_high <- data.frame(parameter = all_intervals_high$parameter,
                                   region = all_intervals_high$region,
                                   go = all_intervals_high$go,
                                   median = round(all_intervals_high$m, 2),
                                   ci_95 = paste0(round(all_intervals_high$ll, 2), ", ", round(all_intervals_high$hh, 2)),
                                   ci_80 = paste0(round(all_intervals_high$l, 2), ", ", round(all_intervals_high$h, 2)))

write_tsv(interval_simple_high, file = "output/biological_pathways/intervals_high.tsv")


