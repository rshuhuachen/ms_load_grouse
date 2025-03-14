get_posterior_load <- function(posteriors, name){
  #name = names(posteriors)
  out = list() # to save all important output
  pacman::p_load(brms, bayesplot, ggridges)
  source("scripts/theme_ggplot.R")
  
  #extract the trait and go as strings from the name
  if(grepl("spots", name)){ #spots contain another underscore (spots_tip)
    trait = "spots"
    go = sub(".*_", "", name)
  }else{
    trait = sub("_.*", "", name)
    go = sub(".*_", "", name)
  }
  
  # output complete interval
  interval <- as.data.frame(mcmc_intervals_data(posteriors, prob =0.8, prob_outer = 0.95))
  interval$trait <- trait
  interval$go <- go
  write.csv(interval, 
            file = paste0("output/biological_pathways/intervals/interval_", name, ".csv"), quote=F, row.names=F)
  
  # get interval for just load
  interval_load <- interval %>% 
    subset(parameter == "b_scaletotal_load")
  
  # save interval
  out[["interval"]] = interval_load
  
  # get areas
  
  area <- mcmc_areas_data(posteriors) %>%
    subset(parameter == "b_scaletotal_load")
  area$trait <- trait
  area$go <- go
  
  # save area
  out[["area"]] = area
  
  ### plotting
  # split by interval
  brms <- split(area, area$interval)
  
  brms$bottom <- brms$outer %>%
    summarise(
      ll = min(.data$x),
      hh = max(.data$x),
      .groups = "drop_last"
    ) %>%
    ungroup()
  
  out[["split"]] = brms
  
  # conditional colour of the plot
  if(interval_load$l < 0 & interval_load$h < 0 | interval_load$l > 0 & interval_load$h > 0){
    clr = clr_gerp
  }else if(interval_load$ll < 0 & interval_load$hh < 0 | interval_load$ll > 0 & interval_load$hh > 0){
    clr = clr_high
  }else{clr = clr_grey}
    
  # plot
  plot <- ggplot(data = brms$outer) +  
    aes(x = .data$x, y = .data$trait) + 
    geom_ridgeline(aes(scale = 0.4, height = scaled_density, fill = trait, col = trait))+
    geom_segment(data=interval_load, aes(x = l, xend = h, yend = trait), col = "black", linewidth=3)+
    geom_segment(data=interval_load, aes(x = ll, xend = hh, yend = trait), col = "black")+
    geom_point(data=interval_load, aes(x = m, y = trait), fill="white",  col = "black", shape=21, size = 6) + 
    geom_vline(xintercept = 0, col = "#ca562c", linetype="longdash")+
    labs(x = "Beta coefficient total load", y = "Density", title = paste0("Trait = ", trait, "; GO term = ", go))+
    scale_fill_manual(values =alpha(c(clr, 0.7))) +
    scale_color_manual(values =c(clr)) +
    theme(axis.text.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = "none") 
  
  png(file = paste0("plots/biological_pathways/posteriors_", name, ".png"), width=600, height=500)
  print(plot)
  dev.off()
  
  out[["plot"]] = plot
  
  return(out)
  
}


