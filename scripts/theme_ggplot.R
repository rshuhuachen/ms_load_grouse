
theme_set(theme_classic() + theme(title = element_text(size=18),
                                  plot.subtitle = element_text(size=14),
                                  axis.title = element_text(size = 20, family = "Arial"),
                                  axis.text = element_text(size = 18, family = "Arial"),
                                  text=element_text(size=18, family = "Arial"),
                                  legend.text =  element_text(size = 18, family = "Arial"),
                                  legend.title = element_text(size = 18, family = "Arial"),
                                  strip.text = element_text(size = 18, family = "Arial"),
                                  axis.title.y = element_text(margin = margin(t = 0, r =10, b = 0, l = 0),
                                                              color = "black"),
                                  plot.margin = margin(0.75,0.75,0.75,0.75, "cm"),
                                  plot.title=element_text(margin=margin(0,0,15,0)),
                                  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                                              color = "black"),
                                  panel.background = element_rect(fill = "white", colour = NA),
                                  plot.background = element_rect(fill = "white", colour = NA),))

#pacman::p_load(viridis, prismatic, ggsci)
#clrs <- viridis::inferno(n = 16) %>% color()
#clrs2 <- viridis::turbo(n = 16) %>% color()
clrs <- c("#284651", "#d34e38", "#54b9c5", "#efefef", "#8e8e99", "#F2CC8F", "#861657", "#CFEE9E",
          "#f49c84", "#ec7c64", "#f4b4a4", "#bccccc", "#3c5587", "#FFBC42", "#8DAA91",
          "#7E6148", "#06a088")  #from ggsci:npg + additions

clrs_hunting <- c("#001524", "#15616D","#78290F", "#FF7D00", "#FFDDAD" )
clrs_related <- c("#3c5587","#f49c84","#FFBC42","#06a088","#861657","#7E6148")
clr_gerp <- "#d34e38" 
clr_high <- "#284651"
clr_froh <- "#8DAA91" 
clr_grey <- "#8e8e99"
clr_highlight <- "#FFBC42" 
clr_2 <- c("#bccccc", "#8e8e99")
