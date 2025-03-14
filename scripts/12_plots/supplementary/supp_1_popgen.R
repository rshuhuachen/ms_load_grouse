##### Create figure 1 popgen #####

pacman::p_load(extrafont, tidyverse, data.table, readr, cowplot, ggpubr, MoMAColors, prismatic,maps, mapdata, maptools, sf, dplyr, ggplot2, cowplot, ggmap, grid, png, rnaturalearth)

#defaults
source("scripts/theme_ggplot.R")

##### Fig 1a map #####

### Load site coordinates
sites <- read.csv("data/siteinfo.csv")
sites[3,3] <- "Nyrölä" #special symbols don't get exported correctly

## Add sample sizes
N <- read.csv("data/phenotypes/phenotypes_lifetime.csv") #samples and sites
gen <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv") #genotypes
N <- subset(N, id %in% gen$id) #only take those with genotypes
N <- N %>% group_by(site) %>% count() %>% as.data.frame() #get samples per site

sites <- left_join(sites, N, by = c("abb" = "site")) #join sites with sample sizes

sites_sf = st_as_sf(sites, coords = c("lon", "lat"), #transform into correct coordinates
                    crs = 4326, agr = "constant")
sites_sf<- st_transform(sites_sf, "EPSG:4326") #transform 

## can skip: just to get coordinates of the middle in 3857 format for GIS
# sites_sf_3857<- st_transform(sites_sf, "EPSG:3857")
# lat_3857 <- c(2756602, 2779510 ,2834944 ,2780091 ,2790460 )
# lon_3857 <- c(8941486, 8920531 ,8937285 ,8937957 ,8905536 )
# 
# center_3857 = paste(min(lat_3857)+(max(lat_3857)-min(lat_3857))/2,
#                min(lon_3857)+(max(lon_3857)-min(lon_3857))/2, sep=" ")

### get world polygon
world <- ne_countries(scale = "medium", returnclass = "sf")
world$finland <- world$admin == "Finland"

# zoom in to scandinavia
world <-  st_transform(world, "EPSG:4326")

scandinavia <- world %>% 
  ggplot() + 
  geom_sf(lwd=0, aes(fill = finland, col="black"))+
  coord_sf(xlim = c(1, 35), ylim = c(58, 70))+
  scale_fill_manual(values = c("grey", clr_highlight))+
  scale_color_manual(values = c("black"))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
  geom_rect(
    xmin = 24.5,
    ymin = 62.1,
    xmax = 25.6,
    ymax = 62.4,
    fill = NA, 
    colour = "black",
    size = 0.8
  )

### make insert to zoom into the sites
## in GIS: add plugin QuickMapServices
## zoom into the center_3857 coordinates
## Project -> New print layout
## Add map -> item properties -> put min max coordinates of X and Y in 3857 (used a converter on google)
## Added legend and arrow, export

gis <- readPNG("plots/qgis_export_focalsamples.png")
gis <- rasterGrob(gis, interpolate=TRUE)

inset <- world %>% 
  ggplot() + 
  annotation_custom(gis, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  geom_sf(data = sites_sf, aes(size = n))+
  geom_text(data = sites, aes(lon-0.05, lat-0.04, label = abb), size = 5.5, family = 'Arial')+
  coord_sf(
    xlim = c(24.5, 25.6),
    ylim = c(62.1, 62.4)) +
  scale_size(range = c(2,5))+
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position=c(0.3, -0.8),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        plot.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        text=element_text(family = "Arial")) +
  guides(size=guide_legend(title="Sample size"))


## add plots together: scandinavia plus inset

ggdraw(scandinavia) +
  draw_plot(
    {
      inset
    },
    x = 0, 
    y = 0.4,
    width = 0.65, 
    height = 0.65) -> fig_1a_map

##### Fig 1b: PCA ####

## load data ###

pca <- fread("data/genomic/intermediate/ltet_snps_filtered_pruned_hwe_maf.eigenvec")
eigenval <- scan("data/genomic/intermediate/ltet_snps_filtered_pruned_hwe_maf.eigenval")
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
pheno <- read.csv("data/phenotypes/phenotypes_lifetime.csv")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca <- right_join(ids[,c(2,9)], pca, by = c("loc" = "ind"))

## add site
pca <- left_join(pca, pheno[,c(1,2,3,7,11)])

#### plot
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

#levels(pca$sitename) <- c("Kummunsuo", "Lehtusuo", "Nyrölä", "Saarisuo", "Teerisuo")
levels(pca$site) <- c("KUM", "LEH", "NYR", "SAA", "TEE")

ggplot(pca, aes(PC1, PC2, col = site)) + geom_point(size = 3) + 
  scale_color_manual(values=alpha(clrs_hunting, 0.8))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
#  labs(col = "Site")
  theme(legend.position="bottom",
        legend.title = element_text(hjust = 0.5),
        legend.key.width = unit(0.5, "cm"),
        legend.spacing.x = unit(0.5, "cm"),
        legend.key.spacing =  unit(0.3, "cm"))+
  guides(col=guide_legend(title="Site", ncol=2, byrow=TRUE)) -> fig_1b_pca

fig_1b_pca

#### Fig 1c: relatedness ####
load(file = "output/pop_gen/relatedness_all.RData")

ggplot(gen, aes(x =R1, y=KING, col = criteria)) +
  geom_point(size = 3) + 
  scale_color_manual(values=alpha(clrs_related, 0.8)) +
  labs(colour = "Criteria") +
  theme(legend.position="bottom")+
  guides(colour = guide_legend(ncol=2, bycol=TRUE))-> fig_1c_relatedness

#### Fig 1d: pairwise origins ####

summary_gen_samesite <- read.csv("output/pop_gen/summary_kinship_output.csv")
summary_gen_samesite <- subset(summary_gen_samesite, criteria != "Unknown")

#make a plot
ggplot(summary_gen_samesite, aes(x = forcats::fct_reorder(criteria, desc(criteria)), y = perc, 
                                 fill = same_site)) + geom_col() + 
  annotate("text", x = c(1:5), y = c(rep(50,5)), size = 7, family = "Arial",
           label = c(unique(prettyNum(rev(summary_gen_samesite$totaln), 
                                      big.mark=",", scientific=F)))) + 
  scale_fill_manual(values=clr_2,
                    labels = c("Different site",  
                                             "Same site"))+
  coord_flip()+
  labs(y = "Percentage of pairwise relationships", 
       x = "Criteria", fill = "Pair origin")+
  guides(fill = guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.position = "bottom") -> fig_1d_relatedness

#### Combine in multi-faceted plot ####
plot_grid(fig_1a_map, fig_1c_relatedness,
          ncol = 1, labels = c("a", "c"), label_fontface = "plain", label_size = 22,
          align="v", axis="r") -> left
plot_grid(fig_1b_pca, fig_1d_relatedness,
          ncol = 1, labels = c("b", "d"), label_fontface = "plain", label_size = 22,
          align="v", axis="l") -> right

png(bg = "white", file = "plots/sup/sup_1_popgen.png", width = 1050, height = 900)
plot_grid(left, right, ncol = 2, align = "hv", axis = "bt")
dev.off()
