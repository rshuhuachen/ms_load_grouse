pacman::p_load(tidyverse, data.table, ape, Biostrings, ggtree, treeio, phylobase, adephylo)
tree <- read.tree("output/cactus/combined_windows.fa.treefile")
source("scripts/theme_ggplot.R")

#take out ancestors remaining

tree <- drop.tip(tree, "birdAnc1")
tree <- drop.tip(tree, "birdAnc322")
tree$tip.label <- gsub("_", " ", tree$tip.label)

branchlengths <- as.data.frame(cophenetic(tree))

sum(distRoot(tree, method="patristic")) #26.32901
sum(tree$edge.length) #5.19451
#plot 
#color: 
# alectura latahmi -> last Galliformes birdAnc340
# anas to chauna -> Anseriform birdAnc341
# crypturelli to EStrutio camelus -> paleaognathae birdAnc342
# Patagioenas tasciata to Mesitornis unicolor -> Columbimorph birdAnc1
# Calvote anna to Nyctiprogne leucopyga -> Strisores birdAnc21
# Corythaeola cristata to Tauraco erythrolophus -> Musophagiformes (order)
# Chlamydotis macqueenii to Lophotis ruficrista-> Otidiformes birdAnc22
#Geococcyx californianus to Piaya cayana -> Cuculiformes birdAnc25

tree2 <- groupClade(tree, c("birdAnc25", "birdAnc22", "birdAnc21", "birdAnc34",
                            "birdAnc342","birdAnc341",  "birdAnc340"))

ggtree(tree2)+ 
  geom_treescale(x = -1, y = 1, offset = 1, label = "Substitutions per site") +  
  geom_tiplab(aes(color=group), size = 5) + theme_tree2() + xlim (0, 0.68) +
  scale_color_manual(values=c(clrs_related, clrs_hunting[1], clrs_hunting[2]),
                     labels = c("Cuculiformes",  "Otidiformes", "Musophagiformes",
                                "Strisores", "Columbimorphae", "Paleaognathae", "Anseriform",
                                "Galliformes")) +
  labs(color = "Clade") +
  theme(text = element_text(size = 16),
        legend.position = "bottom") -> treeplot

  

ggsave(treeplot, file = "plots/Supp_3_phylogenetic.png", width = 12, height = 14)
  

