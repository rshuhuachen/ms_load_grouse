##### Principal Component Analysis #####
### In this script, we first use plink to conduct two additional filtering steps: filter for LD and HWE
### Then, we use plink to conduct a principal component analysis
### We analyse the results in R and plot them

# source of tutorial: https://speciationgenomics.github.io/pca/

#### Additional filtering ####
### Set the directory
repo = "data/genomic/intermediate/" # Directory of the pre-processed VCF file 

### Here we define the VCF variables
VCF= paste0(repo, "rawSNPcalls.vcf")
VCFPRUNED = paste0(repo, "ltet_snps_filtered_prune")
VCFOUT = paste0(repo, "ltet_snps_filtered_pruned_hwe_maf")

### Prune with plink
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ", VCFPRUNED))

#### Run PCA with plink ####

#run pca with filtering for HWE and including only SNPs out of LD
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --maf 0.01 --hwe 0.001 --set-missing-var-ids @:# --extract ", VCFPRUNED, ".prune.in --make-bed --pca --out ", VCFOUT))

#### Analyse PCA in R ####

pacman::p_load(tidyverse, extrafont)

pca <- read_table(paste0(VCFOUT, ".eigenvec"), col_names = FALSE)
eigenval <- scan(paste0(VCFOUT, ".eigenval"))
ids <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
pheno <- read.csv("data/phenotypes/phenotypes_lifetime.csv")

# sort out the pca data
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

ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+ ylab("Percentage variance explained") + theme_light()

levels(pca$site) <- c("Kummunsuo", "Lehtusuo", "Nyrölä", "Saarisuo", "Teerisuo")

ggplot(pca, aes(PC1, PC2, col = site)) + geom_point(size = 4) + 
  rcartocolor::scale_color_carto_d(palette="Geyser")+
 theme_classic() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+labs(title="Principal Component Analysis (PCA) results")+
  guides(col=guide_legend(title="Lekking site"))+
  theme(text = element_text(family = "Times", size = 26),
        legend.text = element_text(family = "Times", size = 26))