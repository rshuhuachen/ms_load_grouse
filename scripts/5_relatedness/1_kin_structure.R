#### Here we analyze relatedness structure in our population based on two earlier manuscripts:
#https://heterozygosity.files.wordpress.com/2021/09/vendrami-et-al-2021.pdf
#https://heterozygosity.files.wordpress.com/2020/08/humble-et-al-2020.pdf

# packages
pacman::p_load(data.table, tidyverse)

#### Kinship Structure ####
## Based on https://github.com/elhumble/Agaz_85K_workflow_2018/blob/master/3.3_IBD.R

#### Additional filtering ####

### Here we define the VCF variables
repo = "data/genomic/intermediate/" 
VCF= paste0(repo, "rawSNPcalls.vcf")
VCFPRUNED = paste0(repo, "ltet_snps_filtered_prune")

### Prune with plink
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out ", VCFPRUNED))

## Prepare file with plink for ngsrelate
### Additional filtering for LD, high MAF, hwe

# format for ngsRelate
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# ",
              "--extract ", VCFPRUNED, ".prune.in ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 --genome ",
              "--out ", repo, "ltet_snps_filtered_kinship_ngsrelate ",
              "--recode vcf-iid"))

## format to assign categories
system(paste0("plink --vcf ", VCF, " --double-id --allow-extra-chr --set-missing-var-ids @:# ",
              "--extract ", VCFPRUNED, ".prune.in ",
              "--geno 0.1 --maf 0.01 --hwe 0.001 --genome ",
              "--out ", repo, "ltet_snps_pruned_kinship_A ", # = causes an additive component file (0/1/2)
              "--recode A"))

### Pt 1: NgsRelate to estimate KING-robust kinship, R0 and R1 for each pair of individuals ####

#https://github.com/elhumble/Agaz_85K_workflow_2018/blob/master/3.2_ngsRelate.R
#https://github.com/elhumble/SHO_reseq_2022/blob/15633040f3a48023e575af5becbca47d289bb05c/workflows/3.0_SNP_filtering_and_analysis.txt

#run ngsrelate
#NgsRelate will estimate the allele frequencies using the individuals provided in the VCF files. 
#Allele frequencies from the INFO field can used be used instead using -A TAG. 
#The TAG usually take the form of AF or AF1 but can be set to anything. 
#By default the PL data (Phred-scaled likelihoods for genotypes) is parsed, however, 
#the called genotypes can also be used instead with the -T GT option. 
#If called genotypes are being used, the software requires an additional argument (-c 1). 
#If using -c 2, ngsRelate calls genotypes assuming hardy-weinberg.

system(paste0("/prj/blackgrouse/bin/ngsRelate/ngsRelate -h ", repo, "ltet_snps_filtered_kinship_ngsrelate.vcf -T GT ",
              "-O ", repo, "ltet_snps_filtered_kinship_ngsrelate.res -c 1"))

### Load in ngsRelate output
ngsrel <- fread(paste0(repo, "ltet_snps_filtered_kinship_ngsrelate.res"), fill = TRUE)

ggplot(ngsrel) +
  geom_point(aes(R1, R0))

ggplot(ngsrel) +
  geom_point(aes(R1, KING))

### Pt 3: Assign categories with plink additive component file ####

gen <- fread(paste0(repo, "data/genomes/ngsrelate/ltet_snps_pruned_kinship_A.genome"), header = T)

summary(gen$PI_HAT)

# assign categories
### Distribution criteria based on https://academic.oup.com/bioinformatics/article/26/22/2867/228512 Manichaikul et al 2010
gen <- gen %>%
  mutate(kinship = (Z1 / 4) + (Z2 / 2)) %>%
  mutate(criteria = case_when(kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 < 0.1 ~ "Parent-offspring",
                              kinship >= 1/2^(5/2) & kinship < 1/2^(3/2) & Z0 > 0.1 & Z0 < 0.365 ~ "Full-sibling",
                              kinship >= 1/2^(7/2) & kinship < 1/2^(5/2) & Z0 > 0.365 & Z0 < 1-(1/(2^(3/2))) ~ "Second-degree",
                              kinship >= 1/2^(9/2) & kinship < 1/2^(7/2) & Z0 > 1-(1/2^(3/2)) & Z0 < 1 -(1/2^(5/2)) ~ "Third-degree",
                              kinship < 1/2^(9/2) & Z0 > 1-(1/2^(5/2)) ~ "Unrelated",
                              TRUE ~ "Unknown"))



### Plotting ###

ggplot(gen, aes(criteria)) + geom_bar() +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "black")
gen %>% group_by(criteria) %>% summarise(n = n())

## Plot IBD
# relatedness coefficients Z0, Z1 and Z2:
# reflect the proportion of the genome where a pair of individuals share zero, one or two alleles identical by descent

gen$R1 <- ngsrel$R1 
gen$R0 <- ngsrel$R0
gen$KING <- ngsrel$KING

ggplot(gen, aes(Z0, Z1)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=1)")+
  theme_classic() 

ggplot(gen, aes(Z0, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=0)") + ylab("Pr(IBD=2)")+
  theme_classic() 

ggplot(gen, aes(Z1, Z2)) + 
  geom_point(aes(colour = factor(criteria)), size = 2) +
  xlab("Pr(IBD=1)") + ylab("Pr(IBD=2)")+
  theme_classic() 

## combine with R1 and King
# individual order is the same
gen$criteria <- as.factor(gen$criteria)

ggplot(gen, aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = criteria))

ggplot(gen, aes(R1, KING)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = criteria))

ggplot(gen, aes(R1, R0)) +
  geom_point(size = 2, alpha = 0.5,
             aes(colour = criteria))



save(gen, file = "output/pop_gen/relatedness_all.RData")
write.csv(gen, file = "output/pop_gen/relatedness_all.csv", quote=F, row.names=F)

### Combine with pop data

files <- read.csv("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")
load("data/phenotypes/phenotypes_lifetime.RData")
id_sites <-pheno_wide %>% select("id", "site")
gen_select <- gen %>% select(-c("IID1", "IID2", "RT", "EZ", "Z2"))

gen_add <- left_join(gen_select, files[,c("loc", "id")], by = c("FID1" = "loc"), suffix = c("_id1", "_id2"))%>%
  left_join(files[,c("loc", "id")], by = c("FID2" = "loc"), suffix = c("_id1", "_id2")) %>%
  left_join(id_sites, by = c("id_id1" = "id"), suffix = c("_id1", "_id2")) %>%
  left_join(id_sites, by = c("id_id2" = "id"), suffix = c("_id1", "_id2"))

gen_add <- gen_add %>% mutate(same_site = case_when(
  site_id1 == site_id2 ~ "same_site",
  TRUE ~ "different_site"))
gen_add$criteria <- as.factor(gen_add$criteria)
gen_add$same_site <- as.factor(gen_add$same_site)

write.csv(gen_add, "output/pop_gen/kinship_output.csv", quote=F, row.names=F)

summary_gen_samesite <- gen_add %>% group_by(criteria, same_site) %>% mutate (n = n()) %>% select("criteria", "same_site", "n") %>% unique() 
summary_gen <- gen %>% group_by(criteria) %>% summarise(totaln = n())
summary_gen_samesite <- left_join(summary_gen_samesite, summary_gen)
summary_gen_samesite$perc <- summary_gen_samesite$n / summary_gen_samesite$totaln * 100
summary_gen_samesite <- summary_gen_samesite %>% arrange(criteria)

write.csv(summary_gen_samesite, "output/pop_gen/summary_kinship_output.csv", quote=F, row.names=F)

#make a plot
ggplot(summary_gen_samesite, aes(x = criteria, y = perc, fill = same_site)) + geom_col() + 
  annotate("text", x = c(1:6), y = c(rep(50,6)), size = 8, family = "Times",
           label = c(unique(prettyNum(summary_gen_samesite$totaln, big.mark=",", scientific=F)))) + 
  rcartocolor::scale_fill_carto_d(palette="Geyser", labels = c("Pair originates from a different site",  "Pair originates from the same site"))+
  coord_flip()+
  labs(title = "Percentage of pairwise relationships that originate 
from the same / a different site", 
       y = "Percentage of pairwise relationships",
       subtitle = "And total number of pairs in each category")+
  theme(text = element_text(family = "Times", size = 26),
        legend.text = element_text(family = "Times", size = 20),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") 
