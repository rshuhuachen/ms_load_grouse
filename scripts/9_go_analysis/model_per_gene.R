###### Here, we can try an alternative modelling approach
###### we will calculate the number of mutations per gene per individual
###### and then create a model per gene predicting their effect on the 6 sexual/behavioural traits 

# packages
pacman::p_load(tidyverse, data.table, lme4, DHARMa, performance, glmmTMB)

### Prepare annotation/genotypes ####
# load annotation
load(file = "output/load/gerp/gerp_annotated_region.RData")
gerp_all$CHROM <- gsub("__", ";", gerp_all$CHROM)
gerp_all$CHROM <- gsub("HRSCAF_", "HRSCAF=", gerp_all$CHROM)

# subset gerp based on annotation
coding_gerp <- subset(gerp_all, is.na(region_noncoding))

#load genotypes
load(file = "output/load/gerp/gerps_cat5_vcfformat.RData") #vcf file format with snps GERP > 4, about 360k

#rename columns in file
filenames <- fread("data/genomic/raw/metadata/idnames.fam")
ids <- fread("data/genomic/raw/metadata/file_list_all_bgi_clean.csv")

#merge
idnames <- left_join(filenames[,c("V1")], ids[,c("loc", "id")], by = c("V1" = "loc"))

names(gerps_cat5_vcf) <- c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format", idnames$id) #rename columns

## subset based on coding region
gerps_cat5_vcf$chr_pos <- paste0(gerps_cat5_vcf$chr, "_", gerps_cat5_vcf$pos)
coding_gerp$chr_pos <- paste0(coding_gerp$CHROM, "_", coding_gerp$POS)

gerps_cat5_vcf_coding <- subset(gerps_cat5_vcf, chr_pos %in% coding_gerp$chr_pos)

#clean env
rm(gerps_cat5_vcf); rm(coding_gerp); rm(gerp_all)
rm(filenames); rm(idnames); rm(ids)

## merge with chicken annotation to have the gene name for each SNP
load(file = "output/load/gerp/chicken_annotated_gerp.RData")

#add column for the annotation from the original gff (ANNXXX)
gerp_annotated_chicken$gene_og_ann <- gsub(".*\\b(ANN\\d+).*", "\\1", gerp_annotated_chicken$gerp_gr.INFO)
gerp_annotated_chicken <- gerp_annotated_chicken %>% mutate(gene_og_ann = case_when(
  grepl("ANN", gerp_annotated_chicken$gene_og_ann) ~ gene_og_ann,
  TRUE ~ NA
))

## merge by chr_pos
gerp_annotated_chicken$chr_pos <- paste0(gerp_annotated_chicken$gerp_gr.seqnames, "_", gerp_annotated_chicken$gerp_gr.start)
gerp_annotated_chicken$chr_pos <- gsub("__", "\\;", gerp_annotated_chicken$chr_pos) #adjust scaf names
gerp_annotated_chicken$chr_pos <- gsub("HRSCAF_", "HRSCAF=", gerp_annotated_chicken$chr_pos) #adjust scaf names
gerp_annotated_chicken <- subset(gerp_annotated_chicken, !is.na(gene)) #only choose those with an annotated gene

gerp <- inner_join(gerps_cat5_vcf_coding, unique(gerp_annotated_chicken[,c("chr_pos", "gene", "gene_og_ann")]), by = "chr_pos")
summary(as.factor(gerp$gene)) #some positions have multiple gene annotations

### only include those with at least 5 mutations in one gene
n_mutation_per_gene <- gerp %>% group_by(gene) %>% tally()
n_mutation_per_gene_min5 <- subset(n_mutation_per_gene, n >= 5)

gerp_sub <- subset(gerp, gene %in% n_mutation_per_gene_min5$gene)

#### Number of mutations per gene ####

#first: recode genotypes to the number of mutations
gt <- c(11:200) #genotypes
select_n3 <- function(x){x = substr(x,1,3)}
gerp_sub[gt] <- lapply(gerp_sub[gt], select_n3)

# transform GT columns to the number of mutations 
gerp_sub[gerp_sub == "1/1"] <- 2
gerp_sub[gerp_sub == "1/0"] <- 1
gerp_sub[gerp_sub == "0/1"] <- 1
gerp_sub[gerp_sub == "0/0"] <- 0
gerp_sub[gerp_sub == "./."] <- NA
gerp_sub[gt] <- lapply(gerp_sub[gt], as.numeric)

#sum the number of mutations per gene
gerp_n_by_gene <- gerp_sub %>% dplyr::select(c(gene, D154259:D229777)) %>% group_by(gene) %>%
  summarise(across(everything(), list(sum = ~sum(., na.rm = TRUE), count = ~sum(!is.na(.)))))

# make df long
# sum of mutations
gerp_n_by_gene_mutation <- cbind(gerp_n_by_gene$gene, gerp_n_by_gene[,grep("_sum", colnames(gerp_n_by_gene))])
gerp_n_by_gene_n <- cbind(gerp_n_by_gene$gene,gerp_n_by_gene[,grep("_count", colnames(gerp_n_by_gene))])

gerp_n_by_gene_mutation_long <- gather(gerp_n_by_gene_mutation, id, n_mutations, D154259_sum:D229777_sum, factor_key=TRUE)
gerp_n_by_gene_n_long <- gather(gerp_n_by_gene_n, id, n_snps, D154259_count:D229777_count, factor_key=TRUE)

#remove _count and _sum from id column and rename columns
gerp_n_by_gene_mutation_long$id <- gsub("_sum", "", gerp_n_by_gene_mutation_long$id)
names(gerp_n_by_gene_mutation_long)[1] <- "gene"
gerp_n_by_gene_n_long$id <- gsub("_count", "", gerp_n_by_gene_n_long$id)
names(gerp_n_by_gene_n_long)[1] <- "gene"

#merge, join to exclude sum with low mutation genes too
gerp_n_by_gene_long <- left_join(gerp_n_by_gene_mutation_long, gerp_n_by_gene_n_long, by = c("gene", "id"))
gerp_n_by_gene_long$mutations_corrected <- gerp_n_by_gene_long$n_mutations / (gerp_n_by_gene_long$n_snps*2) #times 2 because we are counting alleles not snps

gerp_n_by_gene_ls <- gerp_n_by_gene_long %>% group_split(gene)
save(gerp_n_by_gene_ls, file = "output/load/gerp/mutations_per_gene_gerp5.RData")
load(file = "output/load/gerp/mutations_per_gene_gerp5.RData")

#### Modelling ####
### make a function to run a model per gene
function_glmm <- function(df, phenotype, family, generalized) {
  pacman::p_load(lme4, performance, DHARMa)
  # fit the model and null-model, do lrt, and get p-value
  df <- as.data.frame(df)
  
  # load pheno data
  load("data/phenotypes/phenotypes_annual.RData") #annual
  
  core <- subset(pheno_long, core == "core")
  
  #merge with gene data, only core individuals
  df <- left_join(df, core, by = "id")
  
  formula_alt <- formula(paste0(phenotype, " ~ mutations_corrected + age + sqrt(age) + (1|year) + (1|site/id)"))
  
  if(generalized == "lmer"){
    model <- lmerTest::lmer(formula_alt, data=df)
    
    coef <- as.data.frame(summary(model)$coefficients)
    
    mutations_estimate <- coef["mutations_corrected", "Estimate"]
    mutations_se <- coef["mutations_corrected", "Std. Error"]
    mutations_df <- coef["mutations_corrected", "df"]
    mutations_val <- coef["mutations_corrected", "t value"]
    mutations_pval <- coef["mutations_corrected", "Pr(>|t|)"]
    
    r2_conditional <- as.vector(performance::r2(model)[1]$R2_conditional)
    r2_marginal <- as.vector(performance::r2(model)[2]$R2_marginal)
    singular <- check_singularity(model)
    uniformity <- testUniformity(model, plot=F)$p.value # test of overall uniformity, KS test
    dispersion <- testDispersion(model, plot=F)$p.value # overdispersion
    outlier <- testOutliers(model, plot=F)$p.value # outlier test
    zeroinflate <- testZeroInflation(model, plot=F)$p.value # zero inflation
    quantile_dev <- testQuantiles(model, plot=F)$p.value #quantile deviations of residuals
  }
  
  if(generalized == "glmer"){
    model <- glmmTMB::glmmTMB(formula_alt, data=df, family = "poisson", ziformula = ~1)
    coef <- as.data.frame(summary(model)$coefficients$cond)
    
    mutations_estimate <- coef["mutations_corrected", "Estimate"]
    mutations_se <- coef["mutations_corrected", "Std. Error"]
    mutations_df <- NA
    mutations_val <- coef["mutations_corrected", "z value"]
    mutations_pval <- coef["mutations_corrected", "Pr(>|z|)"]
    
    r2_conditional <- as.vector(performance::r2(model)[1]$R2_conditional)
    r2_marginal <- as.vector(performance::r2(model)[2]$R2_marginal)
    singular <- check_singularity(model)
    uniformity <- testUniformity(model, plot=F)$p.value # test of overall uniformity, KS test
    dispersion <- testDispersion(model, plot=F)$p.value # overdispersion
    outlier <- testOutliers(model, plot=F)$p.value # outlier test
    zeroinflate <- testZeroInflation(model, plot=F)$p.value # zero inflation
    quantile_dev <- testQuantiles(model, plot=F)$p.value #quantile deviations of residuals
  }
  
  return(data.frame(gene=as.character(df$gene[1]), 
                    mutations_estimate = mutations_estimate,
                    mutations_se = mutations_se,
                    mutations_df = mutations_df,
                    mutations_val = mutations_val,
                    mutations_pval = mutations_pval,
                    r2_conditional = r2_conditional,
                    r2_marginal = r2_marginal,
                    singular = singular,
                    uniformity = uniformity,
                    dispersion = dispersion,
                    outlier = outlier,
                    zeroinflate = zeroinflate,
                    quantile_dev = quantile_dev))
}

#for each trait
out_lyre <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "lyre", family = "gaussian", 
                                 generalized = "lmer", mc.cores = 12)

out_eyec <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "eyec", family = "gaussian", 
                                 generalized = "lmer", mc.cores = 12)

out_blue <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "blue", family = "gaussian", 
                               generalized = "lmer", mc.cores = 12)

out_attend <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "attend", family = "gaussian", 
                               generalized = "lmer", mc.cores = 12)

out_fight <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "fight", family = "gaussian", 
                               generalized = "lmer", mc.cores = 12)

out_dist <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "dist", family = "gaussian", 
                               generalized = "lmer", mc.cores = 12)

out_ms <- parallel::mcmapply(function_glmm, df = gerp_n_by_gene_ls, phenotype = "MS", family = "poisson", 
                                 generalized = "glmer", mc.cores = 12)

list_models <- list(out_lyre = out_lyre,
                    out_eyec = out_eyec,
                    out_blue = out_blue,
                    out_attend = out_attend,
                    out_fight = out_fight,
                    out_dist = out_dist,
                    out_ms = out_ms)


for (trait in 1:length(list_models)){
  #weird workaround but it works
  gather <- vector()
  for (i in 1:length(list_models[[trait]])){
    element <- as.character(list_models[[trait]][[i]])
    gather <- c(gather, element)
  }
  
  out_df <- as.data.frame(matrix(gather, ncol = 14, byrow = TRUE))
  names(out_df) <- c("gene", "mutations_estimate" , "mutations_se",
                          "mutations_df", "mutations_val", "mutations_pval",
                          "r2_conditional", "r2_marginal", "singular", "uniformity",
                          "dispersion", "outlier", "zeroinflate", "quantile_dev")
  out_df$gene <- as.factor(out_df$gene)
  out_df[2:14] <- lapply(out_df[2:14], as.character)
  out_df[2:14] <- lapply(out_df[2:14], as.numeric)
  
  save(out_df, file = paste0("output/models/per_gene/model_per_gene_", names(list_models[trait]), ".RData"))
}

## load in results ##
load("output/models/per_gene/model_per_gene_out_lyre.RData")
out_lyre <- out_df
load("output/models/per_gene/model_per_gene_out_blue.RData")
out_blue <- out_df
load("output/models/per_gene/model_per_gene_out_eyec.RData")
out_eyec <- out_df
load("output/models/per_gene/model_per_gene_out_dist.RData")
out_dist <- out_df
load("output/models/per_gene/model_per_gene_out_fight.RData")
out_fight <- out_df
load("output/models/per_gene/model_per_gene_out_attend.RData")
out_attend <- out_df
load("output/models/per_gene/model_per_gene_out_ms.RData")
out_ms <- out_df

rm(out_df)

#subset based on model performance
pacman::p_load(qqman)

#subset by dispersion, outliers

out_lyre_sub <- subset(out_lyre, dispersion > 0.05 &
                         outlier > 0.05)

out_blue_sub <- subset(out_blue, dispersion > 0.05 &
                         outlier > 0.05)

out_eyec_sub <- subset(out_eyec, dispersion > 0.05 &
                         outlier > 0.05)

out_dist_sub <- subset(out_dist, dispersion > 0.05 &
                         outlier > 0.05)

out_fight_sub <- subset(out_fight, dispersion > 0.05 &
                         outlier > 0.05)

out_attend_sub <- subset(out_attend, dispersion > 0.05 &
                         outlier > 0.05)

## fdr correction
out_lyre_sub$mutations_qval <- p.adjust(out_lyre_sub$mutations_pval, method = "fdr", 
                                             n = nrow(out_lyre_sub))

out_blue_sub$mutations_qval <- p.adjust(out_blue_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_blue_sub))

out_eyec_sub$mutations_qval <- p.adjust(out_eyec_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_eyec_sub))

out_dist_sub$mutations_qval <- p.adjust(out_dist_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_dist_sub))

out_fight_sub$mutations_qval <- p.adjust(out_fight_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_fight_sub))

out_attend_sub$mutations_qval <- p.adjust(out_attend_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_attend_sub))

out_ms_sub$mutations_qval <- p.adjust(out_ms_sub$mutations_pval, method = "fdr", 
                                        n = nrow(out_ms_sub))

## attend is the only one with significant ones
summary(out_attend_sub$mutations_qval)
qq(out_attend_sub$mutations_pval)

sig1 <- gerp_n_by_gene_ls[[2448]] #LOC107052945
sig2 <- gerp_n_by_gene_ls[[3839]]#OXR1
sig3 <- gerp_n_by_gene_ls[[39]] #ACBD3, almost sig
sig4 <- gerp_n_by_gene_ls[[1395]] #FAM98A, almost sig

sig1 %>% left_join(pheno_long_models_ly_core, by = "id")%>%
  ggplot(aes(x = mutations_corrected, y = attend)) + geom_point() + 
  geom_smooth(method='lm')

sig2 %>% left_join(pheno_long_models_ly_core, by = "id")%>%
  ggplot(aes(x = mutations_corrected, y = attend)) + geom_point() + 
  geom_smooth(method='lm')


sig3 %>% left_join(pheno_long_models_ly_core, by = "id")%>%
  ggplot(aes(x = mutations_corrected, y = attend)) + geom_point() + 
  geom_smooth(method='lm')


sig4 %>% left_join(pheno_long_models_ly_core, by = "id")%>%
  ggplot(aes(x = mutations_corrected, y = attend)) + geom_point() + 
  geom_smooth(method='lm')

## export list of genes

out_attend_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_attend_list_genes_goanalysis.csv", quote=F, row.names = F, col.names = F)
out_lyre_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_lyre_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)
out_blue_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_blue_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)
out_eyec_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_eyec_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)
out_dist_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_dist_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)
out_fight_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_fight_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)
out_ms_sub %>% arrange(mutations_qval) %>% select(c(gene))%>%write.csv("output/models/per_gene/trait_ms_list_genes_goanalysis.csv", quote=F, row.names = F,col.names = F)

