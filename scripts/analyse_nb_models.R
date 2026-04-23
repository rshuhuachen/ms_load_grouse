##### Import all models #####
pacman::p_load(bayesplot, brms, dplyr, data.table, extraDistr)

#### Total load #####
###### GERP ######
load("output/models/total_hom_het/lms_total_gerp45.RData")
m_gerp <- brm_load_t
rm(brm_load_t)

load("output/models/total_hom_het/lms_total_gerp45_nb.RData")
m_gerp_nb <- brm_load_t
rm(brm_load_t)

load("output/models/total_hom_het/lms_total_gerp45_zinb.RData")
m_gerp_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_gerp)
summary(m_gerp_nb)
summary(m_gerp_zinb)

mcmc_intervals_data(m_gerp, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_gerp, type = "bars", ndraws=100)
pp_check(m_gerp_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_zinb, type = "bars", ndraws = 100)

# plot
mcmc_areas(m_gerp,pars=c("b_scaletotal_load"))

# compare
m_gerp <- add_criterion(m_gerp, c("loo"))
m_gerp_nb <- add_criterion(m_gerp_nb, c("loo"))
m_gerp_zinb <- add_criterion(m_gerp_zinb, c("loo"))
loo_compare(m_gerp, m_gerp_nb, m_gerp_zinb, criterion="loo")

###### SnpEff ######
load("output/models/total_hom_het/lms_total_high.RData")
m_high <- brm_load_t
rm(brm_load_t)

load("output/models/total_hom_het/lms_total_high_nb.RData")
m_high_nb <- brm_load_t
rm(brm_load_t)

load("output/models/total_hom_het/lms_total_high_zinb.RData")
m_high_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_high)
summary(m_high_nb)
summary(m_high_zinb)

# intervals
mcmc_intervals_data(m_high, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_high, type = "bars", ndraws=100)
pp_check(m_high_nb, type = "bars", ndraws = 100)
pp_check(m_high_zinb, type = "bars", ndraws = 100)

# compare
m_high <- add_criterion(m_high, c("loo"))
m_high_nb <- add_criterion(m_high_nb, c("loo"))
m_high_zinb <- add_criterion(m_high_zinb, c("loo"))
loo_compare(m_high, m_high_nb, m_high_zinb, criterion="loo")

##### Hom and het ######
###### GERP ######
load("output/models/total_hom_het/lms_het_hom_gerp45.RData")
m_gerp_hethom <- brm_load_het_hom
rm(brm_load_het_hom)

load("output/models/total_hom_het/lms_het_hom_gerp45_nb.RData")
m_gerp_hethom_nb <- brm_load_het_hom
rm(brm_load_het_hom)

load("output/models/total_hom_het/lms_het_hom_gerp45_zinb.RData")
m_gerp_hethom_zinb <- brm_load_het_hom
rm(brm_load_het_hom)

## compare
summary(m_gerp_hethom)
summary(m_gerp_hethom_nb)
summary(m_gerp_hethom_zinb)

mcmc_intervals_data(m_gerp_hethom, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))
mcmc_intervals_data(m_gerp_hethom_nb, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))
mcmc_intervals_data(m_gerp_hethom_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))

# fit raw data to model
pp_check(m_gerp_hethom, type = "bars", ndraws=100)
pp_check(m_gerp_hethom_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_hethom_zinb, type = "bars", ndraws = 100)

# compare
m_gerp_hethom <- add_criterion(m_gerp_hethom, c("loo"))
m_gerp_hethom_nb <- add_criterion(m_gerp_hethom_nb, c("loo"))
m_gerp_hethom_zinb <- add_criterion(m_gerp_hethom_zinb, c("loo"))
loo_compare(m_gerp_hethom, m_gerp_hethom_nb, m_gerp_hethom_zinb, criterion="loo")

###### SnpEff ######
load("output/models/total_hom_het/lms_het_hom_high.RData")
m_high_hethom <- brm_load_het_hom
rm(brm_load_het_hom)

load("output/models/total_hom_het/lms_het_hom_high_nb.RData")
m_high_hethom_nb <- brm_load_het_hom
rm(brm_load_het_hom)

load("output/models/total_hom_het/lms_het_hom_high_zinb.RData")
m_high_hethom_zinb <- brm_load_het_hom
rm(brm_load_het_hom)

## compare
summary(m_high_hethom)
summary(m_high_hethom_nb)
summary(m_high_hethom_zinb)

mcmc_intervals_data(m_high_hethom, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))
mcmc_intervals_data(m_high_hethom_nb, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))
mcmc_intervals_data(m_high_hethom_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scalehet_load", "b_scalehom_load"))

# fit raw data to model
pp_check(m_high_hethom, type = "bars", ndraws=100)
pp_check(m_high_hethom_nb, type = "bars", ndraws = 100)
pp_check(m_high_hethom_zinb, type = "bars", ndraws = 100)

# compare
m_high_hethom <- add_criterion(m_high_hethom, c("loo"))
m_high_hethom_nb <- add_criterion(m_high_hethom_nb, c("loo"))
m_high_hethom_zinb <- add_criterion(m_high_hethom_zinb, c("loo"))
loo_compare(m_high_hethom, m_high_hethom_nb, m_high_hethom_zinb, criterion="loo")

##### Regions ######
###### GERP ######
######## Exons ########
load("output/models/per_gene_region/lms_total_gerp45_exon.RData")
m_gerp_exon <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_exon_nb.RData")
m_gerp_exon_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_exon_zinb.RData")
m_gerp_exon_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_gerp_exon)
summary(m_gerp_exon_nb)
summary(m_gerp_exon_zinb)

mcmc_intervals_data(m_gerp_exon, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_exon_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_exon_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_gerp_exon, type = "bars", ndraws=100)
pp_check(m_gerp_exon_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_exon_zinb, type = "bars", ndraws = 100)

# compare
m_gerp_exon <- add_criterion(m_gerp_exon, c("loo"))
m_gerp_exon_nb <- add_criterion(m_gerp_exon_nb, c("loo"))
m_gerp_exon_zinb <- add_criterion(m_gerp_exon_zinb, c("loo"))
loo_compare(m_gerp_exon, m_gerp_exon_nb, m_gerp_exon_zinb, criterion="loo")

######## promoters ########
load("output/models/per_gene_region/lms_total_gerp45_promoter.RData")
m_gerp_promoter <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_promoter_nb.RData")
m_gerp_promoter_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_promoter_zinb.RData")
m_gerp_promoter_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_gerp_promoter)
summary(m_gerp_promoter_nb)
summary(m_gerp_promoter_zinb)

mcmc_intervals_data(m_gerp_promoter, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_promoter_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_promoter_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_gerp_promoter, type = "bars", ndraws=100)
pp_check(m_gerp_promoter_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_promoter_zinb, type = "bars", ndraws = 100)

# compare
m_gerp_promoter <- add_criterion(m_gerp_promoter, c("loo"))
m_gerp_promoter_nb <- add_criterion(m_gerp_promoter_nb, c("loo"))
m_gerp_promoter_zinb <- add_criterion(m_gerp_promoter_zinb, c("loo"))
loo_compare(m_gerp_promoter, m_gerp_promoter_nb, m_gerp_promoter_zinb, criterion="loo")

######## tss ########
load("output/models/per_gene_region/lms_total_gerp45_tss.RData")
m_gerp_tss <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_tss_nb.RData")
m_gerp_tss_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_tss_zinb.RData")
m_gerp_tss_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_gerp_tss)
summary(m_gerp_tss_nb)
summary(m_gerp_tss_zinb)

mcmc_intervals_data(m_gerp_tss, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_tss_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_tss_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_gerp_tss, type = "bars", ndraws=100)
pp_check(m_gerp_tss_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_tss_zinb, type = "bars", ndraws = 100)

# compare
m_gerp_tss <- add_criterion(m_gerp_tss, c("loo"))
m_gerp_tss_nb <- add_criterion(m_gerp_tss_nb, c("loo"))
m_gerp_tss_zinb <- add_criterion(m_gerp_tss_zinb, c("loo"))
loo_compare(m_gerp_tss, m_gerp_tss_nb, m_gerp_tss_zinb, criterion="loo")

######## intron ########
load("output/models/per_gene_region/lms_total_gerp45_intron.RData")
m_gerp_intron <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_intron_nb.RData")
m_gerp_intron_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_gerp45_intron_zinb.RData")
m_gerp_intron_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_gerp_intron)
summary(m_gerp_intron_nb)
summary(m_gerp_intron_zinb)

mcmc_intervals_data(m_gerp_intron, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_intron_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_gerp_intron_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_gerp_intron, type = "bars", ndraws=100)
pp_check(m_gerp_intron_nb, type = "bars", ndraws = 100)
pp_check(m_gerp_intron_zinb, type = "bars", ndraws = 100)

# compare
m_gerp_intron <- add_criterion(m_gerp_intron, c("loo"))
m_gerp_intron_nb <- add_criterion(m_gerp_intron_nb, c("loo"))
m_gerp_intron_zinb <- add_criterion(m_gerp_intron_zinb, c("loo"))
loo_compare(m_gerp_intron, m_gerp_intron_nb, m_gerp_intron_zinb, criterion="loo")

###### SnpEff ######
######## Exons ########
load("output/models/per_gene_region/lms_total_high_exon.RData")
m_high_exon <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_exon_nb.RData")
m_high_exon_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_exon_zinb.RData")
m_high_exon_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_high_exon)
summary(m_high_exon_nb)
summary(m_high_exon_zinb)

mcmc_intervals_data(m_high_exon, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_exon_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_exon_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_high_exon, type = "bars", ndraws=100)
pp_check(m_high_exon_nb, type = "bars", ndraws = 100)
pp_check(m_high_exon_zinb, type = "bars", ndraws = 100)

# compare
m_high_exon <- add_criterion(m_high_exon, c("loo"))
m_high_exon_nb <- add_criterion(m_high_exon_nb, c("loo"))
m_high_exon_zinb <- add_criterion(m_high_exon_zinb, c("loo"))
loo_compare(m_high_exon, m_high_exon_nb, m_high_exon_zinb, criterion="loo")

######## promoters ########
load("output/models/per_gene_region/lms_total_high_promoter.RData")
m_high_promoter <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_promoter_nb.RData")
m_high_promoter_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_promoter_zinb.RData")
m_high_promoter_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_high_promoter)
summary(m_high_promoter_nb)
summary(m_high_promoter_zinb)

mcmc_intervals_data(m_high_promoter, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_promoter_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_promoter_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_high_promoter, type = "bars", ndraws=100)
pp_check(m_high_promoter_nb, type = "bars", ndraws = 100)
pp_check(m_high_promoter_zinb, type = "bars", ndraws = 100)

# compare
m_high_promoter <- add_criterion(m_high_promoter, c("loo"))
m_high_promoter_nb <- add_criterion(m_high_promoter_nb, c("loo"))
m_high_promoter_zinb <- add_criterion(m_high_promoter_zinb, c("loo"))
loo_compare(m_high_promoter, m_high_promoter_nb, m_high_promoter_zinb, criterion="loo")

######## tss ########
load("output/models/per_gene_region/lms_total_high_tss.RData")
m_high_tss <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_tss_nb.RData")
m_high_tss_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_tss_zinb.RData")
m_high_tss_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_high_tss)
summary(m_high_tss_nb)
summary(m_high_tss_zinb)

mcmc_intervals_data(m_high_tss, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_tss_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_tss_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_high_tss, type = "bars", ndraws=100)
pp_check(m_high_tss_nb, type = "bars", ndraws = 100)
pp_check(m_high_tss_zinb, type = "bars", ndraws = 100)

# compare
m_high_tss <- add_criterion(m_high_tss, c("loo"))
m_high_tss_nb <- add_criterion(m_high_tss_nb, c("loo"))
m_high_tss_zinb <- add_criterion(m_high_tss_zinb, c("loo"))
loo_compare(m_high_tss, m_high_tss_nb, m_high_tss_zinb, criterion="loo")

######## intron ########
load("output/models/per_gene_region/lms_total_high_intron.RData")
m_high_intron <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_intron_nb.RData")
m_high_intron_nb <- brm_load_t
rm(brm_load_t)

load("output/models/per_gene_region/lms_total_high_intron_zinb.RData")
m_high_intron_zinb <- brm_load_t
rm(brm_load_t)

## compare
summary(m_high_intron)
summary(m_high_intron_nb)
summary(m_high_intron_zinb)

mcmc_intervals_data(m_high_intron, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_intron_nb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))
mcmc_intervals_data(m_high_intron_zinb, prob_outer=0.95, prob=0.9, pars=c("b_scaletotal_load"))

# fit raw data to model
pp_check(m_high_intron, type = "bars", ndraws=100)
pp_check(m_high_intron_nb, type = "bars", ndraws = 100)
pp_check(m_high_intron_zinb, type = "bars", ndraws = 100)

# compare
m_high_intron <- add_criterion(m_high_intron, c("loo"))
m_high_intron_nb <- add_criterion(m_high_intron_nb, c("loo"))
m_high_intron_zinb <- add_criterion(m_high_intron_zinb, c("loo"))
loo_compare(m_high_intron, m_high_intron_nb, m_high_intron_zinb, criterion="loo")
