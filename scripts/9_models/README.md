# Modelling

In this sub-directory, we build various Bayesian models to answer whether deleterious mutations affect (lifetime) reproductive success, and whether this is mediated through the sexual traits. 

In the `lms` folder, we use the various types of genetic load as predictor variables and lifetime mating success (LMS) as the response variable. 

In the `traits` folder, we model the effect of genetic load on the various sexual and behavioural traits, and the effect of the traits on annual mating success. 

As both of these modelling approaches are repetitive for the various types of load (GERP vs SnpEff, different GERP categories,  all regions vs specific gene regions), these models were implemented using Snakemake to iterate over similar model structures. The R scripts are called within the snakemake file. The models can be executed using the snakemake calls given at the start of the snakefile (e.g. `snakemake --snakefile scripts/8_models/annual/snakefile_megamodels_traits`, but ideally this is done on a cluster).

First, the models have to be run (done in the separate `annual` and `lms` folders) and then additional scripts can be run. First, we calculate allele frequencies of the derived deleterious alleles (`1_allele_frequencies.R`, doesn't use load though) and of derived neutral alleles (`1b_allele_frequencies_neutral.R`), then we loop over the model output to diagnose the models (e.g. calculate Rhat, plot the Neff etc.) in the script `2_loop_diagnose_models`. For the annual models, we use multiple model outputs to calculate direct and indirect effects of load on mating success. This is executed in the script `3_direct_indirect.R`. In addition to the model diagnosis, we tested for prior sensitivity in the fourth script (`4_test_prior_sensitivity.R`). Note that only the lms models are tested for prior sensitivity in this script, however all models were tested using a lower number of iterations using the same steps. 
