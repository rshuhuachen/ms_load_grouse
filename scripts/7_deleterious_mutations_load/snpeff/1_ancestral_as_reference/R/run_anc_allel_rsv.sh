#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=R_anc_all_scaf
#SBATCH --nodes=4
#SBATCH --error=scripts/snpeff/1_ancestral_as_reference/R/anc_all_scaf.err
#SBATCH --output=scripts/snpeff/1_ancestral_as_reference/R/anc_all_scaf.out

Rscript /vol/cluster-data/rchen/git/genetic_load_ltet/scripts/snpeff/1_ancestral_as_reference/R/anc_allele_tsv.R
