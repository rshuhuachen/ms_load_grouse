#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=prep_zooroh
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/logs/prep_zooroh.err
#SBATCH --output=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/logs/prep_zooroh.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/scripts/inbreeding/1_zooroh.sh

# #get sample names
# vcfsamplenames output/ancestral/ltet_filtered_ann_aa.vcf.gz > output/inbreeding/adults_samples.txt
#  
# #convert vcf to geno format, exclude sex chr, recode 0.33 to 0 (missing data), rename scaffolds
# 
# bcftools convert -t ^'ScEsiA3_16870;HRSCAF=19426' -g output/ancestral/ltet_filtered_ann_aa_init --chrom --tag GT output/ancestral/ltet_filtered_ann_aa.vcf
# 
# zcat output/ancestral/ltet_filtered_ann_aa_init.gen.gz | \
#     sed -e 's/ 0.33/ 0/g' | \
#     sed 's/;HRSCAF=/__HRSCAF_/g' |
#     gzip > output/ancestral/ltet_filtered_ann_aa.gen.gz
#         

# run zooroh
Rscript scripts/inbreeding/run_zooroh.R \
    output/ancestral/ltet_filtered_ann_aa.gen.gz \
    output/inbreeding &> logs/run_zooroh.log

