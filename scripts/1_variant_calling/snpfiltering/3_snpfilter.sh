#!/bin/bash
#$ -P fair_share
#$ -cwd
#$ -l idle=1
#$ -l vf=8G
#$ -N snpfilter_wgr.sh
#$ -e logs/snpfilter_wgr.err
#$ -o logs/snpfilter_wgr.out
#$ -m e

### Specify all parameters ###

## filtering indels
minallele=2
maxallele=2

## min DP, max missing, min QC, max DP
minDP=20
maxmissing=0.7
maxDP=60 #=2x mean DP
minQ=30
set="minDP20_maxmis0.7_maxdp60_minQC30"

## one liner filtering step

# make sure to change the paths where you have the unfiltered vcf file as well as the filtered one


vcftools --vcf data/genomic/intermediate/rawSNPcalls.vcf \
    -min-alleles 2 --max-alleles 2 \
    --minDP 20 --max-missing 0.7 \
    --max-meanDP 60 --minQ 30 \
    --recode --stdout | gzip -c > data/genomic/intermediate/ltet_snps_filtered.vcf
