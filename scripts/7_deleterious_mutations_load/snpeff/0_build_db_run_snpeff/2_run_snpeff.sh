#!/bin/bash
#$ -P fair_share
#$ -cwd
#$ -l idle=1
#$ -l vf=8G
#$ -N snpeff.sh
#$ -e log/snpeff_run.err
#$ -o log/snpeff_run.out
#$ -m e

java -Xmx8g -jar src/snpEff.jar ann -stats  \
  -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -v -formatEff \
  lyrurus_tetrix data/genomic/intermediate/ltet_snps_filtered.vcf \
  > data/genomic/intermediate/snpef/ltet_ann_snp_output.vcf
  


