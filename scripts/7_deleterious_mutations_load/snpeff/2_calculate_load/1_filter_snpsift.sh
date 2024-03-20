#!/bin/bash
#$ -P fair_share
#$ -cwd
#$ -l idle=1
#$ -l vf=8G
#$ -N snpsift.sh
#$ -e log/snpsift_run.err
#$ -o log/snpsift_run.out
#$ -m e

### Filter for high and moderate impact mutations ####

## High impact
cat output/ancestral/ltet_filtered_ann_aa.vcf | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'HIGH' )" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf

## Moderate
cat output/ancestral/ltet_filtered_ann_aa.vcf | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'MODERATE')" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf

## LOF
cat output/ancestral/ltet_filtered_ann_aa.vcf | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'LOW') " > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf
