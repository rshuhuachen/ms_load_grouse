#!/bin/bash
#$ -P fair_share
#$ -cwd
#$ -l idle=1
#$ -l vf=8G
#$ -N snpsift
#$ -e log/snpsift_run.err
#$ -o log/snpsift_run.out
#$ -m e

### Filter for high and moderate impact mutations ####

## High impact
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'HIGH' )" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_HIGH.vcf

## Moderate
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'MODERATE')" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_moderate.vcf

## Low
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'LOW') " > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_low.vcf

## LOF
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " (exists LOF[*].PERC )" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_LOF.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_LOF.vcf

## Missense
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " ANN[*].EFFECT has 'missense_variant'" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_missense.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_missense.vcf

## Synonymous
cat output/ancestral/ltet_filtered_ann_aa.vcf | java -jar src/SnpSift.jar filter "( ANN[*].EFFECT has 'synonymous_variant' )" > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_synonymous.vcf
gzip  data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_synonymous.vcf

## Non-synonymous: missense, nonsense, NMD
cat output/ancestral/ltet_filtered_ann_aa.vcf | java -jar src/SnpSift.jar filter "ANN[*].IMPACT = 'MODERATE' | ANN[*].IMPACT = 'MODIFIER' | ANN[*].IMPACT = 'HIGH' " > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_non_synonymous.vcf
gzip  data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_non_synonymous.vcf

## Modifier
zcat output/ancestral/ltet_filtered_ann_aa.vcf.gz | java -jar src/SnpSift.jar filter " ( ANN[*].IMPACT = 'MODIFIER') " > data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_modifier.vcf
gzip data/genomic/intermediate/snpef/ltet_ann_aa_snp_output_modifier.vcf

