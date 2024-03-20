#!/bin/bash
#$ -P fair_share
#$ -cwd
#$ -l idle=1
#$ -l vf=8G
#$ -N mcmcroh.sh
#$ -e /vol/cluster-data/rchen/wgr/analyses/rohcalling/runmcmcroh.err
#$ -o /vol/cluster-data/rchen/wgr/analyses/rohcalling/runmcmcroh.out
#$ -m e


bcftools roh -G30 --AF-dflt 0.4 data/genomic/intermediate/ltet_snps_filtered.vcf -o output/inbreeding/mcmc_roh.txt
