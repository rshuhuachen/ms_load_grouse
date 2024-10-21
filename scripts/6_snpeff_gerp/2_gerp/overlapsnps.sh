#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=overlap
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/log/overlap.err
#SBATCH --output=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/log/overlap.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/scripts/6_snpeff_gerp/2_gerp/overlapsnps.sh


for i in {4..30};
    do bedtools intersect -a output/gerp/beds/gerp_scaf_$i.bed -b output/ancestral/ltet_filtered_ann_aa.bed -wa -wb | cut -f 6-10 --complement > output/gerp/beds/gerp_overlapSNP_scaf_$i.tsv
    gzip output/gerp/beds/gerp_overlapSNP_scaf_$i.tsv;
    done
