#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=countgerp2
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/log/countgerp2.err
#SBATCH --output=/vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/log/countgerp2.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de
# launch with: sbatch /vol/cluster-data/rchen/git/ms_load_inbreeding_grouse/scripts/6_snpeff_gerp/2_gerp/count_gerp.sh


for i in {5..30};
    do 
    Rscript --vanilla scripts/6_snpeff_gerp/2_gerp/count_gerp.R output/gerp/beds/gerp_overlapSNP_scaf_$i.tsv.gz $i output/gerp/count_mutations_gerp_$i.tsv 
        
    done
