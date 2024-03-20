#!/bin/bash
#SBATCH --account=fair_share
#SBATCH --mem=8000
#SBATCH --job-name=blast
#SBATCH --nodes=10
#SBATCH --error=/vol/cluster-data/rchen/git/genetic_load_ltet/log/blast_ggal.err
#SBATCH --output=/vol/cluster-data/rchen/git/genetic_load_ltet/log/blast_ggal.out
#SBATCH --mail-type=END
#SBATCH --mail-user=rebecca.chen@uni-bielefeld.de

### use liftoff instead

liftoff data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta data/genomic/annotation/gallus_gallus/ncbi_dataset/data/GCF_016699485.2/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna -g data/genomic/annotation/gallus_gallus/ncbi_dataset/data/GCF_016699485.2/genomic.gff -o data/genomic/annotation/gallus_gallus/liftoff_gallus_ltet.gff -u data/genomic/annotation/gallus_gallus/unmapped_features.txt -p 4
