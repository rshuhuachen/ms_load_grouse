# Editing the 363-avian file 

In this repository, we edit the 363-avian multi-alignment file to our specific needs/research questions.

First, we have to reduce the avian multiple alignment file to exclude species of the Neoaves clade (1_remove_subtrees.sh). 

Next, we add two genomes: *Lyrurs tetrix* and *Lagopus lecura*, using the `cactus prepare` function (2_cactus_prepare.txt). To add the *Lagopus lecura* genome to the hal file, the assembly needs to be downloaded from NCBI which can be found on NCBI: https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/datasets/genome/GCF_019238085.1/. This file should be placed under: `data/genomic/intermediate/cactus`. The two files `input_ltet.txt` and `input_lleu.txt` specify the names and file locations of the fasta files of both assemblies.

The two files that get outputted by the cactus preparation step can then be run in sequence to add the new genomes to the phylogenetic tree (` 3_cactus_update_lyrurus_steps.sh` and `4_cactus_update_lagopus_steps.sh`).

For subsequent steps, the resulting hal file is converted to wiggle format (`5_hal_to_wiggle.R`) to investigate the coverage per scaffold, and in maf format per scaffold (`6_split_hal_per_scaf.R`) for both GERP++ and the neutral tree calculation. Note that we only focus on the 30 largest scaffolds of the black grouse genome which over >95% of the genome.

Lastly, the final phylogenetic tree has to be recalculated according to the updated tree, and the branch lengths have to be calculated in substitutions/site (rather than million years ago). The phylogenetic tree from the edited hal file can be found in `output/cactus/topology.tree` which can be generated with the halStats command. The `snakefile_neutraltree` is a snakemake file containing the workflow to estimate the updated branch lengths with the added genomes, and uses the R script `collapse_bed_coverage.R` as well as the `maffilter_templ.txt` file. The resulting updated phylogenetic tree can be found in `output/cactus/combined_windows.fa.treefile`.