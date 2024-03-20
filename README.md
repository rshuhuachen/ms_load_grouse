# ms_load_inbreeding_grouse

This is the repository that contains all scripts used for the analysis in the manuscript titled XXX by Chen et al. (in progress). Genomic data, including the reference genome, the genome annotation and the 190 resequenced grouse genomes can be found on NCBI (see links below) and phenotypic data are stored in this repository. 

Below you will find an explanation of which data files can be found where, and the general structure of the workflow.

# Data

## Phenotypes
There are two files stored on GitHub/Dryad containing phenotypes: one for lifetime phenotypic traits (including lifetime mating success) which can be found under `data/phenotypes/phenotypes_lifetime` both in .csv and .RData format. The second file is in long format, where each row represents a caught male at a lek, one data entry per male per year. This file thus contains the sexual and behavioural traits for each male, and can be found under `data/phenotypes/phenotyes_annual` again in both .csv and .RData format.

## Reference genome
The reference genome can be found on NCBI XXX. Please download the file and store it in the folder `data/genomic/refgenome/` under name PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta or edit the scripts that use the reference genome according to the updated file name.

## Annotation
The genome annotation can be found on NCBI XXX. Please download the file and store it in the folder `data/genomic/annotation/` under name PO2979_Lyrurus_tetrix_black_grouse.annotation.gff or edit the scripts that use the annotation according to the updated file name.

## Resequencing
The raw resequencing data can be found on NCBI XXX. Please download the files and store it in the folder `data/genomic/raw/resequencing/`.

## HAL file
The publicly available 363-avian multiple alignment file can be downloaded from  https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020.hal and should be stored under `data/genomic/intermediate/cactus`. Note that this file is a few hundred Gb's in size!

# Scripts

The scripts used for this study can be found in chronological order under `scripts/`. Note that a large portion of the analyses are computationally intensive and are recommended to be run on a cluster. A brief explanation of each step, divided into different folders under `scripts/`, can be found below. Additionally, you can find a README.md in some of the separate folders for a more thorough explanation the more complex steps. Please note that pathnames should be edited according to your local device. Most processed files are not stored in this repository (such as .vcf files, GERP results) due to their large file sizes.

1_variant_calling: here we call genotypes from the raw whole genome resequencing files and filter the resulting .vcf file for quality as well as HWE, LD, MAF, etc.

2_cactus: here we manipulate the publicly available 363-avian multiple alignment file according to our needs, which is a file required to determine the ancestral state of the genome as well and for GERP++ to infer evolutionary conservation

3_sex_chr: here we align the chicken Z chromosome to our genome to infer which scaffolds are most likely part of the black grouse sex chromosome

4_pop_structure: here we conduct a PCA to inspect population structure within our population 

5_relatedness: here we assess relatedness coefficients based on genotypes in our population using NgsRelate

6_inbreeding: here we identify runs of homozygosity (ROHs) in the genome and calculate the realised inbreeding coefficient per individual

7_deleterious_mutations: here we calculate deleteriousness of mutations using two complementary approaches: 1) using GERP++ to estimate evolutionary conservation and 2) using SnpEff to estimate the effect of a mutation on protein function. We also compare the overlap between the two approaches.

8_genetic_load: based on the two approaches to infer deleteriousness of mutations, we can then calculate a proxy for total, expressed and masked load. Within this folder, we repeat the above per gene region and explore the relationship between the different approaches

9_models: here we use Bayesian mixed effect modelling to understand the relationship between genetic quality (inbreeding, genetic load) and lifetime mating success (LMS), as well as their effects on sexual and behavioural traits

10_plots: in this folder, we produce the plots used in the main manuscript as well as in the supplementary materials