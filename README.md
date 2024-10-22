# Predicted deleterious mutations reveal the genomic mechanisms underlying fitness variation in a lekking bird

# General

This is the repository that contains all scripts used for the analysis in the manuscript titled "Predicted deleterious mutations reveal the genomic mechanisms underlying fitness variation in a lekking bird" by Chen et al. (in preparation). Genomic data, including the reference genome, and the 190 resequenced grouse genomes can be found on NCBI (see links below) and phenotypic data as well as the black grouse annotation are stored in this repository.

Below you will find an explanation of which data files can be found where, and the general structure of the workflow.
You will find a brief overview of the scripts with an explanation here: https://rshuhuachen.github.io/ms_load_grouse/

# Data

## Phenotypes

There are two files stored on GitHub containing phenotypes: one for lifetime phenotypic traits (including lifetime mating success) which can be found under `data/phenotypes/phenotypes_lifetime` both in .csv and .RData format. 
The second file is in long format, where each row represents a caught male at a lek, one data entry per male per year. This file contains the sexual and behavioural traits for each male, as well as annual mating success, and can be found under `data/phenotypes/phenotyes_annual` again in both .csv and .RData format.

## Reference genome

The reference genome can be found on NCBI BioProject PRJNA1085187. Please download the file and store it in the folder `data/genomic/refgenome/` under name PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta or edit the scripts that use the reference genome according to the updated file name.

## Annotation

The genome annotation can be found in the folder `data/genomic/annotation/` under name PO2979_Lyrurus_tetrix_black_grouse.annotation.gff.

## Resequencing

The raw resequencing data can be found on NCBI BioProject PRJNA1085187 under SRA BioAccession Numbers 40722954 - 40723143. Please download the files and store it in the folder `data/genomic/raw/resequencing/`.

## HAL file

The publicly available 363-avian multiple alignment file used for ancestral genome reconstruction and calculating GERP scores can be downloaded from https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020.hal and should be stored under `data/genomic/intermediate/cactus`. Note that this file is a few hundred Gb's in size!

# Scripts

The scripts used for this study can be found in chronological order under `scripts/`. Note that a large portion of the analyses are computationally intensive and are recommended to be run on a cluster. Several scripts include headers for submitting through a job scheduler, such as SLURM/Sun Grid Engine which should be edited according to your local cluster/infrastructure.

A brief explanation of each step, divided into different folders under `scripts/`, can be found below. Additionally, you can find a README.md in some of the separate folders for a more thorough explanation of the more complex steps. Please note that pathnames should be edited according to your local device. R scripts should be opened from the .Rproj to allow easier work within the current directory, as the Rstudio Project will assume the root working directory as the default location. Note that some but not all processed files are stored in this repository (e.g. excluding .vcf files and GERP results) due to their large file sizes. However, all data used for plots are stored within this repository.

To allow easier reproducibility, you can use the conda environment listed in `src/envs/environment_load.yml`. Ensure you have conda installed on your local device, and then use `conda env create --name load --file src/envs/environment_load.yml` to install all the software and their correct versions as used in this manuscript. Ensure to activate the environment (`conda activate load`) before starting the analyses.

## Brief explanation of each script's subdirectory

1_variant_calling: here we align the sequences to the reference genome and call SNPs. We next filter the resulting .vcf file for quality as well as HWE, LD, MAF, etc. 

2_cactus: here we manipulate the publicly available 363-avian multiple alignment file according to our needs (i.e. reducing the phylogenetic tree and adding focal genomes) which is a file required to determine the ancestral state of the genome and for GERP++ to infer evolutionary conservation

3_sex_chr: here we align the chicken Z chromosome to our genome to infer which scaffolds are most likely part of the black grouse sex chromosome

4_pop_structure: here we conduct a PCA to inspect population structure within our population

5_relatedness: here we assess relatedness coefficients based on genotypes in our population using NgsRelate and plink

6_snpeff_gerp: here we use two variant prioritisation tools to infer deleterious mutations: we run snpeff (first build a data base and then we polarise the genome to recode 1 as derived and 0 as ancestral) and gerp (first we run gerp across the genome and then select those that overlap with snps). Note the analyses should be performed with snpeff first, as the polarised genome is also used for the gerp - snp intersection to result in a vcf file with both snpeff annotations, gerp scores, and ancestral-derived recoded genotypes.

7_calculate_load: here we use the snpeff and gerp results to filter for the most deleterious mutations and calculated a proxy for total load, as well as homozygous and heterozygous load. Within this folder, we repeat the above per gene region (promoter, TSS, intron, exon) using the genome annotation.

8_models: here we use Bayesian mixed effect modelling to understand the relationship between mutation load and lifetime mating success (LMS), as well as their effects on sexual and behavioural traits. We also calculate allele frequencies in this directory

9_random_draws: in this subfolder, we use a subsampling approach to select a random set of mutations to control for the number of mutations that contribute to the mutation load estimate

10_plots: in this folder, we produce the plots used in the main manuscript as well as in the supplementary materials

Lastly, please see our full Methods for further details.
