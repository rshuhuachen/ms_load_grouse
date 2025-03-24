# Genetic architecture of male reproductive success in a lekking bird: insights from predicted deleterious mutations

# General

This is the repository that contains all scripts used for the analysis in the manuscript titled "Genetic architecture of male reproductive success in a lekking bird: insights from predicted deleterious mutations", Chen et al. 2024 ([in review](https://doi.org/10.21203/rs.3.rs-5579350/v1)). Genomic data, including the reference genome, and the 190 resequenced grouse genomes can be found on NCBI (see links below) and phenotypic data as well as the black grouse annotation are stored in this repository.

Below you will find an explanation of which data files can be found where, and the general structure of the workflow. You will find a brief overview of the scripts with an explanation here: https://rshuhuachen.github.io/ms_load_grouse/

# Data

## Phenotypes

There are two files stored on GitHub containing phenotypes: one for lifetime phenotypic traits (including lifetime mating success) which can be found under `data/phenotypes/phenotypes_lifetime` both in .csv and .RData format. The second file is in long format, where each row represents a caught male at a lek, one data entry per male per year. This file contains the sexual and behavioural traits for each male, as well as annual mating success, and can be found under `data/phenotypes/phenotyes_annual` again in both .csv and .RData format.

## Reference genome

The reference genome can be found on NCBI BioProject PRJNA1085187. Please download the file and store it in the folder `data/genomic/refgenome/` under name PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta or edit the scripts that use the reference genome according to the updated file name.

## Annotation

The genome annotation can be found in the folder `data/genomic/annotation/` under name PO2979_Lyrurus_tetrix_black_grouse.annotation.gff. The RNA seq data generated to annotate the genome can be found under NCBI Bio Accession [SRX24353984](https://www.ncbi.nlm.nih.gov/sra/SRX24353984%5Baccn%5D%5D)

## Resequencing

The raw resequencing data can be found on NCBI BioProject [PRJNA1085187](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1085187) under SRA study [SRP499251](https://trace.ncbi.nlm.nih.gov/Traces/?view=study&acc=SRP499251) with BioAccession Numbers SRR28526036 – SRR28526225. Please download the files and store it in the folder `data/genomic/raw/resequencing/`.

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

8_inbreeding: here, we calculate genomic inbreeding using BCFtools where we first identify runs of homozygosity, and then calculate FROH. We further model the effect of FROH on lifetime mating success to quantify inbreeding depression.

9_models: here we use Bayesian mixed effect modelling to understand the relationship between mutation load and lifetime mating success (LMS), as well as their effects on sexual and behavioural traits. We also calculate allele frequencies in this directory

10_biological_processes: here we utilise gene ontology annotations to select genes that are important for hypothesized relevant biological processes. We calculate the total load based on subsets of mutations and model their effect on LMS.

11_random_draws: in this subfolder, we use a subsampling approach to select a random set of mutations to control for the number of mutations that contribute to the mutation load estimate

12_plots: in this folder, we produce the plots used in the main manuscript as well as in the supplementary materials

Lastly, please see our full Methods in the main manuscript for further details.

## Dependencies
All command-line based software used and their versions can be found in the file `src/envs/environment_load.yml`
All R packages and their versions can be found in the main manuscript (Methods section) or alternatively, see the list below:

R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] tolerance_3.0.0        ggridges_0.5.6         cowplot_1.1.3         
 [4] rtracklayer_1.64.0     GenomicFeatures_1.56.0 AnnotationDbi_1.66.0  
 [7] Biobase_2.64.0         GenomicRanges_1.56.1   GenomeInfoDb_1.40.1   
[10] IRanges_2.38.1         S4Vectors_0.42.1       BiocGenerics_0.50.0   
[13] bayesplot_1.11.1.9000  brms_2.21.0            Rcpp_1.0.13           
[16] DHARMa_0.4.6           lmerTest_3.1-3         lme4_1.1-35.5         
[19] Matrix_1.7-0           glmmTMB_1.1.10         data.table_1.16.0     
[22] lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1         
[25] dplyr_1.1.4            purrr_1.0.2            readr_2.1.5           
[28] tidyr_1.3.1            tibble_3.2.1           ggplot2_3.5.1         
[31] tidyverse_2.0.0       

loaded via a namespace (and not attached):
  [1] tensorA_0.36.2.1            rstudioapi_0.16.0          
  [3] jsonlite_1.8.8              magrittr_2.0.3             
  [5] TH.data_1.1-2               estimability_1.5.1         
  [7] nloptr_2.1.1                BiocIO_1.14.0              
  [9] zlibbioc_1.50.0             vctrs_0.6.5                
 [11] memoise_2.0.1               minqa_1.2.8                
 [13] Rsamtools_2.20.0            RCurl_1.98-1.16            
 [15] htmltools_0.5.8.1           S4Arrays_1.4.1             
 [17] distributional_0.4.0        curl_5.2.2                 
 [19] SparseArray_1.4.8           StanHeaders_2.35.0.9000    
 [21] htmlwidgets_1.6.4           sandwich_3.1-0             
 [23] plotly_4.10.4               emmeans_1.10.4             
 [25] zoo_1.8-12                  cachem_1.1.0               
 [27] TMB_1.9.15                  GenomicAlignments_1.40.0   
 [29] lifecycle_1.0.4             pkgconfig_2.0.3            
 [31] R6_2.5.1                    fastmap_1.2.0              
 [33] MatrixGenerics_1.16.0       GenomeInfoDbData_1.2.12    
 [35] rbibutils_2.2.16            digest_0.6.37              
 [37] numDeriv_2016.8-1.1         colorspace_2.1-1           
 [39] RSQLite_2.3.7               fansi_1.0.6                
 [41] timechange_0.3.0            httr_1.4.7                 
 [43] abind_1.4-8                 mgcv_1.9-1                 
 [45] compiler_4.4.1              bit64_4.0.5                
 [47] withr_3.0.1                 backports_1.5.0            
 [49] inline_0.3.19               BiocParallel_1.38.0        
 [51] DBI_1.2.3                   QuickJSR_1.3.1             
 [53] pkgbuild_1.4.4              MASS_7.3-61                
 [55] DelayedArray_0.30.1         rjson_0.2.22               
 [57] loo_2.8.0.9000              tools_4.4.1                
 [59] glue_1.7.0                  restfulr_0.0.15            
 [61] nlme_3.1-166                grid_4.4.1                 
 [63] checkmate_2.3.2             generics_0.1.3             
 [65] gtable_0.3.5                tzdb_0.4.0                 
 [67] hms_1.1.3                   utf8_1.2.4                 
 [69] XVector_0.44.0              pillar_1.9.0               
 [71] posterior_1.6.0             splines_4.4.1              
 [73] lattice_0.22-6              survival_3.7-0             
 [75] bit_4.0.5                   tidyselect_1.2.1           
 [77] Biostrings_2.72.1           reformulas_0.3.0           
 [79] gridExtra_2.3               V8_5.0.0                   
 [81] SummarizedExperiment_1.34.0 bridgesampling_1.1-2       
 [83] matrixStats_1.4.1           rstan_2.35.0.9000          
 [85] stringi_1.8.4               UCSC.utils_1.0.0           
 [87] lazyeval_0.2.2              yaml_2.3.10                
 [89] pacman_0.5.1                boot_1.3-31                
 [91] codetools_0.2-20            cli_3.6.3                  
 [93] RcppParallel_5.1.9          xtable_1.8-4               
 [95] Rdpack_2.6.1                munsell_0.5.1              
 [97] coda_0.19-4.1               png_0.1-8                  
 [99] XML_3.99-0.17               parallel_4.4.1             
[101] rstantools_2.4.0.9000       blob_1.2.4                 
[103] Brobdingnag_1.2-9           bitops_1.0-8               
[105] viridisLite_0.4.2           mvtnorm_1.3-0              
[107] scales_1.3.0                crayon_1.5.3               
[109] rlang_1.1.4                 KEGGREST_1.44.1            
[111] multcomp_1.4-26  
