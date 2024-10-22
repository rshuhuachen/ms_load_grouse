
In this script, we calculate mutation load based on GERP and SnpEff deleterious mutations.

First, we filter for various SnpEff impact classess in the script `1_filter_snpeff_snpsift.sh`. Then, we calculate load within the R script `2_calculate_load_snpeff_gerp.R` which calls the functions defined in the script `function_calculate_load.R`. 

## Gene regions

As we are interested in the effect of mutations across different genetic regions, we next annotate deleterious mutations with gene regions using the genome annotation in the file `3_load_per_gene_region.R`. First, we split the gff file in different gene regions, followed by the annotation of each deleterious mutation. Lastly, we calculate genetic load based on the subset of mutations in a specific gene regions (e.g. TSS, intron, exon, promoter).

## Combine

Lastly, we combine all genetic load scores in one file, which is outputted in `output/load/all_loads_combined_da_nosex_29scaf_plus_per_region.RData`. 