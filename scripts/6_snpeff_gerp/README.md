# Inferring deleterious mutations

In this directory, we infer deleterious mutations using GERP++ and SnpEff and calculate genetic load.

The reconstructed genome of the last common ancestor between *Lyrurus tetrix* and *Lagopus leucura* is used to infer the ancestral 'wild-type' state of the genome. Although the manuscript presents the GERP scores before the SnpEff results, in the scripts this is done in reverse as the SnpEff annotated vcf file was used to then append GERP scores, ensuring all results were attached to a single genotype file rather than working with multiple files. Therefore, the workflow in this sub-directory is explained starting with SnpEff.

## SnpEff

In the `scripts/6_snpeff_gerp/snpeff` directory, you can find the workflow to build a SnpEff database from our genome annotation and subsequently, run SnpEff and SnpSift to filter. Make sure to install SnpEff from source https://pcingola.github.io/SnpEff/download/ and put the source code under `src/`. First, the database is built and SnpEff is run `0_build_db_run_snpeff`. Then, using the generated .vcf as input, the ancestral allele is inferred (`1_ancestral_as_reference`). The reconstructed ancestral genome in .fa format has been outputted in `2_cactus/scratch` and should be moved/linked to `output/ancestral`.

The snakemake integrated workflow in `scripts/6_snpeff_gerp/snpeff/1_ancestral_as_reference` results in the SnpEff-annotated, ancestral 'reference' allele-converted .vcf file stored in `output/ancestral/` which is then also used for downstream analyses with GERP.

## GERP++

Next, the MAF files per scaffold are used as input for GERP, and using a snakemake workflow, GERP/RS-scores are calculated per scaffold in the `scripts/6_snpeff_gerp/gerp/1_snakemake_run_gerp` script, resulting in 30 .rates files (only looking at the 30 largest scaffolds). We convert the .rates files as well as the .vcf file to .bed format to use bedtools for overlapping the two types of datafiles: we only calculate genetic load based on variant sites. A snakemake workflow is used to calculate different types of load: expressed, masked, and total (codominant and additive) load for each GERP score category. Since we have the GERP scores per scaffold, an additional script is required to sum up all the scores across scaffolds

