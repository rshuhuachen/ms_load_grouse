library(data.table, R.utils)
setwd("/vol/cluster-data/rchen/git/genetic_load_ltet/scripts/snpeff/1_ancestral_as_reference")
vcf <- fread("processed_data/snpeff_annotated/ltet_ann_snp_output.vcf.gz")
bed <- fread("results/ancestral/anc_allele_assignment.bed.gz")

str(vcf)
str(bed)

annotated <- left_join(vcf, bed, by = c())

rule annotate_vcf:
  input:
#   vcf = "processed_data/snpeff_annotated/ltet_ann_snp_output.vcf.gz",
    vcf = "processed_data/ltet_snps_noindel_minDP20_maxmis0.7_maxdp60_minQC30.recode.vcf.gz",
    bed = "results/ancestral/anc_allele_assignment.bed.gz"
  output:
    vcf = temp( "results/ancestral/ltet_filtered_ann_aa.vcf" )
  conda: "popgen_basics"
  shell:
    """
    zcat {input.vcf} | \
      vcf-annotate -a {input.bed} \
        -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
        -c CHROM,FROM,TO,INFO/AA > {output.vcf}
    """
