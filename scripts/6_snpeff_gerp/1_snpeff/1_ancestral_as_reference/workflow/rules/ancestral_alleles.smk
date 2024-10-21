
"""
from Kosmas https://github.com/k-hench/elephant_seals/blob/master/code/workflow/rules/ancestral_alleles.smk

snakemake --rerun-triggers mtime -n all_anc_allele

snakemake --jobs 100 \
  --latency-wait 30 --rerun-incomplete \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind /vol/cluster-data/rchen/git/genetic_load_ltet" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_aa.{name}.{jobid}.sh \
  -R anc_allele_tsv

"""

rule all_anc_allele:
    input: 
      vcf = "processed_data/snpeff_annotated/ltet_ann_snp_output.vcf.gz"

# skip instead do this with R script per scaffold in R/run_anc_allel_rsv.sh
#rule anc_allele_tsv:
#    input:
#      hal = "processed_data/Anc66.hal"
#    output:
#      tsv = "processed_data/ancestral_66_ltet_snps.tsv"
#    resources:
#      mem_mb=15360
#    container: c_cactus
#    shell:
#      """
#      halSnps {input.hal} Lyrurustetrix Anc66 --tsv {output.tsv}
#      """

# after R script, do awk FNR!=1 processed_data/tmp/tmp* > output/ancestral/ancestral_66_ltet_snps.tsv (and add headers manually)

rule pack_anc_allele_tsv:
    input:
      tsv = "output/ancestral/ancestral_66_ltet_snps.tsv"
    output:
      gz = "output/ancestral/ancestral_66_ltet_snps.tsv.gz"
    resources:
      mem_mb=15360
    shell:
      """
      gzip {input.tsv}
      """


rule anc_tree:
    input:
      hal = "output/ancestral/Anc66.hal"
    output:
      nwk = "output/ancestral/ltet_anc_tree.nwk"
    container: c_cactus
    shell:
      """
      halStats {input.hal} --tree  > {output.nwk}
      """

rule snp_pos_and_alleles:
    input:
      vcf = "data/genomic/intermediate/snpef/ltet_ann_snp_output.vcf"
    output:
      snps_tsv = "output/ancestral/ltet_snps_vcf.tsv.gz"
    shell:
      """
      grep -v "^##" {input.vcf} | cut -f 1,2,4,5 | gzip > {output.snps_tsv}
      """

rule determine_anc_ref:
    input:
      tsv = "output/ancestral/ancestral_66_ltet_snps.tsv.gz",
      snps_tsv = "output/ancestral/ltet_snps_vcf.tsv.gz"
    output:
      anc_tsv = "output/ancestral/anc_allele_assignment.tsv.gz",
      anc_bed = "output/ancestral/anc_allele_assignment.bed",
      tex_miss = "output/ancestral/ancestral_allele_mismatches.tex"
    log:
      "logs/r_ancestral_ref_proposal.log"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla scripts/7_deleterious_mutations/snpeff/1_ancestral_as_reference/R/ancient_alleles_assignment.R 2> {log} 1> {log}
      """

rule remove_string_ids:
  input:
    vcf = "data/genomic/intermediate/snpef/ltet_ann_snp_output.vcf.gz"
  output:
    vcf = "data/genomic/intermediate/snpef/ltet_ann_snp_output_cleanid.vcf.gz"
  shell:
    """
    zcat {input.vcf} | sed 's/\/vol\/cluster-data\/rchen\/wgr\/data\/processed\///g' | sed 's/.sorted.bam//g' | sed 's/;HRSCAF=/__HRSCAF_/' | gzip > {output.vcf}
    """

rule clean_scaf_name_bed:
  input:
    bed = "output/ancestral/anc_allele_assignment.bed"
  output:
    bed = "output/ancestral/anc_allele_assignment_cleanname.bed"
  shell:
    """
    cat {input.bed} | sed 's/;HRSCAF=/__HRSCAF_/' > {output.bed}
    """

rule pack_aa_bed:
    input:
      bed = "output/ancestral/anc_allele_assignment_cleanname.bed"
    output:
      gzbed = "output/ancestral/anc_allele_assignment_cleanname.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      bgzip {input.bed}
      tabix -s 1 -b 2 -e 3 {output.gzbed}
      """

rule annotate_vcf:
  input:
    vcf = "data/genomic/intermediate/snpef/ltet_ann_snp_output_cleanid.vcf.gz",
    bed = "output/ancestral/anc_allele_assignment_cleanname.bed.gz"
  output:
    vcf = temp( "output/ancestral/ltet_filtered_ann_aa.vcf" )
  conda: "popgen_basics"
  shell:
    """
    zcat {input.vcf} | \
      vcf-annotate -a {input.bed} \
        -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
        -c CHROM,FROM,TO,INFO/AA > {output.vcf}
    """

rule convert_vcf_alleles:
  input:
    vcf = "output/ancestral/ltet_filtered_ann_aa.vcf"
  output:
    gzvcf = "output/ancestral/ltet_filtered_ann_aa.vcf_wrongscaf.gz"
  resources:
      mem_mb=25600
  container: c_jvar
  shell:
    """
    java -jar /opt/jvarkit/dist/jvarkit.jar \
      vcffilterjdk \
      -f workflow/containers/script.js {input.vcf} | \
      bgzip > {output.gzvcf}
    """


rule revert_scaf_name:
  input:
    gzvcf = "output/ancestral/ltet_filtered_ann_aa_wrongscaf.vcf.gz"
  output:
    gzvcf = "output/ancestral/ltet_filtered_ann_aa.vcf.gz"
  shell:
    """
    zcat {input.gzvcf} | sed 's/__HRSCAF_/;HRSCAF=/' | bgzip > {output.gzvcf}
    """


rule index_aa_vcf:
    input:
      gzvcf = "output/ancestral/ltet_filtered_ann_aa.vcf.gz"
    output:
      idx = "output/ancestral/ltet_filtered_ann_aa.vcf.gz.tbi"
    conda: "popgen_basics"
    shell:
      """
      tabix -p vcf {input.gzvcf}
      """
