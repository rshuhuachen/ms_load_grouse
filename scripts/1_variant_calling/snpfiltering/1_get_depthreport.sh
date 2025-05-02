vcftools --gzvcf data/genomic/intermediate/rawSNPcalls.vcf.gz --site-mean-depth --out output/genotyping/DC_report
vcftools --gzvcf data/genomic/intermediate/rawSNPcalls.vcf.gz --depth --out output/genotyping/DC_report_id
