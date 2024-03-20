fai <- data.table::fread("scripts/snpeff/1_ancestral_as_reference/processed_data/tom-uni3242-mb-hirise-s8mkv__08-23-2022__hic_output.fasta.fai", fill=T)

fai$scaf <- gsub(";", ":", fai$scaf)
fai$scaf <- gsub("=", ".", fai$scaf)

write.csv(fai, file = "scripts/snpeff/1_ancestral_as_reference/processed_data/tom-uni3242-mb-hirise-s8mkv__08-23-2022__hic_output.csv", row.names=F, quote=F)
