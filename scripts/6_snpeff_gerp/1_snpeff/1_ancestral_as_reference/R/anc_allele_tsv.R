# the hal file will be output in 2_cactus/scratch and should be moved/linked to output/cactus
hal =  "output/ancestral/Anc66.hal"

load(file="data/genomic/raw/metadata/scaffold_names_dovetail.RData")
output_dir = "scripts/snpeff/1_ancestral_as_reference/processed_data/tmp/"

run_anc_allele <- function(hal, scaf, outdir, i){
    system(paste0('apptainer run --cleanenv --fakeroot --overlay scratch --bind ./scratch/tmp:/tmp,./scratch --env PYTHONNOUSERSITE=1 src/containers/cactus_v2.5.1.sif halSnps ', hal, ' Lyrurustetrix Anc66 --refSequence "', scaf, '" --tsv ', output_dir, 'tmp_scaf_', i))
}

for (i in 1:nrow(genome)){
  run_anc_allele(hal = hal,
  scaf = genome[i, 1],
  outdir = output_dir,
  i = i)
}

