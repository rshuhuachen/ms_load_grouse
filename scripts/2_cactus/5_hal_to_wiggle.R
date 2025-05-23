
hal =  "data/genomic/intermediate/cactus/363-avian-reduced.hal"

library(dplyr); library(data.table)
scafs <- fread("data/genomic/refgenome/30_largest.scafs.tsv")

output_dir = "output/cactus/wiggle"
sif = "src/containers/cactus_v2.6.12.sif"
scratch = "scripts/2_cactus/scratch"
tmp_js = "scripts/2_cactus/scratch/tmp/js/wiggle"

hal_to_wiggle_per_scaf <- function(hal, scaf, outdir, scratch, i, sif, tmp_js){
    system(paste0('/vol/apptainer/bin/apptainer run --cleanenv --fakeroot --overlay ', scratch, ' --bind ', scratch, '/tmp:/tmp,', scratch, 
    ' --env PYTHONNOUSERSITE=1 ', sif, ' halAlignmentDepth ', hal, ' Lyrurus_tetrix --noAncestors --refSequence ', scaf, ' --outWiggle ', outdir, '/wig_', i, '.wig'))
}

for (i in 1:30){
  hal_to_wiggle_per_scaf(hal = hal,
  scaf = scafs$scaf[i],
  outdir = output_dir,
  scratch = scratch,
  sif = sif,
  i = i,
  tmp_js = tmp_js)
}

