# In this file, we will use the cactus-update-prepare function to create two scripts that will allow us to add two genomes to our dataset

# set scratch directory
CACTUS_SCRATCH=$(pwd)/scratch/

#enter container
apptainer shell --cleanenv \
  --fakeroot --overlay ${CACTUS_SCRATCH} \
  --bind ${CACTUS_SCRATCH}/tmp:/tmp,$(pwd) \
  --env PYTHONNOUSERSITE=1 \
  src/containers/cactus_v2.5.1.sif 

sh

#first add Lyrurus tetrix, branchlengths will be corrected in a later step

cactus-update-prepare \
  add branch \
  --parentGenome birdAnc334 \
  --childGenome Tympanuchus_cupido \
  data/genomic/intermediate/cactus/363-avian-reduced.hal \
  scripts/2_cactus/input_ltet.txt \
  --cactus-prepare-options \
  '--alignCores 4' \
  --topBranchLength 0.01 \ 
  --outDir scratch/tmp/steps-output \
  --jobStore scratch/tmp/js \
  --ancestorName AncX > scripts/2_cactus/3_cactus_update_lyrurus_steps.sh 

# then add Lagopus leucura

cactus-update-prepare \
  add branch \
  --parentGenome AncX \
  --childGenome Lyrurus_tetrix \
  data/genomic/intermediate/cactus/363-avian-reduced.hal \
  scripts/2_cactus/input_lleu.txt \
  --cactus-prepare-options \
  '--alignCores 4' \
  --topBranchLength 0.01 \
  --outDir scratch/tmp/steps-output \
  --jobStore scratch/tmp/js \
  --ancestorName AncY > scripts/2_cactus/4_cactus_update_lagopus_steps.sh 