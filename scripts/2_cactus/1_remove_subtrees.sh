### Remove subtrees that are not needed for ltet analysis

#set cactus scratch directory
CACTUS_SCRATCH=$(pwd)/scratch/

# enter the container
apptainer shell --cleanenv \
  --fakeroot --overlay ${CACTUS_SCRATCH} \
  --bind ${CACTUS_SCRATCH}/tmp:/tmp,$(pwd) \
  --env PYTHONNOUSERSITE=1 \
  docker:quay.io/comparative-genomics-toolkit/cactus:v2.5.1 

# get stats on the original 363-avian multi-alignment file
halStats data/genomic/intermediate/cactus/363-avian-2020.hal > output/cactus/stats_original_363_hal.txt

# copy the original file to then edit it
cp data/genomic/intermediate/cactus/363-avian-2020.hal data/genomic/intermediate/cactus/363-avian-reduced.hal

# remove subtrees to exclude neoaves
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc1
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc57 
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc69 
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc318 #starts with Heliornis_fulica
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc319 #starts with Psophia_crepitans
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc320 #starts with Charadrius_vociferus
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc321 #starts with Opisthocomus_hoazin
halRemoveSubtree data/genomic/intermediate/cactus/363-avian-reduced.hal birdAnc322 #stars with birdAnc57, so the big chunk of passerines but also a 

#some individual ancestral genomes left to exclude
halRemoveGenome data/genomic/intermediate/cactus//363-avian-reduced.hal birdAnc322
halRemoveGenome data/genomic/intermediate/cactus//363-avian-reduced.hal birdAnc1

# get stats of our subset of genomes

halStats data/genomic/intermediate/cactus/363-avian-reduced.hal > output/cactus/stats_reduced_363_hal.txt

#extract the reduced file 
halExtract data/genomic/intermediate/cactus/363-avian-reduced.hal data/genomic/intermediate/cactus/363-avian-reduced.hal 