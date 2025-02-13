## generated by : /home/cactus/cactus_env/bin/cactus-update-prepare add branch --parentGenome AncX --childGenome Lyrurus_tetrix data/363-avian-reduced.hal data/genomes/input_lleu.txt --cactus-prepare-options '--alignCores 4' --topBranchLength 0.01 --outDir scratch/tmp/steps-output --jobStore scratch/tmp/js --ancestorName AncY
## date : 2023-11-10 04:31:50.894280
## cactus commit : 4823037fe5038c8d2a9dd29a159e153e74a1aab5
## wrapping : /home/cactus/cactus_env/bin/cactus-prepare scratch/tmp/steps-output/seq_file.in --outDir scratch/tmp/steps-output --outSeqFile scratch/tmp/steps-output/seq_file.out --jobStore scratch/tmp/js --alignCores 4 --outHal scratch/tmp/steps-output/AncX.hal

## Preprocessor
cactus-preprocess scratch/tmp/js/0 scratch/tmp/steps-output/seq_file.in scratch/tmp/steps-output/seq_file.out --inputNames Lagopus_leucura --realTimeLogging --logInfo --retryCount 0 --maskMode none 

## Alignment

### Round 0
cactus-blast scratch/tmp/js/1 scratch/tmp/steps-output/seq_file.out scratch/tmp/steps-output/AncY.cigar --root AncY  
cactus-align scratch/tmp/js/2 scratch/tmp/steps-output/seq_file.out scratch/tmp/steps-output/AncY.cigar scratch/tmp/steps-output/AncY.hal --root AncY  --maxCores 4 
hal2fasta scratch/tmp/steps-output/AncY.hal AncY --hdf5InMemory > scratch/tmp/steps-output/AncY.fa

### Round 1
cactus-blast scratch/tmp/js/3 scratch/tmp/steps-output/seq_file.out scratch/tmp/steps-output/AncX.cigar --root AncX --includeRoot  
cactus-align scratch/tmp/js/4 scratch/tmp/steps-output/seq_file.out scratch/tmp/steps-output/AncX.cigar scratch/tmp/steps-output/AncX.hal --root AncX  --maxCores 4 --includeRoot 

## Alignment update
halAddToBranch data/genomic/intermediate/cactus/363-avian-reduced.hal scratch/tmp/steps-output/AncY.hal scratch/tmp/steps-output/AncX.hal AncX AncY Lyrurus_tetrix Lagopus_leucura 0.01 1.0 --hdf5InMemory 

## Alignment validation
halValidate --genome AncX data/genomic/intermediate/cactus/363-avian-reduced.hal --hdf5InMemory
halValidate --genome AncY data/genomic/intermediate/cactus/363-avian-reduced.hal --hdf5InMemory
halValidate --genome Lyrurus_tetrix data/genomic/intermediate/cactus/363-avian-reduced.hal --hdf5InMemory
halValidate --genome Lagopus_leucura data/genomic/intermediate/cactus/363-avian-reduced.hal --hdf5InMemory

