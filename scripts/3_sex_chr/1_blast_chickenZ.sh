### first split lyrurus tetrix fasta in only the 30 largest scafs

# enter samtools 

/vol/apptainer/bin/apptainer shell --cleanenv \
  --fakeroot \
  --bind $(pwd) \
  --env PYTHONNOUSERSITE=1 src/containers/samtools_1.14--hb421002_0.sif

chr="ScEsiA3_18278;HRSCAF=21663 ScEsiA3_19038;HRSCAF=23783  ScEsiA3_18316;HRSCAF=21772 ScEsiA3_16870;HRSCAF=19426 ScEsiA3_18461;HRSCAF=22153 ScEsiA3_17641;HRSCAF=20530 ScEsiA3_15486;HRSCAF=17393  ScEsiA3_16766;HRSCAF=19082 ScEsiA3_17007;HRSCAF=19631  ScEsiA3_17655;HRSCAF=20552 ScEsiA3_16641;HRSCAF=18713  ScEsiA3_19552;HRSCAF=24415 ScEsiA3_12205;HRSCAF=13534  ScEsiA3_13761;HRSCAF=15525  ScEsiA3_669;HRSCAF=1222 ScEsiA3_292;HRSCAF=720 ScEsiA3_17616;HRSCAF=20466  ScEsiA3_16785;HRSCAF=19165  ScEsiA3_18290;HRSCAF=21698  ScEsiA3_18386;HRSCAF=21885  ScEsiA3_21975;HRSCAF=26924 ScEsiA3_21978;HRSCAF=26928 ScEsiA3_21977;HRSCAF=26927 ScEsiA3_21979;HRSCAF=26929 ScEsiA3_16771;HRSCAF=19097  ScEsiA3_21976;HRSCAF=26926  ScEsiA3_18752;HRSCAF=22883  ScEsiA3_12113;HRSCAF=13291  ScEsiA3_15456;HRSCAF=17358  ScEsiA3_18236;HRSCAF=21586"
gen="data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta"

for i in $chr ; do 
    samtools faidx $gen $i > output/genotyping/chicken_blast/scaf_$i.fa
done

### then blast against chicken Z chromosome (males are ZZ, females ZW)

makeblastdb -in data/genomic/intermediate/gallus_gallus_Zchr.fasta -dbtype nucl -parse_seqids -out chicken_Z_db

for i in output/genotyping/chicken_blast/*.fa; do 
    blastn -query "$i" -db chicken_Z_db -evalue 1e-5 -outfmt 6 -out output/genotyping/chicken_blast/out_scaf_$i.out
    done

