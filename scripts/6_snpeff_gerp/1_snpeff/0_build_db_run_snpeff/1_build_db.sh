#### In this script we build a lyrurus tetrix database using gene annotation and ref sequence
#### from https://pcingola.github.io/SnpEff/se_build_db/#add-a-genome-to-the-configuration-file

#step 1: configure genome - edit config file: lyrurus_tetrix.genome: blackgrouse

#step 2 option 3: build database from gff

#use different toolkits to extract cds, genes, protein

/vol/cluster-data/rchen/annotation/gff2fasta/gff3_to_fasta -g data/genomic/annotation/PO2979_Lyrurus_tetrix_black_grouse.annotation.gff \
    -f data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta -st cds -d complete -o data/genomic/intermediate/snpef/cds.fa 

/vol/cluster-data/rchen/annotation/gff2fasta/gff3_to_fasta -g data/genomic/annotation/PO2979_Lyrurus_tetrix_black_grouse.annotation.gff \
    -f data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta -st gene -d complete -o data/genomic/intermediate/snpef/genes.fa

#protein conversion with AGAT
agat_sp_extract_sequences.pl --gff data/genomic/annotation/PO2979_Lyrurus_tetrix_black_grouse.annotation.gff -f \
    data/genomic/refgenome/PO2979_Lyrurus_tetrix_black_grouse.RepeatMasked.fasta -p -o \
    data/genomic/intermediate/snpef/protein.fa

## build database, make sure to install SnpEff on your local device
java -jar src/snpEff.jar build -gff3 -v lyrurus_tetrix
