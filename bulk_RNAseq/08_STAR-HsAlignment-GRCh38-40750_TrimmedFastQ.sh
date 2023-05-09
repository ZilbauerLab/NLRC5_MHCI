#!/bin/bash

REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/GRCh38.p13_ens99_STAR"
FASTQ="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/K.NAYAK-40750_AF"

for FILE in $FASTQ/*.fastq.gz; do

ID="${FILE##*/}"
NAME=${ID%.fastq.gz}

STAR --runThreadN 16 --genomeDir $REFERENCE --readFilesIn $FILE --readFilesCommand gunzip -c --outFileNamePrefix GRCh38_ALIGNED_BAMS/${NAME}. --outSAMtype BAM SortedByCoordinate

#echo $FILE
#echo $ID
#echo $NAME

done


