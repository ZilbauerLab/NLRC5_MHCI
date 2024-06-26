#!/bin/bash

REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/GRCh38.p13_ens99_STAR"
FASTQ="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/HOWELL"

for FILE in $FASTQ/*_1.fastq.gz; do

ID="${FILE##*/}"
NAME=${ID%_1.fastq.gz}

STAR --runThreadN 16 --genomeDir $REFERENCE --readFilesIn $FASTQ/${NAME}_1.fastq.gz $FASTQ/${NAME}_2.fastq.gz --readFilesCommand gunzip -c --outFileNamePrefix GRCh38_ALIGNED_BAMS/${NAME}. --outSAMtype BAM SortedByCoordinate

#echo $FILE
#echo $ID
#echo $NAME

done


