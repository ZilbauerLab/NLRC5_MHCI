#!/bin/bash

INPUT="/home/fp215/rds/rds-fp215-working/HABERMAN/RAW_DATA"
OUTPUT="/home/fp215/rds/rds-fp215-working/HABERMAN/TRIMMED_DATA"
REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/ADAPTERS"

for FILE in $INPUT/*.fastq.gz; do

NAME=$(basename "$FILE")
ID="${NAME%.fastq.gz}"

/home/fp215/SOFTWARE/BBMap_38.26/bbduk.sh in=$FILE out=$OUTPUT/Trimmed_$NAME ref=$REFERENCE/illumina.truseq.sequences.fa ktrim=r k=23 mink=11 tpe tbo -Xmx12g &> $OUTPUT/LOGS/$ID.txt

#echo $FILE
#echo $NAME
#echo $ID

done


