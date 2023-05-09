#!/bin/bash

INPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/NGS-K.NAYAK-40750_AF"
OUTPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/K.NAYAK-40750_AF"
REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/ADAPTERS"

for FILE in $INPUT/*.fastq.gz; do

NAME=$(basename "$FILE")
ID="${NAME%.fastq.gz}"

/home/fp215/SOFTWARE/BBMap_38.26/bbduk.sh in=$FILE out=$OUTPUT/Trimmed_$NAME ref=$REFERENCE/illumina.truseq.sequences.fa ktrim=r k=23 mink=11 bhist=$OUTPUT/LOGS/$ID.bhist.txt qhist=$OUTPUT/LOGS/$ID.qhist.txt gchist=$OUTPUT/LOGS/$ID.gchist.txt aqhist=$OUTPUT/LOGS/$ID.aqhist.txt lhist=$OUTPUT/LOGS/$ID.lhist.txt gcbins=auto tpe tbo -Xmx12g &> $OUTPUT/LOGS/$ID.txt

#echo $FILE
#echo $NAME
#echo $ID

done

