#!/bin/bash

INPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/HOWELL"
OUTPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/HOWELL"
REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/ADAPTERS"

for FILE in $INPUT/*_1.fastq.gz; do

NAME=$(basename "$FILE")
ID="${NAME%_1.fastq.gz}"

/home/fp215/SOFTWARE/BBMap_38.26/bbduk.sh in1=$INPUT/${ID}_1.fastq.gz in2=$INPUT/${ID}_2.fastq.gz out1=$OUTPUT/Trimmed_${ID}_1.fastq.gz out2=$OUTPUT/Trimmed_${ID}_2.fastq.gz ref=$REFERENCE/illumina.truseq.sequences.fa ktrim=r k=23 mink=11 hdist=1 bhist=$OUTPUT/LOGS/$ID.bhist.txt qhist=$OUTPUT/LOGS/$ID.qhist.txt gchist=$OUTPUT/LOGS/$ID.gchist.txt aqhist=$OUTPUT/LOGS/$ID.aqhist.txt lhist=$OUTPUT/LOGS/$ID.lhist.txt gcbins=auto tpe tbo -Xmx12g &> $OUTPUT/LOGS/$ID.txt

#echo $FILE
#echo $NAME
#echo ${ID}

done

