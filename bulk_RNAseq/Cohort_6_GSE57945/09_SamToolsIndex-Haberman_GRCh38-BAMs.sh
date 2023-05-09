#!/bin/bash

INPUT="/home/fp215/rds/rds-fp215-working/HABERMAN/ALIGNED_BAMS"

for FILE in $INPUT/*.Aligned.sortedByCoord.out.bam; do

samtools index $FILE

done

