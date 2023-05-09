!/bin/bash

REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/GRCh38.p13_ens99_STAR"

for FILE in TRIMMED_DATA/*.fastq.gz; do

ID="${FILE##*/}"
NAME=${ID%.fastq.gz}

STAR --runThreadN 16 --genomeDir $REFERENCE --readFilesIn $FILE --readFilesCommand gunzip -c --outFileNamePrefix ALIGNED_BAMS/${NAME}. --outSAMtype BAM SortedByCoordinate

done
