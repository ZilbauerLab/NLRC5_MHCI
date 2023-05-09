#!/bin/bash

# Sample libraries have been run across multiple lanes

for FILE in NGS-K.NAYAK-40706_RE/run1_fastqs/*.fastq.gz; do

ID="${FILE##*/}"
NAME=${ID%.fastq.gz}

zcat NGS-K.NAYAK-40706_RE/run1_fastqs/$ID NGS-K.NAYAK-40706_RE/run2_fastqs/$ID > NGS-K.NAYAK-40706_RE/MERGED-40706_RE/Merged_$NAME.fastq;
gzip NGS-K.NAYAK-40706_RE/MERGED-40706_RE/Merged_$NAME.fastq

done

