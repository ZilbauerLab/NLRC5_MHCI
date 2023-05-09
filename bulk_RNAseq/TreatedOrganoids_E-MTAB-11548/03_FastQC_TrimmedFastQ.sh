#!/bin/bash

# Submit this from /home/fp215/SOFTWARE/FastQC_v0.11.9

AF="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/K.NAYAK-40750_AF"
RE="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/K.NAYAK-40706_RE"
HOWELL="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/HOWELL"
OUTPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQC/"


./fastqc $AF/*.fastq.gz --outdir $OUTPUT ;
./fastqc $RE/*.fastq.gz --outdir $OUTPUT ;
./fastqc $HOWELL/*.fastq.gz --outdir $OUTPUT




