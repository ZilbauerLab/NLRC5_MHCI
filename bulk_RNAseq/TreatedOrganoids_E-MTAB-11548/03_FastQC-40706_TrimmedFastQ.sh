#!/bin/bash

# Submit this from /home/fp215/SOFTWARE/FastQC_v0.11.9

RE="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/K.NAYAK-40706_RE"
OUTPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQC/"

./fastqc $RE/*.fastq.gz --outdir $OUTPUT




