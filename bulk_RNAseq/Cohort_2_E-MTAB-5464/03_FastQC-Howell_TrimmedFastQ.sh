#!/bin/bash

# Submit this from /home/fp215/SOFTWARE/FastQC_v0.11.9

HOWELL="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQS/HOWELL"
OUTPUT="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/TRIMMED_FASTQC/"


./fastqc $HOWELL/*.fastq.gz --outdir $OUTPUT




