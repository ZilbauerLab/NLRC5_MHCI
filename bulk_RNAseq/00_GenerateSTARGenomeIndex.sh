#!/bin/bash

FASTA="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh38.99.gtf"
OUTPUT="/home/fp215/rds/hpc-work/REFERENCE/GRCh38.p13_ens99_STAR"

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $OUTPUT --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 100
