#!/bin/bash

FASTA="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
GTF="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh37.87.gtf"
OUTPUT="/home/fp215/rds/hpc-work/REFERENCE/GRCh37.p13.release-107_STAR"

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir $OUTPUT --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 100
