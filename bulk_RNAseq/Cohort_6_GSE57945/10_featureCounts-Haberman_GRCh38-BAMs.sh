#!/bin/bash

REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh38.99.gtf"
SUBREAD="/home/fp215/SOFTWARE/subread-2.0.3-source/bin"

$SUBREAD/featureCounts -t exon -g gene_id -C -a $REFERENCE -o Haberman_GSE57945_per-gene-counts.txt ALIGNED_BAMS/*.bam


