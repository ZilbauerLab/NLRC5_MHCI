#!/bin/bash

REFERENCE="/home/fp215/rds/hpc-work/REFERENCE/Homo_sapiens.GRCh38.99.gtf"
SUBREAD="/home/fp215/SOFTWARE/subread-2.0.3-source/bin"

$SUBREAD/featureCounts -t exon -g gene_id -p -C -T 4 -a $REFERENCE -o Howell_GRCh38_per-gene-counts.txt GRCh38_ALIGNED_BAMS/HOWELL/*.bam ;
$SUBREAD/featureCounts -t exon -g gene_id -C -T 4 -a $REFERENCE -o RE_GRCh38_per-gene-counts.txt GRCh38_ALIGNED_BAMS/RE/*.bam ;
$SUBREAD/featureCounts -t exon -g gene_id -C -T 4 -a $REFERENCE -o AF_GRCh38_per-gene-counts.txt GRCh38_ALIGNED_BAMS/AF/*.bam


