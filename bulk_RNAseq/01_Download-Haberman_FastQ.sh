#!/bin/bash

for FILE in `cat GSE57945_sra_explorer_fastq_urls.txt`; do

wget $FILE;

done


