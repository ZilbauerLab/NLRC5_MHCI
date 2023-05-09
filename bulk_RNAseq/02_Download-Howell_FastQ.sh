#!/bin/bash

cd /rds/project/fp215/rds-fp215-working/RNASEQ/NLRC5_PAPER/HOWELL

for FILE in `cat ../Howell_ftp.txt`; do

wget $FILE;

done
