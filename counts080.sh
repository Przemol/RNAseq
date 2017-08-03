#!/bin/bash                               
#$ -l rmem=35G
#$ -l mem=35G                       
#$ -M n.bulgakova@sheffield.ac.uk
#$ -m bea

cd /mnt/fastdata/md1nbu/PRJNA369618/fastq/SRR5224080

module load apps/R/3.3.1               

R CMD BATCH count_reads_batch.R count_reads_batch.R.o$JOB_ID
